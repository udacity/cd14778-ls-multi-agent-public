"""
Demo: Cancer Mutation Briefing

- Parallel retrieval (ClinVar, ClinicalTrials.gov v2, optional PubMed)
- Access control levels + basic PHI redaction (for logs/prompts)
- Deterministic knobs (model, temperature, seed)
- Trials filtering (recruiting/active), interventional-only by default
- Defensive parsing, graceful failure messages
- Provenance block + optional JSON sidecar output
- Simple in-memory caching for identical requests within a run

Output: crisp one-page briefing with explicit [ClinVar], [PMID:#######], [NCT:########] citations

Environment:
    .env file containing OPENAI_API_KEY
    Uses Vocareum base URL:
        https://openai.vocareum.com/v1
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import re
import time
import random
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple, Union

# --- Env & LLM client (LangChain) --------------------------------------------------------------
try:
    from dotenv import load_dotenv
except ImportError:
    print("Please install python-dotenv: pip install python-dotenv", file=sys.stderr)
    raise

try:
    from pydantic import BaseModel, Field
except ImportError:
    print("Please install pydantic: pip install pydantic", file=sys.stderr)
    raise

try:
    from langchain_openai import ChatOpenAI
    from langchain.schema.runnable import RunnableParallel, RunnableLambda, RunnableMap
    from langchain_core.messages import SystemMessage, HumanMessage
except ImportError:
    print("Please install langchain and langchain-openai: pip install langchain langchain-openai", file=sys.stderr)
    raise

# Course-provided read-only helpers
try:
    from ls_action_space.action_space import query_clinvar, query_pubmed, query_clinicaltrials
except Exception as e:
    print(
        "Could not import ls_action_space.action_space. "
        "Ensure the course-provided package is installed and importable.",
        file=sys.stderr,
    )
    raise

load_dotenv("../.env")

# --- Access control & redaction ("Good Stuff" from original) -----------------------------------
class PrivacyLevel:
    PUBLIC = "public"
    CUSTOMER = "customer"
    AGENT = "agent"
    FINANCIAL = "financial"
    ADMIN = "admin"

class AccessControl:
    @staticmethod
    def can_access(requester_level: str, data_level: str) -> bool:
        if requester_level == PrivacyLevel.ADMIN:
            return True
        if requester_level == data_level:
            return True
        if requester_level == PrivacyLevel.FINANCIAL and data_level == PrivacyLevel.AGENT:
            return True
        if requester_level == PrivacyLevel.AGENT and data_level == PrivacyLevel.CUSTOMER:
            return True
        if data_level == PrivacyLevel.PUBLIC:
            return True
        return False

class Redactor:
    """Very lightweight PHI redactor for logs/prompts."""
    EMAIL_RE = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
    PHONE_RE = re.compile(r"\b(?:\+?\d{1,3}[\s-]?)?(?:\(?\d{3}\)?[\s-]?)?\d{3}[\s-]?\d{4}\b")
    MRN_RE = re.compile(r"\bMRN[:\s]*\d{5,}\b", re.IGNORECASE)

    @classmethod
    def redact(cls, text: str) -> str:
        text = cls.EMAIL_RE.sub("[REDACTED_EMAIL]", text)
        text = cls.PHONE_RE.sub("[REDACTED_PHONE]", text)
        text = cls.MRN_RE.sub("MRN:[REDACTED]", text)
        return text

def pretty_json(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2, default=str)

def compact_json(obj: Any, limit: int = 2000) -> str:
    s = json.dumps(obj, ensure_ascii=False, separators=(",", ":"), default=str)
    return (s[: limit - 3] + "...") if len(s) > limit else s

# --- Determinism knobs --------------------------------------------------------------------------
def seed_everything(seed: Optional[int]):
    if seed is None:
        return
    random.seed(seed)

# --- Variant utilities --------------------------------------------------------------------------
def normalize_variant(user_str: str) -> Dict[str, Any]:
    """
    Best-effort normalizer to help queries:
      - Detect rsID pattern
      - Extract gene symbol and protein/codon if present (very light heuristic)
    """
    rsid = None
    gene = None
    protein = None

    m = re.search(r"(rs\d+)", user_str, re.IGNORECASE)
    if m:
        rsid = m.group(1)

    # Gene symbol heuristic: uppercase letters/numbers, 2-8 chars, at start
    gm = re.match(r"\s*([A-Z0-9]{2,8})\b", user_str)
    if gm:
        gene = gm.group(1)

    # Simple protein change like V600E, L858R, R132H, etc.
    pm = re.search(r"\b([A-Z]\d{2,5}[A-Z])\b", user_str)
    if pm:
        protein = pm.group(1)

    return {"raw": user_str, "rsid": rsid, "gene": gene, "protein": protein}

# --- Caching wrappers around read-only helpers --------------------------------------------------
# (Speeds up repeated demo runs within a single process)
@lru_cache(maxsize=128)
def cached_clinvar(variant: str) -> Dict[str, Any]:
    return query_clinvar(variant)

@lru_cache(maxsize=128)
def cached_pubmed(query: str, k: int) -> List[Dict[str, Any]]:
    return query_pubmed(query, max_results=k, include_mesh=True, include_citations=False)

@lru_cache(maxsize=128)
def cached_clinicaltrials(expr: str, max_results: int, fields_key: str) -> Dict[str, Any]:
    fields = json.loads(fields_key)
    return query_clinicaltrials(expr=expr, max_results=max_results, fields=fields)

# --- Tool I/O schemas (Pydantic) ----------------------------------------------------------------
class ClinVarInput(BaseModel):
    variant: str = Field(..., description="Variant string (rsID or HGVS).")

class TrialsInput(BaseModel):
    variant: str = Field(..., description="Variant string to search.")
    disease_hint: Optional[str] = Field(None, description="Optional disease or tumor context.")
    max_n: int = Field(30, ge=1, le=200, description="Max studies to retrieve.")
    interventional_only: bool = Field(True, description="Keep only interventional trials.")
    active_only: bool = Field(True, description="Keep only recruiting/active/enrolling.")

class PubMedInput(BaseModel):
    variant: str = Field(..., description="Variant string.")
    disease_hint: Optional[str] = Field(None, description="Optional disease context.")
    k: int = Field(3, ge=1, le=10, description="Number of articles to fetch.")

# --- Tools as Python callables ------------------------------------------------------------------
RECRUITING_KEEP = {"Recruiting", "Active, not recruiting", "Enrolling by invitation"}

TRIAL_FIELDS_DEFAULT = [
    "NCTId","BriefTitle","OverallStatus","StudyType","Condition","InterventionName",
    "PrimaryOutcomeMeasure","EnrollmentCount","StartDate","CompletionDate",
    "Phase","LocationCountry","LeadSponsorName"
]

def trials_expr(variant: str, disease_hint: Optional[str]) -> str:
    if disease_hint:
        return f'("{variant}") OR ({disease_hint})'
    return f'"{variant}"'

def _ensure_list(v: Any) -> List[Any]:
    if v is None:
        return []
    return v if isinstance(v, list) else [v]

def dedupe_by_key(items: List[Dict[str, Any]], key: str) -> List[Dict[str, Any]]:
    seen = set()
    out = []
    for it in items:
        k = it.get(key)
        if not k or k in seen:
            continue
        seen.add(k)
        out.append(it)
    return out

def clinvar_tool(args: ClinVarInput) -> Dict[str, Any]:
    try:
        data = cached_clinvar(args.variant)
        return data or {}
    except Exception as e:
        return {"error": f"ClinVar lookup failed: {e}", "query": args.variant}

def trials_tool(args: TrialsInput) -> Dict[str, Any]:
    fields_key = json.dumps(TRIAL_FIELDS_DEFAULT, ensure_ascii=False)
    expr = trials_expr(args.variant, args.disease_hint)
    try:
        out = cached_clinicaltrials(expr=expr, max_results=args.max_n, fields_key=fields_key)
        studies = out.get("studies", []) or []

        # Filter by interventional?
        if args.interventional_only:
            studies = [s for s in studies if (s.get("StudyType") or "").lower().startswith("interventional")]

        # Filter by active status?
        if args.active_only:
            studies = [s for s in studies if s.get("OverallStatus") in RECRUITING_KEEP]

        # Normalize types
        normed = []
        for s in studies:
            ns = dict(s)
            for lk in ("Condition", "InterventionName", "PrimaryOutcomeMeasure", "LocationCountry"):
                ns[lk] = _ensure_list(ns.get(lk))
            # Cast enrollment (best-effort)
            try:
                if ns.get("EnrollmentCount") is not None:
                    ns["EnrollmentCount"] = int(ns["EnrollmentCount"])
            except Exception:
                pass
            normed.append(ns)

        normed = dedupe_by_key(normed, "NCTId")
        return {"count": len(normed), "studies": normed, "expr": expr}
    except Exception as e:
        return {"error": f"ClinicalTrials query failed: {e}", "expr": expr}

def pubmed_tool(args: PubMedInput) -> List[Dict[str, Any]]:
    q = args.variant + (f" AND {args.disease_hint}" if args.disease_hint else "")
    try:
        # Light sanity: require variant or gene token presence in title/abstract after fetching
        raw = cached_pubmed(q, args.k) or []
        # No heavy filtering—just ensure relevance by token presence if variant is not rsID-only
        norm = normalize_variant(args.variant)
        tokens = [t for t in [norm["gene"], norm["protein"], norm["rsid"], args.variant] if t]
        tokens_lower = [t.lower() for t in tokens]
        filtered = []
        for rec in raw:
            text = ((rec.get("title") or "") + " " + (rec.get("abstract") or "")).lower()
            if not tokens_lower or any(t in text for t in tokens_lower):
                filtered.append(rec)
        return filtered[: args.k]
    except Exception as e:
        return [{"error": f"PubMed query failed: {e}", "query": q}]

# --- LLM prompt & synthesis --------------------------------------------------------------------
SYSTEM = (
    "You draft concise, evidence-aware clinical briefings for experts. "
    "Be precise, cautious, and cite which source each claim comes from "
    "using [ClinVar], [PMID:########], [NCT:########]. Do not provide medical advice."
)

def render_prompt(inputs: Dict[str, Any]) -> str:
    """Build a strict, sectioned prompt. Redact any accidental PHI in logs."""
    clinvar = inputs.get("clinvar") or {}
    trials = inputs.get("trials") or {}
    pubmed = inputs.get("pubmed") or []
    variant = inputs["variant"]
    disease_hint = inputs.get("disease_hint")

    # Redact only for safety in prompt context (usually not needed for these sources)
    clinvar_s = Redactor.redact(pretty_json(clinvar))
    trials_s  = Redactor.redact(pretty_json(trials))
    pubmed_s  = Redactor.redact(pretty_json(pubmed[:5]))  # cap display

    return f"""
Variant: {variant}
Disease hint: {disease_hint or "None"}

ClinVar (context):
{clinvar_s}

Trials (context):
{trials_s}

PubMed (context):
{pubmed_s}

Write a one-page briefing with EXACTLY these sections:

1) Summary
   - 2–3 sentences; include gene if present; neutral tone.

2) Clinical significance
   - ClinVar verdict and review_status; associated conditions if any.
   - Note conflicts/uncertainty.
   - Cite as [ClinVar].

3) Prevalence / allele frequency
   - If rsID frequencies exist, report allele -> frequency succinctly.
   - Otherwise state 'Not available from ClinVar rsID frequency aggregation'.
   - Cite as [ClinVar].

4) Evidence highlights
   - 1–3 bullets.
   - Each bullet: 1-sentence takeaway + [PMID:########].
   - Prefer recency and relevance to disease context.

5) Trial opportunities
   - Up to 3 entries.
   - Each: NCTId, status, phase, one-line rationale.
   - Cite as [NCT:########].

6) Caveats & next steps
   - Data limitations (sample size, review conflicts).
   - Tumor/tissue specificity considerations.
   - Sensible next steps (confirmatory testing, pathway-level search).
   
Formatting rules:
- Use short paragraphs and bullets; be crisp.
- Use only bracketed citations: [ClinVar], [PMID:########], [NCT:########].
- Do NOT fabricate PMIDs or NCT IDs—only cite those present in the context above.
- Avoid clinical directives (no dosing, no patient-specific recommendations).
""".strip()

# --- Provenance / audit -------------------------------------------------------------------------
def provenance_block(inputs: Dict[str, Any]) -> str:
    clinvar = inputs.get("clinvar") or {}
    trials = inputs.get("trials") or {}
    pubmed = inputs.get("pubmed") or []
    lines = [
        "---- Provenance --------------------------------------------------",
        f"Variant: {inputs.get('variant')}",
        f"Disease hint: {inputs.get('disease_hint') or 'None'}",
        f"ClinVar: {'OK' if 'error' not in clinvar else 'ERROR'}",
        f"Trials: count={len(trials.get('studies', [])) if isinstance(trials, dict) else 0}, expr={trials.get('expr') if isinstance(trials, dict) else None}",
        f"PubMed: count={len(pubmed) if isinstance(pubmed, list) else 0}",
        "-----------------------------------------------------------------",
    ]
    return "\n".join(lines)

# --- LangChain graph: parallel retrieval -> join -> synthesize ---------------------------------
def build_graph(llm: ChatOpenAI, include_pubmed: bool) -> Tuple[Any, Any]:
    # fan-out: run tools in parallel
    fan = RunnableParallel(
        clinvar=RunnableLambda(lambda x: clinvar_tool(ClinVarInput(variant=x["variant"]))),
        trials=RunnableLambda(
            lambda x: trials_tool(
                TrialsInput(
                    variant=x["variant"],
                    disease_hint=x.get("disease_hint"),
                    max_n=x.get("max_trials", 30),
                    interventional_only=x.get("interventional_only", True),
                    active_only=x.get("active_only", True),
                )
            )
        ),
        pubmed=RunnableLambda(
            lambda x: pubmed_tool(
                PubMedInput(
                    variant=x["variant"],
                    disease_hint=x.get("disease_hint"),
                    k=x.get("k_pubmed", 3),
                )
            )
        ) if include_pubmed else RunnableLambda(lambda _: [])
    )

    # join inputs and fan results
    join = RunnableMap({
        "variant": lambda x: x["args"]["variant"],
        "disease_hint": lambda x: x["args"].get("disease_hint"),
        "clinvar": lambda x: x["fan"]["clinvar"],
        "trials": lambda x: x["fan"]["trials"],
        "pubmed": lambda x: x["fan"]["pubmed"],
        # pass-through for later reporting
        "max_trials": lambda x: x["args"].get("max_trials", 30),
        "interventional_only": lambda x: x["args"].get("interventional_only", True),
        "active_only": lambda x: x["args"].get("active_only", True),
        "k_pubmed": lambda x: x["args"].get("k_pubmed", 3),
    })

    # compose pipeline
    graph = (
            RunnableLambda(lambda x: x)
            | RunnableMap({"args": lambda x: x, "fan": fan})
            | join
            | RunnableLambda(render_prompt)
            | RunnableLambda(lambda prompt: [SystemMessage(content=SYSTEM), HumanMessage(content=prompt)])
            | llm
    )

    return graph, join  # return join for provenance inputs if needed

# --- CLI / Runner --------------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Cancer Mutation Briefing — LangChain Demo (Life Sciences)"
    )
    # Make variant/disease OPTIONAL with good defaults
    p.add_argument("--variant", default="BRAF V600E",
                   help='Variant string (default: "BRAF V600E")')
    p.add_argument("--disease", default="melanoma",
                   help='Disease hint (default: "melanoma")')

    p.add_argument("--include-pubmed", action="store_true",
                   help="Include PubMed retrieval (default: off unless no args)")
    p.add_argument("--model", default="gpt-4o-mini")
    p.add_argument("--temperature", type=float, default=0.2)
    p.add_argument("--max-trials", type=int, default=30)
    p.add_argument("--k-pubmed", type=int, default=3)
    p.add_argument("--interventional-only", action="store_true",
                   help="Keep only interventional trials")
    p.add_argument("--active-only", action="store_true",
                   help="Keep only recruiting/active/enrolling")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--emit-json", action="store_true")

    # If user provided no CLI args, we’ll flip on a couple flags for a richer demo
    if argv is None and len(sys.argv) == 1:
        # No flags passed → run a friendly default example
        ns = p.parse_args([])
        ns.include_pubmed = True            # show PubMed by default in the zero-arg path
        ns.interventional_only = True
        ns.active_only = True
        ns._zero_arg_banner = True          # internal flag just to print a banner
        return ns

    ns = p.parse_args(argv)
    ns._zero_arg_banner = False
    return ns


def main(argv=None) -> int:
    args = parse_args(argv)

    if getattr(args, "_zero_arg_banner", False):
        print("No CLI arguments detected — running the default example:\n"
              f'  variant="{args.variant}", disease="{args.disease}", '
              f"include_pubmed={args.include_pubmed}, interventional_only={args.interventional_only}, "
              f"active_only={args.active_only}\n"
              "Tip: run `--help` to see all options.\n")
              
    # Determinism across random ops (not strictly needed here but useful for grading)
    seed_everything(args.seed)

    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        print("OPENAI_API_KEY not found. Please set it in .env or environment.", file=sys.stderr)
        return 2

    # Build LLM with Vocareum endpoint
    llm = ChatOpenAI(
        model=args.model,
        base_url="https://openai.vocareum.com/v1",
        temperature=args.temperature,
    )

    # Build graph
    graph, _ = build_graph(llm=llm, include_pubmed=bool(args.include_pubmed))

    # Inputs
    user_inputs = {
        "variant": args.variant,
        "disease_hint": args.disease,
        "max_trials": args.max_trials,
        "k_pubmed": args.k_pubmed,
        "interventional_only": bool(args.interventional_only),
        "active_only": bool(args.active_only),
    }

    # Invoke graph
    t0 = time.time()
    try:
        completion = graph.invoke(user_inputs)
        # LangChain ChatOpenAI returns an AIMessage; .content holds text
        briefing_text = getattr(completion, "content", str(completion))
    except Exception as e:
        print(f"Error during synthesis: {e}", file=sys.stderr)
        return 1
    dt = time.time() - t0

    # For provenance, recompute tool outputs once (fast via cache) to print counts/expr
    # We directly call the tools (which hit lru_cache) for a clean provenance summary.
    clinvar_out = clinvar_tool(ClinVarInput(variant=args.variant))
    trials_out  = trials_tool(
        TrialsInput(
            variant=args.variant,
            disease_hint=args.disease,
            max_n=args.max_trials,
            interventional_only=bool(args.interventional_only),
            active_only=bool(args.active_only),
        )
    )
    pubmed_out = pubmed_tool(
        PubMedInput(variant=args.variant, disease_hint=args.disease, k=args.k_pubmed)
    ) if args.include_pubmed else []

    # Assemble final printed output
    print(f"=== Mutation Briefing: {args.variant} ===")
    if args.disease:
        print(f"(disease hint: {args.disease})")
    print()
    print(briefing_text.strip())
    print()
    print(provenance_block({
        "variant": args.variant,
        "disease_hint": args.disease,
        "clinvar": clinvar_out,
        "trials": trials_out,
        "pubmed": pubmed_out,
    }))
    print(f"Generated in {dt:.2f}s")
    print("\n[Note] This content is informational and not medical advice.")

    # Optional JSON sidecar (useful for grading / downstream UI)
    if args.emit_json:
        sidecar = {
            "inputs": {
                "variant": args.variant,
                "disease_hint": args.disease,
                "max_trials": args.max_trials,
                "k_pubmed": args.k_pubmed,
                "interventional_only": bool(args.interventional_only),
                "active_only": bool(args.active_only),
                "include_pubmed": bool(args.include_pubmed),
                "model": args.model,
                "temperature": args.temperature,
                "seed": args.seed,
            },
            "normalized_variant": normalize_variant(args.variant),
            "clinvar": clinvar_out,
            "trials": trials_out,
            "pubmed": pubmed_out if args.include_pubmed else [],
            "provenance": {
                "notes": "All retrievals are read-only wrappers around public APIs provided by the course.",
                "filters": {
                    "trials": {
                        "interventional_only": bool(args.interventional_only),
                        "active_only": bool(args.active_only),
                        "status_kept": sorted(list(RECRUITING_KEEP)),
                    }
                }
            }
        }
        path = "mutation_briefing_sidecar.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(sidecar, f, ensure_ascii=False, indent=2)
        print(f"\n(JSON sidecar written to ./{path})")

    return 0

if __name__ == "__main__":
    main()

