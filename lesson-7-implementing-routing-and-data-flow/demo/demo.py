#!/usr/bin/env python3
"""
Pharmacovigilance Inbox — end-to-end demo (single file)

- Single intake stream (stdin or JSONL) → LLM classifier → priority queue → specialist agents
- Optional read-only enrichments (PubMed/ClinicalTrials/ClinVar/FAERS) if available
- Minimal, dependency-light, yet realistic structure

Usage examples:
  # Read reports from JSONL (each line: {"id": "...", "text": "...", "metadata": {...}})
  python pharmacovigilance_inbox.py --input reports.jsonl --enrich

  # Interactive stdin mode (type/paste a report; Ctrl-D to end)
  python pharmacovigilance_inbox.py --stdin

  # If you run with *no arguments*, an example auto-runs:
  python pharmacovigilance_inbox.py
"""

from __future__ import annotations

import argparse
import heapq
import json
import os
import sys
import time
from dataclasses import dataclass
from typing import Any, Dict, Generator, Iterable, List, Literal, Optional, Tuple

# --- LLM client setup (Vocareum) ------------------------------------------------
from dotenv import load_dotenv  # pip install python-dotenv
from openai import OpenAI       # pip install openai

load_dotenv("../.env")  # expects OPENAI_API_KEY
client = OpenAI(base_url="https://openai.vocareum.com/v1")

# --- Optional read-only life-science tools --------------------------------------
HAS_LS_TOOLS = True
try:
    from ls_action_space.action_space import (
        query_pubmed,
        query_clinicaltrials,
        query_clinvar,
        get_faers_contingency_counts,
    )
except Exception:
    HAS_LS_TOOLS = False


# --- Data structures -------------------------------------------------------------
Label = Literal["adverse_event", "manufacturing_issue", "general_feedback"]
Severity = Literal["fatal", "serious", "non_serious", "unknown"]

@dataclass
class Report:
    id: str
    text: str
    metadata: Dict[str, Any]

    @staticmethod
    def from_dict(obj: Dict[str, Any]) -> "Report":
        return Report(
            id=str(obj.get("id") or obj.get("report_id") or obj.get("uuid") or ""),
            text=str(obj.get("text") or obj.get("body") or obj.get("content") or ""),
            metadata=dict(obj.get("metadata") or {})
        )


# --- LLM classification ----------------------------------------------------------
SYSTEM_PROMPT = (
    "You are a pharmacovigilance triage assistant. "
    "Read the input report text and return a STRICT JSON object with keys:\n"
    '{\n'
    '  "label": one of ["adverse_event","manufacturing_issue","general_feedback"],\n'
    '  "severity": one of ["fatal","serious","non_serious","unknown"],\n'
    '  "entities": {\n'
    '      "product": string[] (suspected or mentioned products/drugs/devices),\n'
    '      "events": string[] (clinical events; prefer MedDRA PT if clear),\n'
    '      "lot": string[] (batch/lot identifiers if present)\n'
    "  }\n"
    "}\n"
    "Output ONLY JSON. No extra text."
)

DEFAULT_CLASSIFICATION: Dict[str, Any] = {
    "label": "general_feedback",
    "severity": "unknown",
    "entities": {"product": [], "events": [], "lot": []},
}


def _safe_json_parse(text: str) -> Optional[Dict[str, Any]]:
    try:
        return json.loads(text)
    except Exception:
        start = text.find("{")
        end = text.rfind("}")
        if start != -1 and end != -1 and end > start:
            try:
                return json.loads(text[start : end + 1])
            except Exception:
                return None
    return None


def classify_report_llm(text: str, model: str = "gpt-4o-mini", max_retries: int = 3, timeout_s: float = 20.0) -> Dict[str, Any]:
    """
    Ask the LLM to classify report text and extract entities.
    Returns a dict with keys: label, severity, entities (never None).
    """
    deadline = time.time() + timeout_s
    last_err = None

    for attempt in range(1, max_retries + 1):
        try:
            resp = client.chat.completions.create(
                model=model,
                temperature=0,
                messages=[
                    {"role": "system", "content": SYSTEM_PROMPT},
                    {"role": "user", "content": text},
                ],
            )
            content = (resp.choices[0].message.content or "").strip()
            data = _safe_json_parse(content) or {}
            label = str(data.get("label", DEFAULT_CLASSIFICATION["label"]))
            severity = str(data.get("severity", DEFAULT_CLASSIFICATION["severity"]))
            entities = data.get("entities") or {}
            if not isinstance(entities, dict):
                entities = {}
            entities.setdefault("product", [])
            entities.setdefault("events", [])
            entities.setdefault("lot", [])
            return {"label": label, "severity": severity, "entities": entities}
        except Exception as e:
            last_err = e
            if time.time() > deadline:
                break
            time.sleep(min(0.4 * attempt, 2.0))

    sys.stderr.write(f"[WARN] Falling back to default classification: {last_err}\n")
    return DEFAULT_CLASSIFICATION.copy()


# --- Specialist agents -----------------------------------------------------------
class Agent:
    def handle(self, report: Report, info: Dict[str, Any], *, enrich: bool = False) -> str:
        raise NotImplementedError


class AdverseEventAgent(Agent):
    def handle(self, report: Report, info: Dict[str, Any], *, enrich: bool = False) -> str:
        entities = info.get("entities", {})
        products = [p for p in entities.get("product", []) if isinstance(p, str)]
        events = [e for e in entities.get("events", []) if isinstance(e, str)]

        enrich_notes: List[str] = []
        if enrich and HAS_LS_TOOLS:
            try:
                if products and events:
                    counts = get_faers_contingency_counts(products[0], events[0])
                    enrich_notes.append(
                        f"FAERS 2x2 (a,b,c,d)=({counts['a']},{counts['b']},{counts['c']},{counts['d']})"
                    )
                q = " AND ".join(
                    [f'"{products[0]}"'] if products else [] +
                    [events[0]] if events else []
                ) or None
                if q:
                    pubs = query_pubmed(q, max_results=3, include_mesh=True)
                    titles = [p.get("title") for p in pubs if p.get("title")]
                    if titles:
                        enrich_notes.append(f"PubMed hits: {len(titles)} (e.g., {titles[0][:120]}...)")
            except Exception as e:
                enrich_notes.append(f"Enrichment skipped: {e!s}")

        summary = {
            "agent": "AdverseEventAgent",
            "report_id": report.id,
            "label": info.get("label"),
            "severity": info.get("severity"),
            "products": products,
            "events": events,
            "enrichment": enrich_notes or None,
        }
        return json.dumps(summary, ensure_ascii=False)


class ManufacturingIssueAgent(Agent):
    def handle(self, report: Report, info: Dict[str, Any], *, enrich: bool = False) -> str:
        entities = info.get("entities", {})
        lots = [x for x in entities.get("lot", []) if isinstance(x, str)]

        enrich_notes: List[str] = []
        if enrich and HAS_LS_TOOLS:
            try:
                products = [p for p in entities.get("product", []) if isinstance(p, str)]
                if products:
                    ct = query_clinicaltrials(products[0], max_results=5)
                    enrich_notes.append(f"ClinicalTrials.gov studies ~ {ct.get('count')}")
            except Exception as e:
                enrich_notes.append(f"Enrichment skipped: {e!s}")

        summary = {
            "agent": "ManufacturingIssueAgent",
            "report_id": report.id,
            "label": info.get("label"),
            "severity": info.get("severity"),
            "lots": lots,
            "enrichment": enrich_notes or None,
        }
        return json.dumps(summary, ensure_ascii=False)


class GeneralFeedbackAgent(Agent):
    def handle(self, report: Report, info: Dict[str, Any], *, enrich: bool = False) -> str:
        summary = {
            "agent": "GeneralFeedbackAgent",
            "report_id": report.id,
            "label": info.get("label"),
            "severity": info.get("severity"),
            "note": "Routed to CX/Insights queue",
        }
        return json.dumps(summary, ensure_ascii=False)


# --- Router with priority queue --------------------------------------------------
class Router:
    """
    Content-based routing (by LLM label) + priority-based dispatch (by severity).
    """
    def __init__(self):
        self.ae_agent = AdverseEventAgent()
        self.mfg_agent = ManufacturingIssueAgent()
        self.fb_agent = GeneralFeedbackAgent()
        self._pq: List[Tuple[int, float, Report, Dict[str, Any]]] = []

    @staticmethod
    def _priority_from_severity(severity: Severity) -> int:
        return {"fatal": 0, "serious": 1, "non_serious": 2, "unknown": 3}.get(severity, 3)

    def ingest(self, report: Report, info: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        info = info or classify_report_llm(report.text)
        priority = self._priority_from_severity(info.get("severity", "unknown"))  # type: ignore[arg-type]
        heapq.heappush(self._pq, (priority, time.monotonic(), report, info))
        return info

    def dispatch_next(self, *, enrich: bool = False) -> Optional[str]:
        if not self._pq:
            return None
        _, _, report, info = heapq.heappop(self._pq)
        label = info.get("label", "general_feedback")
        if label == "adverse_event":
            return self.ae_agent.handle(report, info, enrich=enrich)
        elif label == "manufacturing_issue":
            return self.mfg_agent.handle(report, info, enrich=enrich)
        else:
            return self.fb_agent.handle(report, info, enrich=enrich)

    def empty(self) -> bool:
        return not self._pq


# --- Input streaming (no hardcoded data) ----------------------------------------
def read_jsonl(path: str) -> Generator[Report, None, None]:
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            rpt = Report.from_dict(obj)
            if rpt.id and rpt.text:
                yield rpt


def read_stdin() -> Generator[Report, None, None]:
    """
    Interactive mode: user types/pastes report text one per line.
    IDs are generated incrementally (no static content).
    """
    print("# Enter one report per line (Ctrl-D to end):", file=sys.stderr)
    for idx, line in enumerate(sys.stdin, start=1):
        text = line.strip()
        if not text:
            continue
        yield Report(id=f"stdin-{idx}", text=text, metadata={})


# --- Auto example (runs when script has no argv) --------------------------------
EXAMPLE_GEN_SYSTEM = (
    "You are generating synthetic, anonymized pharmacovigilance-style reports "
    "for a demo. Create diverse short texts with invented generic product names "
    "(e.g., Drug-A, Device-B), invented lots (e.g., LOT-1234), and varied concerns. "
    "Do not use real brand names or identifiable persons. Return STRICT JSONL: "
    'each line is an object with keys: id (string), text (string), metadata (object).'
)

def generate_reports_via_llm(n: int = 3, model: str = "gpt-4o-mini") -> List[Report]:
    """
    Uses the LLM to synthesize example reports at runtime (no hardcoded data).
    Falls back to a minimal synthetic generator if parsing fails.
    """
    try:
        user_prompt = (
            f"Produce {n} demo reports. Vary across adverse events, manufacturing issues, "
            "and general feedback. Keep each text under 80 words."
        )
        resp = client.chat.completions.create(
            model=model,
            temperature=0,
            messages=[
                {"role": "system", "content": EXAMPLE_GEN_SYSTEM},
                {"role": "user", "content": user_prompt},
            ],
        )
        content = (resp.choices[0].message.content or "").strip()
        reports: List[Report] = []
        # Accept either raw JSONL or a JSON array (be forgiving)
        if "\n" in content and content.strip().startswith("{"):
            # Likely JSONL
            for i, line in enumerate(content.splitlines(), start=1):
                obj = _safe_json_parse(line)
                if not obj:
                    continue
                reports.append(Report.from_dict(obj))
        else:
            data = _safe_json_parse(content)
            if isinstance(data, list):
                for obj in data:
                    reports.append(Report.from_dict(obj))
        # Final sanitation
        reports = [r for r in reports if r.id and r.text]
        if reports:
            return reports
    except Exception as e:
        sys.stderr.write(f"[WARN] Example LLM generation failed, falling back: {e}\n")

    # Fallback: generate generic, non-identifying placeholder reports (algorithmic, not hardcoded corpus)
    fallback: List[Report] = []
    for i in range(1, n + 1):
        kind = ["adverse event", "manufacturing issue", "general feedback"][((i - 1) % 3)]
        txt = (
            f"Customer note {i}: possible {kind} after use of Product-{i}. "
            f"Mentions LOT-{1000+i} and follow-up planned. Severity not stated."
        )
        fallback.append(Report(id=f"demo-{i}", text=txt, metadata={"source": "auto"}))
    return fallback


def run_auto_example(enrich: bool = False, model: str = "gpt-4o-mini") -> None:
    """
    Auto-runs when no CLI args are provided. Generates synthetic example reports at runtime,
    classifies them, routes via priority queue, and prints JSON lines of agent outputs.
    """
    sys.stderr.write("[INFO] No arguments provided. Running built-in example...\n")
    reports = generate_reports_via_llm(3, model=model)
    router = Router()
    for rpt in reports:
        info = classify_report_llm(rpt.text, model=model)
        priority = router._priority_from_severity(info.get("severity", "unknown"))  # type: ignore[arg-type]
        heapq.heappush(router._pq, (priority, time.monotonic(), rpt, info))
        sys.stderr.write(
            f"[INGEST] id={rpt.id} label={info.get('label')} severity={info.get('severity')} priority={priority}\n"
        )
    while not router.empty():
        out = router.dispatch_next(enrich=enrich)
        if out is None:
            break
        print(out)


# --- CLI ------------------------------------------------------------------------
def parse_args(argv: Optional[List[str]] = None) -> Optional[argparse.Namespace]:
    """
    Returns argparse.Namespace or None if argv is empty (to allow auto-example).
    """
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) == 0:
        return None

    p = argparse.ArgumentParser(description="Pharmacovigilance Inbox Demo (single-file)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--input", type=str, help="Path to JSONL file with reports")
    src.add_argument("--stdin", action="store_true", help="Read reports from stdin (interactive)")

    p.add_argument("--enrich", action="store_true", help="Use optional read-only LS tools when available")
    p.add_argument("--drain", action="store_true", help="After ingest, drain the queue and print JSON summaries")
    p.add_argument("--model", type=str, default="gpt-4o-mini", help="LLM model name")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    ns = parse_args(argv)

    # Auto-example mode if no args
    if ns is None:
        model = os.getenv("PV_MODEL_NAME", "gpt-4o-mini")
        enrich = os.getenv("PV_EXAMPLE_ENRICH", "0") == "1"
        run_auto_example(enrich=enrich, model=model)
        return

    model = os.getenv("PV_MODEL_NAME", ns.model)
    router = Router()

    # Ingest stream
    stream: Iterable[Report]
    if ns.input:
        stream = read_jsonl(ns.input)
    else:
        stream = read_stdin()

    for rpt in stream:
        info = classify_report_llm(rpt.text, model=model)
        priority = router._priority_from_severity(info.get("severity", "unknown"))  # type: ignore[arg-type]
        heapq.heappush(router._pq, (priority, time.monotonic(), rpt, info))
        sys.stderr.write(
            f"[INGEST] id={rpt.id} label={info.get('label')} severity={info.get('severity')} priority={priority}\n"
        )

    # Drain queue in priority order
    if ns.drain or True:
        while not router.empty():
            out = router.dispatch_next(enrich=ns.enrich)
            if out is None:
                break
            print(out)


if __name__ == "__main__":
    main()

