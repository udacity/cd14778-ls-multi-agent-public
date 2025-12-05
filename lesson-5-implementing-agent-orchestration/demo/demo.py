"""
Core-Lab Service Desk — Single-File Demo

- A tiny agent framework (BaseAgent + 3 specialists + Intake + Orchestrator)
- Sequential → conditional routing → delegation → follow-ups (quote + checklist)
- Clean seams for tools (e.g., PubMed), pricing, and policies without hardcoded data

How to run:
1) pip install python-dotenv openai
2) Put OPENAI_API_KEY in a .env file next to this script.
   Optional: set PRICE_TABLE (JSON), CURRENCY, and CHECKLIST_TEMPLATE (JSON) in .env
3) python core_lab_service_desk_demo.py
"""

from __future__ import annotations
import json
import os
import sys
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Optional, Protocol, List
from ls_action_space.action_space import (
    query_pubmed, query_clinicaltrials, query_clinvar
)

# --- LLM client (optional but wired) -------------------------------------------------------------
# Uses Vocareum proxy, per instructions. Gracefully degrades if no key is present.
try:
    from dotenv import load_dotenv
    from openai import OpenAI
    load_dotenv("../.env")
    _OPENAI_READY = bool(os.getenv("OPENAI_API_KEY"))
except Exception:
    print("FAILED to load OpenAI key")
    _OPENAI_READY = False
    OpenAI = None  # type: ignore

LLM_MODEL = "gpt-4o-mini"

def get_llm():
    if _OPENAI_READY and OpenAI is not None:
        try:
            return OpenAI(base_url="https://openai.vocareum.com/v1")
        except Exception:
            pass
    return None  # Fallback: no external calls; agents will produce sensible strings.


# --- Tool protocol (read-only) -------------------------------------------------------------------
# In your real environment, you could import these from ls_action_space.action_space:
#   from ls_action_space.action_space import query_pubmed, query_clinicaltrials, query_clinvar
# Here we keep a clean seam (callables in a dict). If absent, agents skip tool calls.
ToolFn = Callable[..., Any]


# --- Minimal agent framework ---------------------------------------------------------------------
class RunnableAgent(Protocol):
    def run(self, task: str, context: Optional[Dict[str, Any]] = None) -> str: ...


@dataclass
class BaseAgent:
    name: str
    llm: Any = None
    tools: Dict[str, ToolFn] = field(default_factory=dict)

    def _llm(self, prompt: str, **kwargs) -> Optional[str]:
        """
        Small helper: If LLM available, get a short completion.
        Kept intentionally minimal; no hardcoded content.
        """
        if not self.llm:
            return None
        try:
            resp = self.llm.chat.completions.create(
                model=LLM_MODEL,
                messages=[{"role": "system", "content": f"You are {self.name}, a helpful specialist."},
                          {"role": "user", "content": prompt}],
                temperature=0.2,
                max_tokens=300,
            )
            return resp.choices[0].message.content or ""
        except Exception:
            return None

    def run(self, task: str, context: Optional[Dict[str, Any]] = None) -> str:
        """
        Default run: forward to LLM if present; otherwise return a concise, template-based response.
        This keeps the demo runnable without external dependencies while avoiding hardcoded domain data.
        """
        ctx = context or {}
        prompt = f"{task}\n\nContext:\n{json.dumps(ctx, indent=2, ensure_ascii=False)}"
        out = self._llm(prompt)
        if out:
            return out.strip()
        # Fallback (no LLM): return a short synthesized message.
        summary_bits = [f"{k}={v}" for k, v in ctx.items()]
        return f"[{self.name}] {task.strip()} :: " + (", ".join(summary_bits) if summary_bits else "no-context")


# --- Intake + Specialists ------------------------------------------------------------------------
class IntakeAgent(BaseAgent):
    def diagnose(self, user_request: str) -> Dict[str, Any]:
        """
        Extract coarse intent and key entities. No hardcoded rules beyond simple keyword cues;
        students can replace with an LLM prompt or a classifier.
        """
        req = user_request.lower()
        if any(w in req for w in ["rna", "seq", "ngs", "illumina", "library"]):
            intent = "sequencing"
        elif any(w in req for w in ["ms", "mass spec", "proteomics", "dia", "dda", "lc-ms"]):
            intent = "mass_spec"
        elif any(w in req for w in ["image", "confocal", "microscop", "spinning", "objective"]):
            intent = "imaging"
        else:
            intent = "unknown"

        # Entity extraction: keep generic; do not embed fixed values.
        entities = {
            "organism": None,
            "modality": None,
            "goal": None,
            "sample_count": None,
            "timing": None,
            "free_text": user_request.strip(),
        }
        return {"intent": intent, "entities": entities}

    def initial_response(self, diagnosis: Dict[str, Any]) -> str:
        intent = diagnosis.get("intent", "unknown")
        if intent == "unknown":
            return ("Thanks for the details—I'll need a bit more info to route this properly. "
                    "If you can, include organism, modality, goal, and approximate sample count.")
        return f"Got it—this looks like a {intent.replace('_', ' ')} request. I’ll route it to the right specialist."


class SequencingAgent(BaseAgent):
    def advise(self, brief: Dict[str, Any]) -> str:
        """
        Provide modality/depth/QC framing. Optionally consult tools if provided.
        Keep guidance abstract and parameterized; do not embed lab-specific numbers.
        """
        # Optional literature grounding (skipped if tool not present):
        if "pubmed" in self.tools and brief:
            try:
                _ = self.tools["pubmed"](query="generated at runtime", max_results=3)  # placeholder pattern
            except Exception:
                pass

        task = ("Draft guidance on assay selection (e.g., RNA-seq/ATAC-seq/WGS as appropriate), "
                "suggest read depth ranges qualitatively, and outline QC checkpoints "
                "(RNA integrity, library metrics, contamination checks). "
                "Parameterize recommendations based on the brief.")
        return self.run(task, {"brief": brief})


class MassSpecAgent(BaseAgent):
    def advise(self, brief: Dict[str, Any]) -> str:
        """
        Recommend proteomics workflow selection (DDA/DIA/SRM), enrichment options, and QC checkpoints.
        No hardcoded gradients/instrument models.
        """
        if "pubmed" in self.tools and brief:
            try:
                _ = self.tools["pubmed"](query="generated at runtime", max_results=3)
            except Exception:
                pass

        task = ("Outline whether DDA, DIA, or targeted SRM/PRM is most suitable, note sample prep considerations "
                "(e.g., lysis, digestion, optional enrichment), and list QC/validation steps. "
                "Keep details adaptable to the brief.")
        return self.run(task, {"brief": brief})


class ImagingAgent(BaseAgent):
    def advise(self, brief: Dict[str, Any]) -> str:
        """
        Guide modality (widefield/confocal/spinning-disk), objectives, channels, environment control,
        and acquisition parameters at a conceptual level—no fixed instrument presets.
        """
        task = ("Suggest an imaging modality and objective class based on the biological question, "
                "recommend channel layout and exposure/laser considerations qualitatively, "
                "and list controls (e.g., fluorophore separation, live-cell conditions).")
        return self.run(task, {"brief": brief})


# --- Orchestrator (sequential + conditional + delegation + follow-ups) ---------------------------
class Orchestrator(BaseAgent):
    def __init__(
        self,
        llm: Any = None,
        tools: Optional[Dict[str, ToolFn]] = None,
        pricing_model: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None,
        checklist_policy: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None,
    ):
        super().__init__("Orchestrator", llm, tools or {})
        self.intake = IntakeAgent("Intake", llm, tools or {})
        self.seq = SequencingAgent("Sequencing", llm, tools or {})
        self.ms = MassSpecAgent("MassSpec", llm, tools or {})
        self.img = ImagingAgent("Imaging", llm, tools or {})

        # Dependencies are injected; no hardcoded pricing or SOPs.
        self._pricing_model = pricing_model or default_pricing_model
        self._checklist_policy = checklist_policy or default_checklist_policy

    def handle_request(self, user_request: str) -> str:
        # 1) Sequential: diagnose first
        diagnosis = self.intake.diagnose(user_request)

        # 2) Sequential: brief acknowledgment / expectation setting
        initial = self.intake.initial_response(diagnosis)

        # 3) Conditional branching: choose specialist
        intent = diagnosis.get("intent", "unknown")
        brief = diagnosis.get("entities", {})

        if intent == "sequencing":
            advice = self.seq.advise(brief)
        elif intent == "mass_spec":
            advice = self.ms.advise(brief)
        elif intent == "imaging":
            advice = self.img.advise(brief)
        else:
            return initial + "\n\nI couldn’t match this to sequencing, mass spectrometry, or imaging. " \
                             "Please add organism, modality, goal, and rough sample count."

        # 4) Follow-ups (delegated helpers): quote + checklist
        quote = self._draft_quote(brief)
        checklist = self._build_checklist(brief)

        return f"{initial}\n\n{advice}\n\n---\nProvisional Quote\n{quote}\n\nSample-Prep Checklist\n{checklist}"

    # --- Follow-ups ------------------------------------------------------------------------------

    def _draft_quote(self, brief: Dict[str, Any]) -> str:
        """
        Turn a dict from the pricing model into a readable quote without hardcoded numbers here.
        """
        pricing_dict = self._pricing_model(brief)  # externalized computation
        return self.run("Summarize this pricing dictionary as a clear, itemized quote.",
                        {"pricing": pricing_dict})

    def _build_checklist(self, brief: Dict[str, Any]) -> str:
        """
        Format a checklist produced by a policy function; optionally cross-reference tools.
        """
        policy = self._checklist_policy(brief)
        # Optionally: verify steps with a literature tool if available (skipped by default).
        return self.run("Format this policy as a concise, numbered checklist with short rationales.",
                        {"policy": policy})


# --- Default dependencies (no hardcoded prices or SOPs in code) ---------------------------------
def default_pricing_model(brief: Dict[str, Any]) -> Dict[str, Any]:
    """
    Pricing is parameterized and externalized:
    - PRICE_TABLE (JSON) may come from env: {"sequencing":{"per_sample":...}, "mass_spec":...}
    - CURRENCY from env (default "USD")
    If PRICE_TABLE absent, we still return a structured, label-only quote scaffold.
    """
    currency = os.getenv("CURRENCY", "USD")
    try:
        price_table = json.loads(os.getenv("PRICE_TABLE", "{}"))
    except Exception:
        price_table = {}

    # Very light, data-driven sketch—no embedded numbers here.
    # The UI or env can supply actual rates at runtime.
    base = {
        "currency": currency,
        "line_items": [],
        "notes": ["Final price depends on validated sample count and run parameters."]
    }

    # Build line items based on intent if the intake provided it in brief (optional).
    intent = brief.get("intent") or "unspecified"
    sample_count = brief.get("sample_count") or "TBD"

    # Pull rates (if provided). We only copy keys; we do not embed concrete values.
    rates = price_table.get(intent, {})

    base["line_items"].append({
        "label": f"{intent} core services",
        "quantity": sample_count,
        "unit_price_from_env": list(rates.keys()) or ["see rate card"],
        "computed_subtotal": "computed externally"
    })
    base["line_items"].append({
        "label": "QC & data delivery",
        "quantity": sample_count,
        "unit_price_from_env": ["see rate card"],
        "computed_subtotal": "computed externally"
    })
    return base


def default_checklist_policy(brief: Dict[str, Any]) -> Dict[str, Any]:
    """
    Checklist template can be supplied via CHECKLIST_TEMPLATE (JSON). If absent, return a generic,
    parameterized skeleton without lab-specific SOPs or kit SKUs.
    """
    try:
        template = json.loads(os.getenv("CHECKLIST_TEMPLATE", "{}"))
        if isinstance(template, dict) and template:
            return template
    except Exception:
        pass

    return {
        "title": "Sample-Prep Checklist (Provisional)",
        "steps": [
            {"step": "Confirm sample metadata",
             "why": "Routing, safety, and downstream parameterization depend on accurate metadata.",
             "inputs": ["organism", "tissue/cell type", "preservation method"]},
            {"step": "Assess input quality",
             "why": "Quality gates prevent wasted runs and improve reproducibility.",
             "inputs": ["concentration", "integrity metrics (qualitative)"]},
            {"step": "Prepare according to selected modality",
             "why": "Prep differs by modality; refer to current facility guidance.",
             "inputs": ["assay-specific requirements (abstract)"]},
            {"step": "Include controls",
             "why": "Controls are critical for interpretation.",
             "inputs": ["negative/positive controls as appropriate"]},
            {"step": "Package & submit samples",
             "why": "Ensure chain-of-custody and safe transport.",
             "inputs": ["labels, form, container requirements"]},
        ]
    }


# --- Simple CLI runner ---------------------------------------------------------------------------
def main():
    llm = get_llm()

    # Tools
    tools: Dict[str, ToolFn] = {
        "pubmed": query_pubmed,               # usage: tools["pubmed"](query="term", max_results=10)
        "clinicaltrials": query_clinicaltrials,  # usage: tools["clinicaltrials"](expr="condition", max_results=25)
        "clinvar": query_clinvar,             # usage: tools["clinvar"]("rs123")
    }

    orch = Orchestrator(llm=llm, tools=tools)

    # Read a single request from stdin or from CLI args
    if len(sys.argv) > 1:
        user_request = " ".join(sys.argv[1:])
    else:
        print("Free text experiment description.")
        print("Examples: e.g.")
        print("""
Sequencing: RNA-seq on 12 human samples; goal: differential expression
Mass spec / Proteomics: LC-MS DIA proteomics, mouse liver n=8; workflow advice
Imaging: Confocal live-cell imaging of mitochondria; 63x with CO₂/temp control
Unknown / Clarification: Not sure if sequencing or proteomics fits; 10 samples—please advise
Edge case / Robustness: RNA library QC?? 12 samples; urgent timeline
        """)
        print("Describe your experiment request (free text), then press Enter:")
        user_request = sys.stdin.readline().strip()
        print("Processing your request: ", user_request)

    # Inject minimal fields back into brief so pricing can reflect intent/sample_count if user provided it later.
    # For the demo, we simply run orchestrator and print.
    result = orch.handle_request(user_request)
    print("\n=== Orchestration Result ===\n")
    print(result)


if __name__ == "__main__":
    main()