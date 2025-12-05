# lesson-2-multi-agent-architecture-implementation/demo/life_sciences_demo.py
# -----------------------------------------------------------------------------
# Life Sciences Evidence Scout — Demo
#
# What this file shows:
# - Minimal multi-agent scaffolding (Agent base, Orchestrator)
# - Explicit JSON serialization via a MessageBus (visible, testable)
# - Real, read-only tools for live data: ClinicalTrials.gov v2, PubMed
# - LLM synthesis via OpenAI (Vocareum base URL) with .env API key
# - Optional stretch agents for FAERS and ClinVar (off by default)
#
# -----------------------------------------------------------------------------

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

# ---------- LLM client (Vocareum endpoint) ----------
from dotenv import load_dotenv
from openai import OpenAI

load_dotenv("../.env")  # expects OPENAI_API_KEY in your .env
OPENAI_BASE_URL = "https://openai.vocareum.com/v1"
client = OpenAI(base_url=OPENAI_BASE_URL)

# ---------- Read-only tool wrappers (live APIs; provided by the course) ----------
from ls_action_space.action_space import (
    query_clinicaltrials,  # ClinicalTrials.gov v2
    query_pubmed,          # PubMed (NCBI)
    # Optional tools — only used if you enable stretch agents:
    # query_clinvar,
    # get_faers_contingency_counts,
)

# =============================================================================
# Infrastructure: Blackboard + MessageBus + Agent base
# =============================================================================

@dataclass
class Blackboard:
    """Shared key-value memory for the overall run."""
    topic: Dict[str, Any] = field(default_factory=dict)     # e.g., {"query": "..."}
    synthesis: str = ""                                     # final LLM brief
    logs: List[Dict[str, Any]] = field(default_factory=list)


class MessageBus:
    """
    Very small message bus that *explicitly serializes* envelopes to JSON strings.
    This makes serialization/deserialization visible and easy to unit-test.
    """
    def __init__(self) -> None:
        self._msgs: List[str] = []

    def send(self, sender: str, topic: str, payload: Dict[str, Any]) -> None:
        envelope = {"sender": sender, "topic": topic, "payload": payload}
        self._msgs.append(json.dumps(envelope))

    def drain(self, topic: str) -> List[Dict[str, Any]]:
        out: List[Dict[str, Any]] = []
        keep: List[str] = []
        for m in self._msgs:
            d = json.loads(m)
            if d["topic"] == topic:
                out.append(d)
            else:
                keep.append(m)
        self._msgs = keep
        return out


class Agent:
    """Abstract base for all agents: implement decide() and act()."""
    def __init__(self, name: str, board: Blackboard, bus: MessageBus) -> None:
        self.name = name
        self.board = board
        self.bus = bus

    def decide(self) -> bool:
        raise NotImplementedError

    def act(self) -> None:
        raise NotImplementedError


# =============================================================================
# Domain Agents (3-Agent Core): Trials, PubMed, LLM Synthesis
# =============================================================================

class TrialScoutAgent(Agent):
    """
    Fetches a compact set of trial fields from ClinicalTrials.gov v2 based on a query.
    Emits a JSON envelope on topic "trials.results".
    """
    def __init__(self, name: str, board: Blackboard, bus: MessageBus, max_results: int = 3) -> None:
        super().__init__(name, board, bus)
        self.max_results = max_results
        self._done = False

    def decide(self) -> bool:
        return bool(self.board.topic.get("query")) and not self._done

    def act(self) -> None:
        expr = self.board.topic["query"]
        fields = [
            "NCTId",
            "BriefTitle",
            "OverallStatus",
            "StudyType",
            "Phase",
            "PrimaryOutcomeMeasure",
            "LeadSponsorName",
            "StartDate",
            "CompletionDate",
        ]
        try:
            resp = query_clinicaltrials(expr=expr, max_results=self.max_results, fields=fields)
            trials = resp.get("studies", []) or []
        except Exception as e:
            trials = []
            self.board.logs.append({"agent": self.name, "level": "warn", "msg": f"clinicaltrials error: {e}"})
        self.bus.send(self.name, "trials.results", {"trials": trials})
        self._done = True


class LitReviewAgent(Agent):
    """
    Fetches a compact set of PubMed articles matching the query.
    Emits a JSON envelope on topic "literature.results".
    """
    def __init__(self, name: str, board: Blackboard, bus: MessageBus, max_results: int = 3) -> None:
        super().__init__(name, board, bus)
        self.max_results = max_results
        self._done = False

    def decide(self) -> bool:
        return bool(self.board.topic.get("query")) and not self._done

    def act(self) -> None:
        q = self.board.topic["query"]
        try:
            arts = query_pubmed(q, max_results=self.max_results, include_mesh=True, include_citations=False) or []
        except Exception as e:
            arts = []
            self.board.logs.append({"agent": self.name, "level": "warn", "msg": f"pubmed error: {e}"})
        self.bus.send(self.name, "literature.results", {"articles": arts})
        self._done = True


class LLMReviewerAgent(Agent):
    """
    Consumes trial and literature envelopes, synthesizes a brief.
    Writes result into board.synthesis.
    """
    SYS_PROMPT = (
        "You are a biomedical brief writer. Synthesize trials and literature into a crisp summary. "
        "Include identifiers (e.g., NCT IDs, DOIs/PMIDs) and years where available. "
        "Flag uncertainties and data gaps succinctly."
    )

    def __init__(self, name: str, board: Blackboard, bus: MessageBus) -> None:
        super().__init__(name, board, bus)
        self._done = False

    def decide(self) -> bool:
        return not self._done

    def _collect_payload(self) -> Dict[str, Any]:
        trials_env = self.bus.drain("trials.results")
        lits_env = self.bus.drain("literature.results")
        trials = []
        for e in trials_env:
            trials.extend(e["payload"].get("trials", []))
        articles = []
        for e in lits_env:
            articles.extend(e["payload"].get("articles", []))
        return {"query": self.board.topic.get("query"), "trials": trials, "articles": articles}

    def act(self) -> None:
        payload = self._collect_payload()
        # The LLM call is minimal and deterministic (temperature 0.2)
        try:
            resp = client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {"role": "system", "content": self.SYS_PROMPT},
                    {"role": "user", "content": f"Summarize succinctly:\n{json.dumps(payload, ensure_ascii=False)}"},
                ],
                temperature=0.2,
            )
            self.board.synthesis = (resp.choices[0].message.content or "").strip()
        except Exception as e:
            self.board.synthesis = ""
            self.board.logs.append({"agent": self.name, "level": "error", "msg": f"LLM error: {e}"})
        self._done = True


# =============================================================================
# Optional Stretch Agents (disabled by default)
#   - SafetySignalAgent (FAERS ROR)
#   - VariantAgent (ClinVar)
# Enable via CLI flags to demonstrate extensibility beyond the core lesson.
# =============================================================================

class SafetySignalAgent(Agent):
    """
    OPTIONAL: Builds a simple disproportionality (ROR) using FAERS counts
    for (drug, adverse_event). Emits envelope on topic "safety.results".
    Requires get_faers_contingency_counts to be available.
    """
    def __init__(self, name: str, board: Blackboard, bus: MessageBus) -> None:
        super().__init__(name, board, bus)
        self._done = False
        # Lazy import to avoid dependency when not used
        try:
            from ls_action_space.action_space import get_faers_contingency_counts  # type: ignore
            self._faers = get_faers_contingency_counts
        except Exception as e:
            self._faers = None
            self.board.logs.append({"agent": self.name, "level": "warn", "msg": f"FAERS helper unavailable: {e}"})

    def decide(self) -> bool:
        t = self.board.topic
        return bool(t.get("drug") and t.get("adverse_event")) and not self._done

    @staticmethod
    def _ror(a: int, b: int, c: int, d: int) -> Dict[str, float]:
        import math
        a, b, c, d = [max(0, int(x)) for x in (a, b, c, d)]
        ror = ((a + 0.5) / (b + 0.5)) / ((c + 0.5) / (d + 0.5))
        se = (1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5)) ** 0.5
        lcl = float(__import__("math").exp(__import__("math").log(ror) - 1.96 * se))
        ucl = float(__import__("math").exp(__import__("math").log(ror) + 1.96 * se))
        return {"ROR": float(ror), "LCL95": lcl, "UCL95": ucl}

    def act(self) -> None:
        if not self._faers:
            self.bus.send(self.name, "safety.results", {"error": "faers_helper_unavailable"})
            self._done = True
            return
        drug = self.board.topic["drug"]
        ae = self.board.topic["adverse_event"]
        try:
            counts = self._faers(drug, ae)  # returns dict a,b,c,d
            metrics = self._ror(**counts)
            self.bus.send(self.name, "safety.results", {"event": ae, "counts": counts, "metrics": metrics})
        except Exception as e:
            self.bus.send(self.name, "safety.results", {"error": f"faers_error: {e}"})
            self.board.logs.append({"agent": self.name, "level": "warn", "msg": f"FAERS error: {e}"})
        self._done = True


class VariantAgent(Agent):
    """
    OPTIONAL: Looks up a ClinVar record for a variant (rsID or HGVS string).
    Emits envelope on topic "variant.results".
    """
    def __init__(self, name: str, board: Blackboard, bus: MessageBus) -> None:
        super().__init__(name, board, bus)
        self._done = False
        try:
            from ls_action_space.action_space import query_clinvar  # type: ignore
            self._clinvar = query_clinvar
        except Exception as e:
            self._clinvar = None
            self.board.logs.append({"agent": self.name, "level": "warn", "msg": f"ClinVar helper unavailable: {e}"})

    def decide(self) -> bool:
        return bool(self.board.topic.get("variant")) and not self._done

    def act(self) -> None:
        if not self._clinvar:
            self.bus.send(self.name, "variant.results", {"error": "clinvar_helper_unavailable"})
            self._done = True
            return
        try:
            result = self._clinvar(self.board.topic["variant"])
        except Exception as e:
            result = {"error": f"clinvar_error: {e}"}
            self.board.logs.append({"agent": self.name, "level": "warn", "msg": f"ClinVar error: {e}"})
        self.bus.send(self.name, "variant.results", {"variant": result})
        self._done = True


# =============================================================================
# Orchestrator
# =============================================================================

class Orchestrator:
    """Simple round-robin orchestrator until no agent wants to act."""
    def __init__(self, agents: List[Agent]) -> None:
        self.agents = agents

    def run(self, max_iters: int = 5) -> None:
        for _ in range(max_iters):
            progressed = False
            for agent in self.agents:
                try:
                    if agent.decide():
                        agent.act()
                        progressed = True
                except Exception as e:
                    # Keep the system resilient: log and continue
                    if hasattr(agent, "board"):
                        agent.board.logs.append({"agent": agent.name, "level": "error", "msg": f"crash: {e}"})
            if not progressed:
                break


# =============================================================================
# Public API: run_demo
# =============================================================================

def run_demo(
    query: str,
    *,
    max_results: int = 3,
    include_safety: bool = False,
    drug_for_safety: Optional[str] = None,
    adverse_event: Optional[str] = None,
    include_variant: bool = False,
    variant: Optional[str] = None,
) -> Blackboard:
    """
    Runs the demo for a given query and optional stretch tasks.

    Args:
        query: Free-text query (e.g., disease, drug) used for trials and PubMed.
        max_results: Cap on results per tool for determinism and speed.
        include_safety: If True, adds SafetySignalAgent (requires drug + adverse_event).
        drug_for_safety: Drug name for FAERS analysis (only if include_safety=True).
        adverse_event: MedDRA PT for FAERS analysis (only if include_safety=True).
        include_variant: If True, adds VariantAgent (requires variant string).
        variant: rsID (e.g., "rs429358") or HGVS string.

    Returns:
        Blackboard containing topic, synthesis, and logs.
    """
    board = Blackboard(
        topic={
            "query": query,
            "drug": drug_for_safety if include_safety else None,
            "adverse_event": adverse_event if include_safety else None,
            "variant": variant if include_variant else None,
        }
    )
    bus = MessageBus()

    agents: List[Agent] = [
        TrialScoutAgent("trials", board, bus, max_results=max_results),
        LitReviewAgent("pubmed", board, bus, max_results=max_results),
        LLMReviewerAgent("llm", board, bus),
    ]

    if include_safety and drug_for_safety and adverse_event:
        agents.insert(2, SafetySignalAgent("faers", board, bus))  # run before LLM to include signal

    if include_variant and variant:
        agents.insert(2, VariantAgent("clinvar", board, bus))     # run before LLM to include variant context

    Orchestrator(agents).run()
    return board


# =============================================================================
# CLI
# =============================================================================

def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Life Sciences Evidence Scout (3-Agent Core with optional stretch agents)"
    )
    # Allow running with no argv by providing a default (env overrides)
    p.add_argument(
        "--query",
        required = False,
        default = os.getenv("DEMO_QUERY", "semaglutide obesity randomized"),
        help = "Free-text query for ClinicalTrials.gov and PubMed. "
         "(Default from DEMO_QUERY env var or a built-in query.)",
    )
    p.add_argument("--max_results", type=int, default=3, help="Max results per tool (default: 3).")

    # Optional FAERS
    p.add_argument("--include_safety", action="store_true", help="Enable FAERS SafetySignalAgent.")
    p.add_argument("--drug_for_safety", type=str, help="Drug name for FAERS.")
    p.add_argument("--adverse_event", type=str, help="MedDRA Preferred Term for FAERS.")

    # Optional ClinVar
    p.add_argument("--include_variant", action="store_true", help="Enable ClinVar VariantAgent.")
    p.add_argument("--variant", type=str, help="Variant query (rsID or HGVS).")

    # Output
    p.add_argument("--print_payloads", action="store_true", help="Print final synthesis and logs as JSON.")
    return p


def main() -> None:
    args = _build_arg_parser().parse_args()

    board = run_demo(
        query=args.query,
        max_results=args.max_results,
        include_safety=args.include_safety,
        drug_for_safety=args.drug_for_safety,
        adverse_event=args.adverse_event,
        include_variant=args.include_variant,
        variant=args.variant,
    )

    # Human-friendly output
    print("\n" + "=" * 80)
    print("LIFE SCIENCES EVIDENCE SCOUT — SYNTHESIS")
    print("=" * 80)
    if board.synthesis:
        print(board.synthesis)
    else:
        print("(No synthesis produced. See logs.)")

    # Optional JSON dump for grading/debugging
    if args.print_payloads:
        print("\n" + "=" * 80)
        print("DEBUG LOGS (JSON)")
        print("=" * 80)
        print(json.dumps({"topic": board.topic, "logs": board.logs}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

