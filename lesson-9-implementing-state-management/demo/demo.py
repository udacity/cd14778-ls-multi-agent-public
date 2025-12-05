#!/usr/bin/env python3
"""
Demo: Lab Notebook Logger 

- A SINGLE shared `lab_state` dict is passed by reference to all agents.
- EntryAgent WRITES experiment notes.
- SummaryAgent READS notes and WRITES a TL;DR back into the same state.
- A tiny CLI and a default "no-args" demo runner.

.env layout:
OPENAI_API_KEY=your_key_here
"""

from __future__ import annotations
import os
import sys
import argparse
from typing import Dict, Any, List, Optional

# --- load API key from .env if present --------------------------------
try:
    from dotenv import load_dotenv
    load_dotenv("../.env")
except Exception:
    # It's fine if dotenv isn't installed; we'll fall back to offline summary.
    pass

# --- OpenAI client (Vocareum base URL). Falls back if unavailable. ----
OpenAIClient = None
try:
    # Using the official OpenAI Python client interface
    from openai import OpenAI  # type: ignore

    def _maybe_make_client() -> Optional["OpenAI"]:
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            return None
        # Client reads key from environment; base_url per instructions.
        return OpenAI(base_url="https://openai.vocareum.com/v1")

    OpenAIClient = _maybe_make_client()

except Exception:
    # No client available (package missing or no key). We'll use offline summary.
    OpenAIClient = None


# ================================================================================
# Shared STATE (single source of truth)
# ================================================================================
lab_state: Dict[str, Any] = {
    "notes": [],          # list[dict]: {"text": str, "author": str, "tag": Optional[str]}
    "daily_summary": "",  # str: produced by SummaryAgent
}


# ================================================================================
# Agents
# ================================================================================
class EntryAgent:
    """Writes free-text experiment notes into the shared notebook."""
    def __init__(self, state: Dict[str, Any]):
        self.state = state  # reference to shared state

    def add_note(self, text: str, author: str = "anonymous", tag: Optional[str] = None) -> None:
        self.state["notes"].append({"text": text, "author": author, "tag": tag})
        print(f"[EntryAgent] Note added by {author}. Total notes: {len(self.state['notes'])}")


class SummaryAgent:
    """
    Reads all notes and writes a concise TL;DR back into the shared state.
    If an OpenAI client is available, uses it. Otherwise, falls back to a
    simple deterministic summarizer so the demo always runs.
    """
    def __init__(self, state: Dict[str, Any], llm_client: Optional[Any] = None):
        self.state = state
        self.client = llm_client

    def _llm_summary(self, notes: List[Dict[str, Any]]) -> str:
        corpus = "\n".join(f"- {n['text']}" for n in notes)
        prompt = (
            "You are a lab assistant. Produce a concise, factual TL;DR "
            "(5–8 bullet points) summarizing the experiment notes. Avoid "
            "speculation; highlight outcomes, parameters, anomalies, and "
            "next steps.\n\nNotes:\n" + corpus
        )
        resp = self.client.responses.create(
            model="gpt-4o-mini",
            input=[{"role": "user", "content": prompt}],
        )
        return resp.output_text

    def _offline_summary(self, notes: List[Dict[str, Any]]) -> str:
        # Very small, deterministic, no-ML fallback to keep demo runnable.
        bullets = []
        if not notes:
            return "- No notes recorded today."
        # Extract up to first 8 concise statements.
        for n in notes[:8]:
            txt = n["text"].strip()
            # Keep the first sentence-ish under ~140 chars.
            cut = txt.split(". ")[0]
            if len(cut) > 140:
                cut = cut[:137].rstrip() + "..."
            bullets.append(f"- {cut}.")
        # Heuristic "next step" if any measurement/observation keywords appear
        joined = " ".join(n["text"] for n in notes).lower()
        if any(k in joined for k in ("observe", "signal", "yield", "growth", "pcr", "gel", "fluorescence", "assay")):
            bullets.append("- Next: replicate key conditions; record controls and parameters.")
        return "\n".join(bullets)

    def generate_daily_summary(self) -> None:
        notes = self.state.get("notes", [])
        if self.client:
            try:
                summary = self._llm_summary(notes)
            except Exception as e:
                print(f"[SummaryAgent] LLM failed ({e}); using offline summary.")
                summary = self._offline_summary(notes)
        else:
            summary = self._offline_summary(notes)

        self.state["daily_summary"] = summary
        print("[SummaryAgent] Daily summary updated.")


# ================================================================================
# CLI + Demo Runner
# ================================================================================
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Lab Notebook Logger — state-sharing multi-agent demo."
    )
    sub = p.add_subparsers(dest="cmd")

    # add-note
    add = sub.add_parser("add-note", help="Append a note to the shared notebook.")
    add.add_argument("--text", required=True, help="Free-text note contents.")
    add.add_argument("--author", default="anonymous", help="Author name.")
    add.add_argument("--tag", default=None, help="Optional tag (e.g., 'cell_culture').")

    # summarize
    sub.add_parser("summarize", help="Generate end-of-day TL;DR.")

    # show
    sub.add_parser("show", help="Print current notes and summary.")

    return p


def _cmd_add_note(args) -> None:
    entry = EntryAgent(lab_state)
    entry.add_note(args.text, author=args.author, tag=args.tag)


def _cmd_summarize(_args) -> None:
    summarizer = SummaryAgent(lab_state, OpenAIClient)
    summarizer.generate_daily_summary()


def _cmd_show(_args) -> None:
    print("\n=== NOTES ===")
    for i, n in enumerate(lab_state["notes"], 1):
        who = n.get("author") or "anonymous"
        tag = f" [{n['tag']}]" if n.get("tag") else ""
        print(f"{i:02d}. {n['text']}  — {who}{tag}")
    print("\n=== DAILY SUMMARY ===")
    print(lab_state.get("daily_summary", "").strip() or "(empty)")
    print()


def demo_run() -> None:
    """
    Runs automatically if the script is executed without arguments.
    Shows WRITE -> READ/WRITE behavior over ONE shared state object.
    """
    print("[Demo] Running default scenario (no CLI args provided).")
    entry = EntryAgent(lab_state)
    summary = SummaryAgent(lab_state, OpenAIClient)

    # During the day: capture notes (example data kept simple for clarity).
    entry.add_note(
        "Seeded HEK293 cells at ~70% confluence; transfected pGFP with "
        "lipid reagent; incubated at 37°C, 5% CO₂.",
        author="A. Kim",
        tag="cell_culture",
    )
    entry.add_note(
        "At 24h post-transfection, observed GFP fluorescence approx. 2x "
        "vs. mock control. No overt cytotoxicity.",
        author="A. Kim",
        tag="readout",
    )

    # End of day: produce TL;DR (READ notes -> WRITE daily_summary)
    summary.generate_daily_summary()

    # Inspect shared state (both agents operated on the SAME dict)
    _cmd_show(None)


def main(argv: List[str]) -> None:
    if not argv:
        demo_run()
        return

    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "add-note":
        _cmd_add_note(args)
    elif args.cmd == "summarize":
        _cmd_summarize(args)
    elif args.cmd == "show":
        _cmd_show(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main(sys.argv[1:])

