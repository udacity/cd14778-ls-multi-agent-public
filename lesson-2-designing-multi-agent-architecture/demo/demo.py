# -*- coding: utf-8 -*-
"""
Variant & Evidence Concierge (Life Sciences)
- Orchestrator pattern with role-specialized agents
- No hardcoded data: all answers come from read-only public APIs via helper functions
- Uses OpenAI (vocareum base URL) for light NLU and optional response drafting/localization

Setup:
  1) pip install python-dotenv biopython
  2) Create a .env file with: OPENAI_API_KEY=sk-...
  3) Ensure ls_action_space.action_space is importable (provided by the course environment)

Run:
  python life_sciences_demo.py
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass
from dotenv import load_dotenv
import os

# --- OpenAI client (vocareum base URL) ---
from openai import OpenAI
load_dotenv("../.env")
client = OpenAI(base_url="https://openai.vocareum.com/v1")

# --- Read-only API wrappers already provided for you ---
# (they call public services under the hood; no state is mutated)
from ls_action_space.action_space import (
    query_clinvar,          # variant -> ClinVar + (rsID) allele frequencies + PubMed links
    query_pubmed,           # query -> list of articles w/ metadata
    query_clinicaltrials    # expr  -> ClinicalTrials.gov v2 projected fields
)

# ----------------------------
# Language + Intent Detection
# ----------------------------

class LanguageAgent:
    """Detects the user's language using a lightweight LLM call."""
    def detect_language(self, text: str) -> str:
        msg = [
            {"role": "system", "content": "Detect language code (ISO 639-1). Reply with only the code, e.g., 'en'."},
            {"role": "user", "content": text}
        ]
        out = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=msg,
            temperature=0
        )
        code = out.choices[0].message.content.strip().lower()
        return code or "en"


class IntentAgent:
    """
    Classifies the user’s question into one of:
      - 'frequency'      (allele frequencies for a variant)
      - 'pathogenicity'  (clinical significance/evidence)
      - 'literature'     (papers about a variant/disease/drug)
      - 'trials'         (active/completed trials related to topic)
    Keep it LLM-based to avoid brittle rules.
    """
    ALLOWED = ["frequency", "pathogenicity", "literature", "trials"]

    def classify(self, text: str) -> str:
        sys = ("Classify the request as exactly one of: "
               "frequency, pathogenicity, literature, trials. Reply with the single word only.")
        out = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "system", "content": sys},
                      {"role": "user", "content": text}],
            temperature=0
        )
        label = (out.choices[0].message.content or "").strip().lower()
        return label if label in self.ALLOWED else "literature"

# ----------------------------
# Specialist Agents
# ----------------------------

@dataclass
class FrequencyAgent:
    """
    Answers allele-frequency questions for a variant.
    Uses ClinVar + NCBI Variation Services via query_clinvar().
    Works best for rsIDs (e.g., 'rs429358'). For HGVS, frequency may be unavailable.
    """
    def answer(self, variant: str, language: str = "en") -> str:
        data = query_clinvar(variant)  # dict with allele_frequencies for rsIDs
        if "error" in data:
            return self._draft(f"Variant not found in ClinVar for query '{variant}'.", language)

        freqs = data.get("allele_frequencies") or {}
        if not freqs:
            return self._draft(
                f"No aggregate allele frequency data available for '{variant}'. "
                f"You may try an rsID if you used HGVS.", language
            )

        # Build a concise report
        lines = [f"Allele frequencies for **{variant}**:"]
        for allele, freq in sorted(freqs.items(), key=lambda kv: kv[0]):
            lines.append(f"- {allele}: {freq}")
        return self._draft("\n".join(lines), language)

    def _draft(self, content: str, language: str) -> str:
        # Optional: localize/clean with LLM. Safe fallback: return as-is.
        try:
            out = client.chat.completions.create(
                model="gpt-4o-mini",
                temperature=0.2,
                messages=[
                    {"role": "system", "content": f"Rewrite this clearly in {language}. Keep technical accuracy."},
                    {"role": "user", "content": content}
                ]
            )
            return out.choices[0].message.content.strip()
        except Exception:
            return content


@dataclass
class PathogenicityAgent:
    """
    Answers pathogenicity/clinical significance questions for a variant using ClinVar.
    Emphasizes significance, review status, associated conditions, and PubMed links.
    """
    def answer(self, variant: str, language: str = "en") -> str:
        data = query_clinvar(variant)
        if "error" in data:
            return f"Variant not found in ClinVar for query '{variant}'."

        cs = data.get("clinical_significance") or "Not provided"
        review = data.get("review_status") or "Not provided"
        gene = data.get("gene") or "N/A"
        conditions = data.get("conditions") or []
        pmids = data.get("pubmed_pmids") or []
        hgvs = data.get("hgvs") or []

        bullet = [
            f"**Variant:** {variant}",
            f"**Gene:** {gene}",
            f"**Clinical Significance:** {cs}",
            f"**Review Status:** {review}",
            f"**Conditions:** {', '.join(conditions) if conditions else '—'}",
            f"**HGVS expressions:** {', '.join(hgvs[:5])}{' …' if len(hgvs) > 5 else ''}",
            f"**Linked PubMed PMIDs:** {', '.join(pmids[:10])}{' …' if len(pmids) > 10 else ''}",
        ]
        return "\n".join(bullet)


@dataclass
class LiteratureAgent:
    """
    Searches PubMed for relevant articles. Returns a compact bibliography.
    """
    max_results: int = 10

    def answer(self, query: str, language: str = "en") -> str:
        results = query_pubmed(query, max_results=self.max_results, include_mesh=True, include_citations=False)
        if not results:
            return f"No PubMed results for '{query}'. Try broader terms or synonyms."

        lines = [f"Top PubMed results for **{query}**:"]
        for i, art in enumerate(results, 1):
            title = art.get("title") or "Untitled"
            year = art.get("year") or "n.d."
            journal = art.get("journal") or ""
            doi = art.get("doi")
            pmid = art.get("pmid")
            lines.append(f"{i}. {title} ({year}). *{journal}*. PMID: {pmid}{'; DOI: ' + doi if doi else ''}")
        return "\n".join(lines)


@dataclass
class TrialsAgent:
    """
    Searches ClinicalTrials.gov v2 for trials matching a term (disease/variant/drug).
    """
    max_results: int = 20

    def answer(self, expr: str, language: str = "en") -> str:
        fields = [
            "NCTId", "BriefTitle", "OverallStatus", "StudyType", "Condition",
            "InterventionName", "Phase", "EnrollmentCount", "LocationCountry",
            "StartDate", "CompletionDate", "LeadSponsorName"
        ]
        out = query_clinicaltrials(expr, max_results=self.max_results, fields=fields)
        studies = out.get("studies", []) if isinstance(out, dict) else []
        if not studies:
            return f"No ClinicalTrials.gov v2 studies found for '{expr}'."

        lines = [f"ClinicalTrials.gov results for **{expr}** (showing up to {self.max_results}):"]
        for s in studies[: self.max_results]:
            title = s.get("BriefTitle") or "Untitled"
            nct = s.get("NCTId") or "NCT?"
            status = s.get("OverallStatus") or "Unknown"
            phase = s.get("Phase") or "—"
            count = s.get("EnrollmentCount")
            lines.append(f"- {nct}: {title} — {status}, Phase: {phase}, Enrollment: {count if count is not None else '—'}")
        return "\n".join(lines)

# ----------------------------
# Orchestrator
# ----------------------------

class BioOrchestrator:
    """
    Central manager that detects language, classifies intent, and routes to the right specialist.
    Add new agents without modifying existing ones (open-closed principle).
    """
    def __init__(self):
        self.language_agent = LanguageAgent()
        self.intent_agent = IntentAgent()
        self.frequency_agent = FrequencyAgent()
        self.pathogenicity_agent = PathogenicityAgent()
        self.literature_agent = LiteratureAgent()
        self.trials_agent = TrialsAgent()

    def handle_query(self, query: str) -> str:
        # Step 1: language detection
        language = self.language_agent.detect_language(query)

        # Step 2: intent classification
        intent = self.intent_agent.classify(query)

        # Step 3: delegate based on intent
        if intent == "frequency":
            # Expect a variant token (rsID or HGVS) somewhere in the text.
            variant = self._extract_variant(query)
            return self.frequency_agent.answer(variant or query, language)
        elif intent == "pathogenicity":
            variant = self._extract_variant(query)
            return self.pathogenicity_agent.answer(variant or query, language)
        elif intent == "trials":
            # Use full query as search expression (disease/drug/variant acceptable).
            return self.trials_agent.answer(query, language)
        else:
            # Default: literature
            return self.literature_agent.answer(query, language)

    # Minimal heuristic; in production you might do better NER or regexp sets
    def _extract_variant(self, text: str) -> Optional[str]:
        import re
        # rsID pattern
        m = re.search(r"\brs\d+\b", text, flags=re.IGNORECASE)
        if m:
            return m.group(0)
        # Simple HGVS token grab (very permissive)
        m = re.search(r"\b(?:[A-Z]{1,3}\d+[A-Z*]|[Cc]\.\d+[A-Z]?>[A-Z*]|[Gg]\.\d+[ACGT]?>[ACGT*])\b", text)
        return m.group(0) if m else None


# ----------------------------
# Example usage (CLI-style)
# ----------------------------
if __name__ == "__main__":
    bot = BioOrchestrator()
    print("Variant & Evidence Concierge — type a question (Ctrl+C to exit)")
    print("Example questions: e.g.")
    print("""
(frequency): How common is the rs429358 variant across different populations?
(pathogenicity): Is the MTHFR variant rs1801133 considered clinically pathogenic, and for which conditions?
(literature): Could you pull recent PubMed papers on lung cancer?
(trials): What clinical trials are currently recruiting for type 1 diabetes?
    """)
    try:
        while True:
            q = input("\nYou: ").strip()
            if not q:
                continue
            print("Processing your query: "+q)
            ans = bot.handle_query(q)
            print("\nAssistant:\n" + ans)
    except KeyboardInterrupt:
        print("\nGoodbye!")

