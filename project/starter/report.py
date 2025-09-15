# report.py

import os
from datetime import datetime
from typing import Dict, Any, List
import pandas as pd

def generate_report_md(run_id: str, seed: int, top_variants: List[Dict[str, Any]],
                       selected_pmids: List[str], pmid_syntheses: Dict[str, str],
                       trial_match: Dict[str, Any]) -> str:
    """Generates the full content for the report.md file."""
    timestamp = datetime.now().isoformat(timespec='seconds')

    # --- Header ---
    lines = [
        "# Orphan Finder Brief",
        f"_Run: {run_id} · Seed: {seed} · Date: {timestamp}_",
        "",
        f"## Top Variants (n={len(top_variants)})",
        "| Gene | Variant | Score | Components (sig, pheno, rarity) | Rationale |",
        "|---|---|---:|---|---|"
    ]

    # --- Top Variants Table ---
    for variant in top_variants:
        components = variant.get("components", {})
        comp_text = f"{components.get('sig', 0)} / {components.get('pheno', 0)} / {components.get('rarity', 0)}"
        lines.append(
            f"| {variant.get('gene', 'N/A')} "
            f"| {variant.get('variant', 'N/A')} "
            f"| {variant.get('score', 0)} "
            f"| {comp_text} "
            f"| {variant.get('rationale', '')} |"
        )

    # --- Evidence Section ---
    lines.extend(["", "## Evidence (PMIDs ≤2)"])
    if selected_pmids:
        for pmid in selected_pmids:
            synthesis = pmid_syntheses.get(pmid, "Synthesis not available.")
            lines.append(f"- **PMID {pmid}** — {synthesis}")
    else:
        lines.append("- No qualifying evidence found.")

    # --- Trial Match Section ---
    lines.extend(["", "## Trial Match"])
    nct_id = trial_match.get("nct_id")
    if nct_id and nct_id != "NCT_NOT_FOUND":
        lines.append(
            f"- **{nct_id}** — {trial_match.get('title', '')} ({trial_match.get('status', '')}). "
            f"**Why:** {trial_match.get('reason', '')}"
        )
    else:
        lines.append("- **None found** — No recruiting trials matched the criteria.")

    # --- Disclaimer ---
    lines.extend([
        "",
        "> **Disclaimer:** Educational exercise. Not for clinical decision-making. Verify all sources."
    ])

    return "\n".join(lines)

def generate_ranked_variants_csv(top_variants: List[Dict[str, Any]]) -> pd.DataFrame:
    """Generates a pandas DataFrame for the ranked_variants.csv file."""
    rows = []
    for variant in top_variants:
        components = variant.get("components", {})
        comp_text = (
            f"sig={components.get('sig', 0)};"
            f"pheno={components.get('pheno', 0)};"
            f"rarity={components.get('rarity', 0)}"
        )
        rows.append({
            "gene": variant.get("gene", ""),
            "variant": variant.get("variant", ""),
            "score": variant.get("score", 0),
            "components": comp_text,
            "rationale": variant.get("rationale", "")
        })
    return pd.DataFrame(rows, columns=["gene", "variant", "score", "components", "rationale"])

def write_artifacts(run_dir: str, run_id: str, seed: int,
                    top_variants: List[Dict[str, Any]],
                    selected_pmids: List[str], pmid_syntheses: Dict[str, str],
                    trial_match: Dict[str, Any]) -> None:
    """
    Writes all output artifacts (report.md, ranked_variants.csv) to the run directory.
    """
    # 1. Generate and write the markdown report
    report_content = generate_report_md(run_id, seed, top_variants, selected_pmids, pmid_syntheses, trial_match)
    report_path = os.path.join(run_dir, "report.md")
    with open(report_path, 'w') as f:
        f.write(report_content)

    # 2. Generate and write the CSV of ranked variants
    df = generate_ranked_variants_csv(top_variants)
    csv_path = os.path.join(run_dir, "ranked_variants.csv")
    df.to_csv(csv_path, index=False)