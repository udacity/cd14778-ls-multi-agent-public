# main.py

import argparse
import os
import random
from datetime import datetime
import report

from logging_utils import JSONLLogger

import sys
sys.path.append("/workspace/code") # set this to the path containing the "ls_action_space" package

# tools for querying ClinVar, PubMed, ClinTrials
from ls_action_space.action_space import query_clinvar, query_pubmed, query_clinicaltrials

def demonstrate_tool_usage():
    print("\n" + "=" * 50)
    print("DEMONSTRATING TOOL USAGE (from tools.py)")
    print("=" * 50)

    # 1. Example: query_clinvar
    print("\n🔍 1. Calling query_clinvar('rs113993960')...")
    clinvar_result = query_clinvar('rs113993960')
    print("ClinVar Result:")
    print(json.dumps(clinvar_result, indent=2))

    # 2. Example: query_pubmed
    print("\n🔍 2. Calling query_pubmed('CFTR cystic fibrosis')...")
    pubmed_result = query_pubmed('CFTR cystic fibrosis')
    print("PubMed Result:")
    print(json.dumps(pubmed_result, indent=2))

    # 3. Example: query_clinicaltrials
    print("\n🔍 3. Calling query_clinicaltrials('CFTR cystic fibrosis')...")
    trials_result = query_clinicaltrials('CFTR cystic fibrosis')
    print("ClinicalTrials.gov Result:")
    print(json.dumps(trials_result, indent=2))
    print("=" * 50)
    print("END OF DEMONSTRATION")
    print("=" * 50 + "\n")


class MockOrchestrator:
    """
    A mock orchestrator to generate dummy data for the report.
    YOUR TASK: Replace this entire class with your multi-agent implementation.
    Your agents should produce data in the same format and use the logger as shown below.
    """

    def __init__(self, logger, seed):
        self.logger = logger
        self.seed = seed
        print(f"MockOrchestrator initialized with seed: {self.seed}")

    def run_workflow(self, variants_path, phenotypes_path, run_dir):
        """
        Simulates the full agent workflow and produces hardcoded results.
        """
        run_id = os.path.basename(run_dir)
        self.logger.log_orchestrator_start(run_id, {"variants": 10, "phenotypes": 3})

        print("MockOrchestrator: Running mock workflow...")

        # --- YOUR AGENT LOGIC WILL REPLACE THIS SECTION ---

        # 1. Mock Variant Prioritization
        self.logger.log_action(
            agent="VariantPrioritizer",
            tool="query_clinvar",
            target="rs113993960",
            inputs={"variant": "rs113993960"},
            outputs={"clinical_significance": "Pathogenic", "gene": "CFTR"},
            latency_ms=150
        )
        time.sleep(0.1)  # Simulate work
        top_variants = [
            {
                'gene': 'CFTR', 'variant': 'p.Phe508del', 'score': 4,
                'components': {'sig': 3, 'pheno': 1, 'rarity': 0},
                'rationale': "ClinVar pathogenic; matches 'cystic fibrosis'; rare allele (<1%)"
            },
            {
                'gene': 'TTR', 'variant': 'p.Val142Ile', 'score': 4,
                'components': {'sig': 2, 'pheno': 2, 'rarity': 0},
                'rationale': "ClinVar likely pathogenic; matches 'amyloidosis'; rare allele (<1%)"
            }
        ]

        # 2. Mock Evidence Synthesis
        selected_pmids = ["20301428", "20301611"]
        pmid_syntheses = {
            "20301428": "Study linking CFTR gene variants to cystic fibrosis.",
            "20301611": "Research on TTR gene function and its role in hereditary amyloidosis."
        }
        self.logger.log_action(
            agent="EvidenceSynthesizer",
            tool="pmid_selection",
            target="global",
            inputs={"top_variants_count": 2},
            outputs={"selected_pmids": selected_pmids},
            latency_ms=50
        )
        time.sleep(0.1)  # Simulate work

        # 3. Mock Trial Matching
        trial_match = {
            'nct_id': 'NCT04858883',
            'title': 'A Study of an Investigational Drug in People With Cystic Fibrosis',
            'status': 'Recruiting',
            'reason': "Mentions cystic fibrosis (RECRUITING, ['PHASE2'])"
        }
        self.logger.log_action(
            agent="TrialMatcher",
            tool="query_clinicaltrials",
            target="CFTR cystic fibrosis",
            inputs={"expr": "CFTR cystic fibrosis"},
            outputs={"matched_nct": trial_match['nct_id']},
            latency_ms=250
        )

        # --- END OF MOCK DATA SECTION ---

        # 4. Generate and write the final report artifacts
        report.write_artifacts(
            run_dir=run_dir, run_id=run_id, seed=self.seed,
            top_variants=top_variants, selected_pmids=selected_pmids,
            pmid_syntheses=pmid_syntheses, trial_match=trial_match
        )

        self.logger.log_orchestrator_finish(run_id, "completed")
        print(f"SUCCESS: Mock artifacts written to {run_dir}")


def main():
    """Main CLI entrypoint for Orphan Finder."""
    parser = argparse.ArgumentParser(description="Orphan Finder - A mock variant-to-therapy pipeline.")
    parser.add_argument("--variants", default="variants.csv", help="Path to variants CSV file")
    parser.add_argument("--phenotypes", default="phenotype_keywords.txt", help="Path to phenotype keywords text file")
    parser.add_argument("--outdir", default="runs", help="Output directory for run results")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()

    # 1. Set up run environment
    run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_id = f"{run_timestamp}_{random.randint(1000, 9999)}"
    run_dir = os.path.join(args.outdir, run_id)
    os.makedirs(run_dir, exist_ok=True)

    print(f"Starting run {run_id}")
    print(f"Output will be saved to: {run_dir}")

    # 2. Initialize logger
    log_filepath = os.path.join(run_dir, "run.jsonl")
    logger = JSONLLogger(log_filepath, seed=args.seed)

    # 3. Initialize and run the orchestrator
    # NOTE FOR STUDENTS: Replace MockOrchestrator with your real Orchestrator class.
    orchestrator = MockOrchestrator(logger=logger, seed=args.seed)
    orchestrator.run_workflow(
        variants_path=args.variants,
        phenotypes_path=args.phenotypes,
        run_dir=run_dir
    )


if __name__ == "__main__":
    main()
