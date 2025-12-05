import os
import re
import requests
import tempfile
from pathlib import Path
from typing import TypedDict, Annotated, List, Dict, Any
from dotenv import load_dotenv

from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.graph import StateGraph, END

# Import the life science tools from the provided helper files
from ls_action_space.vina_helpers_bio import blastp_from_seq, smiles_to_sdf_3d, meeko_ligand, meeko_receptor, vina_dock
from ls_action_space.action_space import query_pubmed


# --- 1. SETUP: LOAD ENVIRONMENT AND INITIALIZE LLM ---
# ---------------------------------------------------

# Load the API key from the .env file
load_dotenv(dotenv_path="../../.env")
openai_api_key = os.getenv("OPENAI_API_KEY")

# Initialize the OpenAI model to be used by the orchestrator
# This model will be used for generating the final summary report.
llm = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=openai_api_key,
    base_url="https://openai.vocareum.com/v1",
    temperature=0.0
)

# --- 2. HELPER FUNCTIONS ---
# ---------------------------

def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Downloads a PDB file from the RCSB database."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    response.raise_for_status()  # Raise an exception for bad status codes
    pdb_file_path = output_dir / f"{pdb_id}.pdb"
    pdb_file_path.write_text(response.text)
    print(f"--- ✅ Downloaded PDB file: {pdb_file_path} ---")
    return pdb_file_path

def parse_sacc_for_pdb_id(sacc: str) -> str:
    """Parses a BLAST subject accession string (e.g., 'pdb|1EQG|A') to get the PDB ID."""
    match = re.search(r"pdb\|(\w+)\|", sacc)
    if not match:
        raise ValueError(f"Could not parse PDB ID from BLAST subject accession: {sacc}")
    return match.group(1)

# --- 3. SPECIALIST TOOLS DEFINITION ---
# --------------------------------------
# Each function is decorated with @tool, turning it into a LangChain Tool.
# These tools act as the "specialist agents" that the orchestrator can call.

@tool
def run_homology_search(sequence: str) -> Dict[str, Any]:
    """
    Finds a structural template for a protein sequence using BLAST.
    
    **PREREQUISITE**: Requires BLAST+ to be installed and the 'pdbaa' database
    to be available in the $BLASTDB environment directory.
    """
    print("--- 🔬 Running Homology Search (BLAST) ---")
    try:
        # Run BLAST against the Protein Data Bank database
        hits = blastp_from_seq(sequence, db="pdbaa", max_targets=1)
        if not hits:
            return {"error": "No significant BLAST hits found."}
        
        top_hit = hits[0]
        pdb_id = parse_sacc_for_pdb_id(top_hit["sacc"])
        identity = float(top_hit["pident"])
        
        print(f"--- ✅ Top BLAST Hit: PDB ID {pdb_id} with {identity}% identity ---")
        return {"pdb_id": pdb_id, "identity": identity}
    except Exception as e:
        return {"error": f"BLAST search failed: {e}"}

@tool
def run_docking_simulation(smiles: str, pdb_id: str) -> Dict[str, Any]:
    """
    Performs computational docking of a small molecule to a protein target using AutoDock Vina.
    
    **PREREQUISITES**: Requires Meeko and AutoDock Vina command-line tools to be installed.
    """
    print(f"--- ⚙️ Starting Docking Simulation for SMILES on PDB ID {pdb_id} ---")
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # 1. Prepare Ligand (Small Molecule)
            smiles_to_sdf_3d(smiles, temp_path / "ligand.sdf")
            meeko_ligand(temp_path / "ligand.sdf", temp_path / "ligand.pdbqt")

            # 2. Prepare Receptor (Protein Target)
            receptor_pdb = download_pdb(pdb_id, temp_path)
            meeko_receptor(receptor_pdb, temp_path / "receptor.pdbqt")

            # 3. Run Docking
            print("--- 🚀 Running AutoDock Vina... (This may take a minute) ---")
            vina_log = vina_dock(
                receptor_pdbqt=temp_path / "receptor.pdbqt",
                ligand_pdbqt=temp_path / "ligand.pdbqt",
                out_pdbqt=temp_path / "poses.pdbqt",
                center=(15, 54, 16), # Example coordinates, should be determined from protein structure
                size=(25, 25, 25)
            )
            
            # 4. Parse result
            match = re.search(r"^\s+1\s+(-?\d+\.\d+)", vina_log, re.MULTILINE)
            if not match:
                return {"error": "Could not parse docking score from Vina output."}
                
            best_score = float(match.group(1))
            print(f"--- ✅ Docking Complete. Best Affinity: {best_score} kcal/mol ---")
            return {"binding_affinity_kcal_mol": best_score}
            
    except Exception as e:
        return {"error": f"Docking simulation failed: {e}"}

@tool
def run_literature_lookup(molecule_name: str, target_name: str) -> Dict[str, Any]:
    """Searches PubMed for existing studies on the molecule-target interaction."""
    query = f'"{molecule_name}" AND "{target_name}"'
    results = query_pubmed(query)
    if not results:
        return {"summary": "No relevant publications found in PubMed."}
    
    summary = f"Found {len(results)} relevant publications. Top result: '{results[0]['title']}' by {results[0]['authors'][0]} et al. ({results[0]['year']})."
    return {"summary": summary, "publications": results}

# --- 4. ORCHESTRATOR GRAPH DEFINITION (using LangGraph) ---
# ---------------------------------------------------------

class GraphState(TypedDict):
    """Defines the state that flows through the orchestrator graph."""
    molecule_smiles: str
    molecule_name: str
    target_sequence: str
    target_name: str
    
    # Fields to be populated by the tools
    template_pdb_id: str
    sequence_identity: float
    docking_score: float
    literature_summary: str
    final_report: str
    error_message: str

# Define the Nodes of the Graph
def homology_search_node(state: GraphState):
    """Node that runs the homology search tool."""
    result = run_homology_search(state["target_sequence"])
    if "error" in result:
        return {"error_message": result["error"]}
    return {"template_pdb_id": result["pdb_id"], "sequence_identity": result["identity"]}

def docking_node(state: GraphState):
    """Node that runs the docking simulation tool."""
    result = run_docking_simulation(smiles=state["molecule_smiles"], pdb_id=state["template_pdb_id"])
    if "error" in result:
        return {"error_message": result["error"]}
    return {"docking_score": result.get("binding_affinity_kcal_mol")}

def literature_lookup_node(state: GraphState):
    """Node that runs the literature lookup tool."""
    result = run_literature_lookup(molecule_name=state["molecule_name"], target_name=state["target_name"])
    return {"literature_summary": result["summary"]}

def generate_report_node(state: GraphState):
    """Node that uses the LLM to generate a final summary report."""
    print("--- 📝 Generating Final Report ---")
    if state.get("error_message"):
        report_text = f"Workflow failed with an error: {state['error_message']}"
    else:
        prompt = f"""
        Generate a concise summary report for the drug-target interaction analysis.
        
        Molecule: {state['molecule_name']}
        Target: {state['target_name']}
        
        Analysis Results:
        - Homology Search: Found template PDB {state['template_pdb_id']} with {state['sequence_identity']:.2f}% identity.
        - Docking Simulation: {'Not performed (sequence identity < 90%).' if not state.get('docking_score') else f'Predicted binding affinity of {state["docking_score"]} kcal/mol.'}
        - Literature Review: {state['literature_summary']}
        
        Conclude with a brief, one-sentence overall assessment.
        """
        response = llm.invoke(prompt)
        report_text = response.content
        
    return {"final_report": report_text}

# Define the Conditional Edge
def should_dock(state: GraphState) -> str:
    """
    This is the conditional logic branch. It decides the next step based on sequence identity.
    If identity is >= 90%, it proceeds to docking. Otherwise, it skips to literature lookup.
    """
    if state.get("error_message"):
        return "generate_report" # Go straight to report on error
        
    identity = state.get("sequence_identity", 0)
    print(f"--- 🤔 Decision Point: Sequence identity is {identity}%. Threshold is 90%. ---")
    if identity >= 90:
        print("--- ✅ Identity >= 90%. Proceeding to Docking. ---")
        return "docking"
    else:
        print("--- ❌ Identity < 90%. Skipping Docking. ---")
        return "literature_lookup"

# Assemble the Graph
workflow = StateGraph(GraphState)
workflow.add_node("homology_search", homology_search_node)
workflow.add_node("docking", docking_node)
workflow.add_node("literature_lookup", literature_lookup_node)
workflow.add_node("generate_report", generate_report_node)

workflow.set_entry_point("homology_search")
workflow.add_conditional_edges(
    "homology_search",
    should_dock,
    {"docking": "docking", "literature_lookup": "literature_lookup", "generate_report": "generate_report"}
)
workflow.add_edge("docking", "literature_lookup")
workflow.add_edge("literature_lookup", "generate_report")
workflow.add_edge("generate_report", END)

# Compile the graph into a runnable application
app = workflow.compile()


# --- 5. EXECUTION ---
# --------------------

if __name__ == "__main__":
    print("\n--- Starting Automated Drug-Target Interaction Workflow ---\n")
    
    # --- INPUTS ---
    # Aspirin targeting Cyclooxygenase-1 (COX-1)
    # Using a sequence with >90% identity to a known PDB structure (1EQG)
    inputs = {
        "molecule_name": "Aspirin",
        "molecule_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", # SMILES for Aspirin
        "target_name": "COX-1",
        "target_sequence": (
            "MSRSLLLWFLFLLLLPPLPVLLADPGAPTPVNPCCYYPCQHQGICVRFGLDRYQCDCTRTGYSGPNCTIPGLWTWLR"
            "NSLRPSPSFTHFLLTHGRWFWEFVNATFIREMIMLHPDSAWLREAGHPPLFEGNPQYCRFLQDRQDDLRGYPGRSGLG"
            "HCNFTSFSYEQLTDQSMKFYQTLYWHFLRFSLPLPSVEDYNLFARGALADVLCTYKIWPLPESTFGGEVGFNIVKTAT"
            "LKKLVCLNTKTCPYVSFRVPDASQDDGPAVERPSTEL"
        )
    }

    # Invoke the workflow with the inputs
    final_state = app.invoke(inputs)

    # Print the final report
    print("\n" + "="*50)
    print("          FINAL WORKFLOW REPORT")
    print("="*50)
    print(final_state['final_report'])
    print("="*50 + "\n")
