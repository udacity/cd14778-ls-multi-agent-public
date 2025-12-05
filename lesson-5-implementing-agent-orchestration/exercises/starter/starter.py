# starter_drug_target_workflow.py

import os
import re
import requests
import tempfile
from pathlib import Path
from typing import TypedDict, Dict, Any
from dotenv import load_dotenv

from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.graph import StateGraph, END

# Import the life science tools from the provided helper files
from ls_action_space.vina_helpers_bio import blastp_from_seq, smiles_to_sdf_3d, meeko_ligand, meeko_receptor, vina_dock
from ls_action_space.action_space import query_pubmed

# --- 1. SETUP (PROVIDED) ---
# ---------------------------
load_dotenv(dotenv_path="../../.env")
openai_api_key = os.getenv("OPENAI_API_KEY")
llm = ChatOpenAI(model="gpt-4o-mini", api_key=openai_api_key, base_url="https://openai.vocareum.com/v1", temperature=0.0)

def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    response.raise_for_status()
    pdb_file_path = output_dir / f"{pdb_id}.pdb"
    pdb_file_path.write_text(response.text)
    return pdb_file_path

def parse_sacc_for_pdb_id(sacc: str) -> str:
    match = re.search(r"pdb\|(\w+)\|", sacc)
    if not match: raise ValueError(f"Could not parse PDB ID from: {sacc}")
    return match.group(1)

# --- 2. SPECIALIST TOOLS (PROVIDED) ---
# --------------------------------------
@tool
def run_homology_search(sequence: str) -> Dict[str, Any]:
    """Finds a structural template for a protein sequence using BLAST."""
    print("--- 🔬 Running Homology Search (BLAST) ---")
    hits = blastp_from_seq(sequence, db="pdbaa", max_targets=1)
    if not hits: return {"error": "No BLAST hits found."}
    top_hit = hits[0]
    pdb_id = parse_sacc_for_pdb_id(top_hit["sacc"])
    identity = float(top_hit["pident"])
    print(f"--- ✅ Top BLAST Hit: PDB ID {pdb_id} with {identity}% identity ---")
    return {"pdb_id": pdb_id, "identity": identity}

@tool
def run_docking_simulation(smiles: str, pdb_id: str) -> Dict[str, Any]:
    """Performs computational docking using AutoDock Vina."""
    print(f"--- ⚙️ Starting Docking Simulation for PDB ID {pdb_id} ---")
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        smiles_to_sdf_3d(smiles, temp_path / "ligand.sdf")
        meeko_ligand(temp_path / "ligand.sdf", temp_path / "ligand.pdbqt")
        receptor_pdb = download_pdb(pdb_id, temp_path)
        meeko_receptor(receptor_pdb, temp_path / "receptor.pdbqt")
        vina_log = vina_dock(
            receptor_pdbqt=temp_path / "receptor.pdbqt",
            ligand_pdbqt=temp_path / "ligand.pdbqt",
            out_pdbqt=temp_path / "poses.pdbqt",
            center=(15, 54, 16), size=(25, 25, 25)
        )
        match = re.search(r"^\s+1\s+(-?\d+\.\d+)", vina_log, re.MULTILINE)
        if not match: return {"error": "Could not parse docking score."}
        best_score = float(match.group(1))
        print(f"--- ✅ Docking Complete. Best Affinity: {best_score} kcal/mol ---")
        return {"binding_affinity_kcal_mol": best_score}

@tool
def run_literature_lookup(molecule_name: str, target_name: str) -> Dict[str, Any]:
    """Searches PubMed for existing studies."""
    print(f"--- 📚 Searching PubMed for {molecule_name} and {target_name} ---")
    query = f'"{molecule_name}" AND "{target_name}"'
    results = query_pubmed(query)
    if not results: return {"summary": "No relevant publications found."}
    summary = f"Found {len(results)} relevant publications. Top result: '{results[0]['title']}'."
    return {"summary": summary, "publications": results}

# --- 3. ORCHESTRATOR GRAPH DEFINITION (STUDENT'S TASK) ---
# --------------------------------------------------------

### YOUR CODE HERE - STEP 1: Define the state ###
# Define the fields that will be passed between nodes in the workflow.
# Think about what information you need to start with, what each tool
# will generate, and what you need for the final report.
class GraphState(TypedDict):
    """Defines the state that flows through the orchestrator graph."""
    pass # Replace with your state fields

### YOUR CODE HERE - STEP 2: Implement the nodes ###
# Each node function takes the current state and calls a tool.
# It should return a dictionary with the fields to update in the state.

def homology_search_node(state: GraphState):
    """Node that runs the homology search tool."""
    pass

def docking_node(state: GraphState):
    """Node that runs the docking simulation tool."""
    pass

def literature_lookup_node(state: GraphState):
    """Node that runs the literature lookup tool."""
    pass

def generate_report_node(state: GraphState):
    """Node that uses the LLM to generate a final summary report."""
    print("--- 📝 Generating Final Report ---")
    prompt = f"""
    Generate a concise summary report for the drug-target interaction analysis.
    Molecule: {state['molecule_name']}, Target: {state['target_name']}
    - Homology Search: Found template {state.get('template_pdb_id')} with {state.get('sequence_identity', 0):.2f}% identity.
    - Docking Simulation: {'Not performed (sequence identity < 90%).' if not state.get('docking_score') else f'Predicted binding affinity of {state.get('docking_score')} kcal/mol.'}
    - Literature Review: {state.get('literature_summary')}
    Conclude with a brief, one-sentence overall assessment.
    """
    response = llm.invoke(prompt)
    return {"final_report": response.content}

### YOUR CODE HERE - STEP 3: Implement the conditional edge ###
def should_dock(state: GraphState) -> str:
    """
    This is the conditional logic. Check the sequence identity in the state.
    - If identity is >= 90, return the string "docking".
    - Otherwise, return the string "literature_lookup" to skip the docking step.
    """
    print(f"--- 🤔 Decision Point: Sequence identity is {state.get('sequence_identity', 0)}%. Threshold is 90%. ---")
    pass # Implement the conditional logic here

# --- 4. ASSEMBLE THE GRAPH (STUDENT'S TASK) ---
# ---------------------------------------------
workflow = StateGraph(GraphState)

### YOUR CODE HERE - STEP 4: Connect the nodes and edges ###
# Use the following methods to build the graph:
# - workflow.add_node("node_name", node_function)
# - workflow.set_entry_point("start_node_name")
# - workflow.add_edge("from_node", "to_node")
# - workflow.add_conditional_edges("source_node", conditional_function, {"path_name": "destination_node", ...})


# Compile the graph into a runnable application
app = workflow.compile()

# --- 5. EXECUTION (PROVIDED) ---
# -------------------------------
if __name__ == "__main__":
    print("\n--- Starting Automated Drug-Target Interaction Workflow ---\n")
    inputs = {
        "molecule_name": "Aspirin",
        "molecule_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "target_name": "COX-1",
        "target_sequence": (
            "MSRSLLLWFLFLLLLPPLPVLLADPGAPTPVNPCCYYPCQHQGICVRFGLDRYQCDCTRTGYSGPNCTIPGLWTWLR"
            "NSLRPSPSFTHFLLTHGRWFWEFVNATFIREMIMLHPDSAWLREAGHPPLFEGNPQYCRFLQDRQDDLRGYPGRSGLG"
            "HCNFTSFSYEQLTDQSMKFYQTLYWHFLRFSLPLPSVEDYNLFARGALADVLCTYKIWPLPESTFGGEVGFNIVKTAT"
            "LKKLVCLNTKTCPYVSFRVPDASQDDGPAVERPSTEL"
        )
    }
    final_state = app.invoke(inputs)
    print("\n" + "="*50 + "\n          FINAL WORKFLOW REPORT\n" + "="*50)
    print(final_state.get('final_report', 'Workflow did not produce a report.'))
    print("="*50 + "\n")
