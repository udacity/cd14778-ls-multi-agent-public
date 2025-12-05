#
# 1. Imports and Setup (Provided)
#
import os
import json
import random

# --- LangChain & Environment ---
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.tools import tool
from langchain.agents import AgentExecutor, create_openai_tools_agent
from langchain_openai import ChatOpenAI
from dotenv import load_dotenv

# --- Bioinformatics Tools ---
# Make sure BLAST+ is installed and the BLASTDB env var is set
from ls_action_space.vina_helpers_bio import blastp_from_seq

#
# 2. Environment and Model Configuration (Provided)
#
load_dotenv("../../.env")
if not os.getenv("OPENAI_API_KEY"):
    raise ValueError("OPENAI_API_KEY not found in .env file.")

llm = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=os.getenv("OPENAI_API_KEY"),
    base_url="https://openai.vocareum.com/v1",
    temperature=0,
)

#
# 3. Agent Roles Defined as Tools
#
# This dictionary is the "shared whiteboard" for our tools to communicate.
PIPELINE_STATE = {}

# --- Tool 1: FASTA Validator (Provided as a complete example) ---
@tool
def validate_and_normalize_fasta(fasta_string: str) -> str:
    """
    Role 1: FASTA Validator.
    Validates a FASTA string. If valid, it normalizes the sequence, stores it
    in the shared pipeline state under a unique header, and returns the header.
    If invalid, it returns a descriptive error message. This is the mandatory first step.
    """
    lines = fasta_string.strip().split('\n')
    if not lines or not lines[0].startswith('>'):
        return "Error: Invalid FASTA format. Must start with a '>' header."

    header = lines[0][1:].strip()
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    valid_chars = "ACDEFGHIKLMNPQRSTVWY"
    if not sequence or not all(c in valid_chars for c in sequence):
        return f"Error: Sequence for '{header}' is empty or contains invalid characters."

    PIPELINE_STATE[header] = {"status": "validated", "sequence": sequence}
    return f"Successfully validated sequence with header: {header}"

# --- Tool 2: BLAST Specialist (Student's Task) ---
@tool
def run_homology_search(header: str) -> str:
    """
    YOUR CODE HERE: Task 1
    Define this tool's role in one sentence. What is its job?
    This docstring is what the agent's "brain" (the LLM) reads to understand the tool.
    """
    if header not in PIPELINE_STATE:
        return f"Error: Cannot find sequence for header '{header}'. Please validate it first."
    
    sequence = PIPELINE_STATE[header]["sequence"]
    
    try:
        # YOUR CODE HERE: Task 2
        # Call the imported `blastp_from_seq` function with the sequence.
        # Store the results and update the status in the PIPELINE_STATE dictionary.
        # Return a success message string, e.g., "Homology search completed for {header}."
        print(f"  [LIVE TOOL RUNNING]: blastp for '{header}'...")
        blast_hits = # ... your call to blastp_from_seq ...
        
        PIPELINE_STATE[header]["blast_hits"] = blast_hits
        PIPELINE_STATE[header]["status"] = "homology_searched"
        return f"Homology search completed for {header}. Found {len(blast_hits)} hits."

    except Exception as e:
        return f"An error occurred during BLAST search: {str(e)}"

# --- Tool 3: Docking Simulator (Student's Task) ---
@tool
def simulate_docking_and_get_score(header: str) -> str:
    """
    YOUR CODE HERE: Task 3
    Define this tool's role in one sentence. What is its job?
    """
    if header not in PIPELINE_STATE:
        return f"Error: Cannot run docking for unknown sequence '{header}'. Please validate first."
    
    # YOUR CODE HERE: Task 4
    # 1. Generate a random docking score (a float, e.g., between -12.0 and -5.0).
    # 2. Create a result dictionary, e.g., {"target": "P38-MAPK", "docking_score_kcal_mol": your_score}.
    # 3. Store this dictionary and update the status in PIPELINE_STATE.
    # 4. Return a success message string with the score, e.g., "Docking simulation complete. Score: {score}"
    
    print(f"  [SIMULATION RUNNING]: Generating docking score for '{header}'...")
    mock_score = round(random.uniform(-12.0, -5.0), 2)
    result = {"target": "P38-MAPK", "docking_score_kcal_mol": mock_score}
    
    PIPELINE_STATE[header]["docking_result"] = result
    PIPELINE_STATE[header]["status"] = "docking_complete"
    return f"Docking simulation for {header} complete. Affinity Score: {mock_score} kcal/mol."


#
# 4. The Coordinator Agent Architecture
#
class PipelineCoordinator:
    """A class to encapsulate the agent's architecture."""
    def __init__(self, llm_model):
        self.tools = [
            validate_and_normalize_fasta,
            run_homology_search,
            simulate_docking_and_get_score,
        ]

        # YOUR CODE HERE: Task 5
        # This is the agent's "job description". Write a clear, step-by-step
        # set of instructions telling the agent HOW and in WHAT ORDER to use its tools.
        # Be explicit about starting with validation.
        prompt_template = """
        You are an automated bioinformatics pipeline coordinator. Your responsibility is to
        process a protein sequence using tools in a specific, non-negotiable order.

        For every sequence, you MUST execute the following steps sequentially:
        1.  Always begin by calling `validate_and_normalize_fasta`.
        2.  If validation succeeds, use the header to call `run_homology_search`.
        3.  After the search, use the same header to call `simulate_docking_and_get_score`.

        Report the final result of the last step you performed.
        """
        
        # --- The rest of the setup is provided ---
        self.prompt = ChatPromptTemplate.from_messages([
            ("system", prompt_template), ("human", "{input}"), ("placeholder", "{agent_scratchpad}"),
        ])
        agent = create_openai_tools_agent(llm_model, self.tools, self.prompt)
        self.agent_executor = AgentExecutor(agent=agent, tools=self.tools, verbose=True)

    def run_pipeline(self, fasta_string: str):
        """Invokes the agent to run the full pipeline."""
        print(f"\n{'='*60}\n🚀 Coordinating pipeline for: {fasta_string.splitlines()[0]}...\n{'='*60}")
        self.agent_executor.invoke({"input": fasta_string})

#
# 5. Simulation Runner (Provided)
#
def main():
    """Initializes the coordinator and runs the simulation."""
    coordinator = PipelineCoordinator(llm_model=llm)

    candidate_1_fasta = """>candidate_001|KRAS_G12C_variant
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS
"""
    coordinator.run_pipeline(candidate_1_fasta)

    print("\n\n" + "="*60)
    print("📋 Final Pipeline State:")
    print(json.dumps(PIPELINE_STATE, indent=2))
    print("="*60)

if __name__ == "__main__":
    main()
