### **Solution: Drug-Discovery Sequence Screening**

#
# 1. Imports and Setup
#
import os
import json
import random
from typing import Dict, Any, List

# --- LangChain Core Components ---
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.tools import tool
from langchain.agents import AgentExecutor, create_openai_tools_agent

# --- LangChain OpenAI Integration ---
from langchain_openai import ChatOpenAI

# --- Environment and Bioinformatics Tools ---
from dotenv import load_dotenv
# These are the real helper functions that call command-line tools.
from ls_action_space.vina_helpers_bio import blastp_from_seq

#
# 2. Environment and Model Configuration
#
load_dotenv("../../.env")
if not os.getenv("OPENAI_API_KEY"):
    raise ValueError("OPENAI_API_KEY not found in .env file.")
# NOTE: For the BLAST tool to work, you must have NCBI BLAST+ installed
# and the BLASTDB environment variable pointing to your databases.
# Example: export BLASTDB=/path/to/your/blast/db
if not os.getenv("BLASTDB"):
    print("⚠️ WARNING: BLASTDB environment variable not set. The homology search tool may fail.")

# Initialize the LLM for the agent's "brain"
llm = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=os.getenv("OPENAI_API_KEY"),
    base_url="https://openai.vocareum.com/v1",
    temperature=0, # Use a low temperature for predictable, tool-following behavior
)

#
# 3. Agent Roles Defined as Tools
#
# This dictionary serves as a shared "whiteboard" for our agents (tools) to
# read from and write to, allowing them to pass complex data between steps.
PIPELINE_STATE = {}

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
        return "Error: Invalid FASTA format. The input must start with a '>' header line."

    header = lines[0][1:].strip()
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    valid_chars = "ACDEFGHIKLMNPQRSTVWY" # Standard amino acids, no ambiguity codes
    if not sequence or not all(c in valid_chars for c in sequence):
        return f"Error: Sequence for '{header}' is empty or contains invalid characters."

    # Store initial state for other tools to access
    PIPELINE_STATE[header] = {"status": "validated", "sequence": sequence}
    return f"Successfully validated sequence with header: {header}"

@tool
def run_homology_search(header: str) -> str:
    """
    Role 2: BLAST Specialist.
    Runs a real BLASTp sequence similarity search using the blastp command-line tool.
    It requires a sequence header that has already been validated. It fetches the
    sequence from the shared state, runs the search, and stores the results back.
    """
    if header not in PIPELINE_STATE:
        return f"Error: Cannot find sequence for header '{header}'. Please validate it first."
    
    sequence = PIPELINE_STATE[header]["sequence"]
    
    try:
        print(f"  [LIVE TOOL RUNNING]: blastp for '{header}'...")
        # This is a REAL call to the NCBI BLAST+ command-line tool.
        blast_hits = blastp_from_seq(seq=sequence, db="pdbaa", threads=4)
        
        # Storing the structured result in our shared state
        PIPELINE_STATE[header]["blast_hits"] = blast_hits
        PIPELINE_STATE[header]["status"] = "homology_searched"
        return f"Homology search completed for {header}. Found {len(blast_hits)} significant hits."
    except FileNotFoundError:
        return "Error: 'blastp' command not found. Is NCBI BLAST+ installed and in your PATH?"
    except Exception as e:
        # Catch other potential errors from the subprocess call
        return f"An error occurred during BLAST search: {str(e)}"

@tool
def simulate_docking_and_get_score(header: str) -> str:
    """
    Role 3: Docking Simulator.
    Simulates a molecular docking calculation for a given sequence header. As docking
    requires 3D structures not available from a sequence, this tool provides a realistic
    placeholder result. It returns a random docking score and stores it in the state.
    """
    if header not in PIPELINE_STATE:
        return f"Error: Cannot run docking for unknown sequence '{header}'. Please validate first."
    
    print(f"  [SIMULATION RUNNING]: Generating docking score for '{header}'...")
    # More negative scores indicate better binding affinity
    mock_score = round(random.uniform(-12.0, -5.0), 2)
    result = {"target": "P38-MAPK", "docking_score_kcal_mol": mock_score}
    
    # Update the shared state with the simulated result
    PIPELINE_STATE[header]["docking_result"] = result
    PIPELINE_STATE[header]["status"] = "docking_complete"
    return f"Docking simulation for {header} complete. Affinity Score: {mock_score} kcal/mol."

#
# 4. The Coordinator Agent Architecture
#
class PipelineCoordinator:
    """
    This class defines the high-level architecture for our coordinator agent.
    It includes its "job description" (the system prompt) and its "skills" (the tools).
    """
    def __init__(self, llm_model):
        self.tools = [
            validate_and_normalize_fasta,
            run_homology_search,
            simulate_docking_and_get_score,
        ]

        # The System Prompt is the agent's core identity and instruction set.
        # It clearly defines its role, responsibilities, and the workflow it must follow.
        prompt_template = """
        You are an automated bioinformatics pipeline coordinator. Your responsibility is to
        process a protein sequence through a standard screening workflow by calling tools
        in a specific, non-negotiable order.

        For every sequence provided by the user, you MUST execute the following steps sequentially:
        1.  Always begin by calling `validate_and_normalize_fasta`.
        2.  If validation succeeds, use the header returned by the tool to call `run_homology_search`.
        3.  After the homology search is complete, use the same header to call `simulate_docking_and_get_score`.

        If any step returns an error, stop the pipeline for that sequence and report the error.
        After the final step, state the result clearly.
        """
        
        self.prompt = ChatPromptTemplate.from_messages([
            ("system", prompt_template),
            ("human", "{input}"),
            # The 'agent_scratchpad' is where LangChain stores the conversation history
            # between the agent and the tools it calls.
            ("placeholder", "{agent_scratchpad}"),
        ])

        # This constructor binds the LLM, tools, and prompt together into a functional agent.
        agent = create_openai_tools_agent(llm_model, self.tools, self.prompt)

        # The AgentExecutor is the runtime environment that executes the agent's chosen tools.
        self.agent_executor = AgentExecutor(
            agent=agent,
            tools=self.tools,
            verbose=True # Set to True to see the agent's "chain of thought".
        )

    def run_pipeline(self, fasta_string: str):
        """Invokes the agent to run the full pipeline for a given FASTA string."""
        print(f"\n{'='*60}\n🚀 Coordinating pipeline for sequence:\n{fasta_string[:50]}...\n{'='*60}")
        self.agent_executor.invoke({"input": fasta_string})
        # After the agent's tool-based workflow is done, we perform the final ranking.
        self.perform_final_ranking()

    def perform_final_ranking(self):
        """
        Role 4 & 5: Ranker and Summarizer.
        This final step uses the LLM's reasoning ability to synthesize all the data
        collected by the tools into a final, human-readable assessment.
        """
        print("\n" + "-"*30)
        print("🧠 Performing Final Analysis and Ranking...")
        
        for header, data in PIPELINE_STATE.items():
            # Only rank sequences that successfully completed the entire pipeline.
            if data.get("status") == "docking_complete":
                blast_hits = data.get("blast_hits", [])
                top_hit = blast_hits[0] if blast_hits else {}
                
                ranking_prompt = f"""
                Analyze the drug discovery data for sequence '{header}':
                - Top Homolog: Species '{top_hit.get('sscinames', 'N/A')}', Identity {top_hit.get('pident', 'N/A')}%
                - Docking Score: {data.get('docking_result', {}).get('docking_score_kcal_mol', 'N/A')} kcal/mol

                Task: Provide a final assessment.
                1. Write a one-sentence summary of the candidate's potential.
                2. Assign a 'priority_score' from 0.0 to 10.0 (higher is better). A strong candidate
                   has high identity (>90%) to 'Homo sapiens' AND a strong docking score (<-9.0).

                Respond ONLY with a valid JSON object with keys "summary" and "priority_score".
                """
                
                response_str = llm.invoke(ranking_prompt).content
                try:
                    assessment = json.loads(response_str)
                    PIPELINE_STATE[header]['final_assessment'] = assessment
                    PIPELINE_STATE[header]['status'] = 'ranking_complete'
                    print(f"✅ Assessment for '{header}': {assessment['summary']} (Score: {assessment['priority_score']})")
                except json.JSONDecodeError:
                    print(f"❌ Failed to parse final ranking JSON for '{header}'.")


#
# 5. Simulation Runner
#
def main():
    """Initializes the coordinator and runs the simulation."""
    coordinator = PipelineCoordinator(llm_model=llm)

    candidate_1_fasta = """>candidate_001|KRAS_G12C_variant
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVIS
"""
    candidate_2_invalid_fasta = ">invalid_seq\nThis contains B and Z which are not valid characters."

    # Run the pipeline for both candidates
    coordinator.run_pipeline(candidate_1_fasta)
    coordinator.run_pipeline(candidate_2_invalid_fasta)

    # Print the final, complete state object showing all collected data
    print("\n\n" + "="*60)
    print("📋 Final State of All Processed Sequences:")
    print(json.dumps(PIPELINE_STATE, indent=2))
    print("="*60)


if __name__ == "__main__":
    main()
