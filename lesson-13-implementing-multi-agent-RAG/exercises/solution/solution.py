# main.py

import os
import json
from typing import Dict, Any, List

# --- Environment and API Setup ---
# Note: Ensure you have a .env file with your OPENAI_API_KEY
from dotenv import load_dotenv
load_dotenv("../../.env")

# --- LangChain Core Imports ---
from langchain_openai import ChatOpenAI
from langchain.agents import AgentExecutor, create_tool_calling_agent
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.tools import tool

# --- Import Provided Life Science Tools ---
# Wrappers for public APIs (PubMed, ClinVar).
from ls_action_space.action_space import query_pubmed, query_clinvar

# ----------------------------------------------------------------------------
# STEP 1: INITIALIZE THE LLM
# The model is configured to use the Udacity-provided OpenAI endpoint.
# ----------------------------------------------------------------------------

llm = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=os.getenv("OPENAI_API_KEY"),
    base_url="https://openai.vocareum.com/v1",
    temperature=0.0 # Set to 0 for deterministic and factual outputs
)

# ----------------------------------------------------------------------------
# STEP 2: CREATE SPECIALIZED RETRIEVER AGENTS
# Each retriever is an agent with a single tool, making it a specialist.
# This modular design is central to the multi-agent RAG pattern.
# ----------------------------------------------------------------------------

## PubMed Retriever Agent ##

# Wrap the imported function in a LangChain tool decorator
@tool
def pubmed_search_tool(query: str) -> List[Dict[str, Any]]:
    """
    Searches PubMed for scientific articles related to a given query,
    which is typically a disease name or medical concept.
    Returns a list of articles with their title, abstract, authors, and year.
    """
    return query_pubmed(query, max_results=5)

# System prompt to define the agent's persona and task
pubmed_retriever_prompt = ChatPromptTemplate.from_messages([
    (
        "system",
        """You are a specialized agent expert in biomedical literature research.
        Your sole purpose is to use the `pubmed_search_tool` to find relevant
        scientific articles for a given query.
        Summarize the key findings from the retrieved abstracts concisely."""
    ),
    ("human", "{input}"),
    ("placeholder", "{agent_scratchpad}"),
])

# Create the agent and its executor
pubmed_retriever_agent = create_tool_calling_agent(llm, [pubmed_search_tool], pubmed_retriever_prompt)
pubmed_agent_executor = AgentExecutor(agent=pubmed_retriever_agent, tools=[pubmed_search_tool], verbose=True)


## ClinVar Retriever Agent ##

@tool
def clinvar_search_tool(variant_or_gene: str) -> Dict[str, Any]:
    """
    Searches the ClinVar database for information about a specific genetic
    variant or gene. Accepts gene symbols (e.g., 'FBN1'), rsIDs (e.g., 'rs121913459'),
    or HGVS strings.
    Returns details on clinical significance, conditions, and review status.
    """
    return query_clinvar(variant_or_gene)

# System prompt for the ClinVar agent
clinvar_retriever_prompt = ChatPromptTemplate.from_messages([
    (
        "system",
        """You are a specialized agent expert in clinical genetics.
        Your sole purpose is to use the `clinvar_search_tool` to find data
        on a specific gene or genetic variant.
        Present the retrieved information in a clear, summarized format."""
    ),
    ("human", "{input}"),
    ("placeholder", "{agent_scratchpad}"),
])

# Create the agent and its executor
clinvar_retriever_agent = create_tool_calling_agent(llm, [clinvar_search_tool], clinvar_retriever_prompt)
clinvar_agent_executor = AgentExecutor(agent=clinvar_retriever_agent, tools=[clinvar_search_tool], verbose=True)


# ----------------------------------------------------------------------------
# STEP 3: CREATE THE SYNTHESIS AGENT (ORCHESTRATOR)
# This agent does not access the low-level tools directly. Instead, it uses
# tools that delegate tasks to the specialized retriever agents created above.
# ----------------------------------------------------------------------------

@tool
def run_pubmed_retrieval(disease_query: str) -> str:
    """
    Delegates the task of searching PubMed for a specific disease to the
    specialized PubMed Retriever Agent. Use this to gather literature context.
    """
    print(f"--- Synthesis Agent delegating to PubMed Retriever for: {disease_query} ---")
    response = pubmed_agent_executor.invoke({"input": disease_query})
    return json.dumps(response['output'], indent=2)

@tool
def run_clinvar_retrieval(variant_query: str) -> str:
    """
    Delegates the task of searching ClinVar for a gene or variant to the
    specialized ClinVar Retriever Agent. Use this to get genetic details.
    """
    print(f"--- Synthesis Agent delegating to ClinVar Retriever for: {variant_query} ---")
    response = clinvar_agent_executor.invoke({"input": variant_query})
    return json.dumps(response['output'], indent=2)

# System prompt for the Synthesis Agent
synthesis_prompt = ChatPromptTemplate.from_messages([
    (
        "system",
        """You are a clinical assistant AI responsible for creating rare disease evidence dossiers.
        Your task is to orchestrate specialized agents to gather information and then synthesize it
        into a concise report for a clinician.

        Follow these steps:
        1.  Receive a user request for a rare disease and an associated gene/variant.
        2.  Use the `run_pubmed_retrieval` tool to find literature about the disease.
        3.  Use the `run_clinvar_retrieval` tool to find information on the gene/variant.
        4.  Once all information is gathered, compile a final report titled "Rare Disease Evidence Dossier".
        5.  The report must be well-structured with the following sections:
            - **Disease Overview**: A brief summary based on the literature.
            - **Literature Summary**: Key findings from the top 2-3 PubMed articles.
            - **Genetic Variant Report**: A summary of the ClinVar data, including clinical significance and associated conditions.
        6.  Your final output should ONLY be the formatted report. Do not include conversational text before or after the dossier.
        """
    ),
    ("human", "{input}"),
    ("placeholder", "{agent_scratchpad}"),
])

# Create the synthesis agent and its executor
synthesis_agent_tools = [run_pubmed_retrieval, run_clinvar_retrieval]
synthesis_agent = create_tool_calling_agent(llm, synthesis_agent_tools, synthesis_prompt)
synthesis_agent_executor = AgentExecutor(
    agent=synthesis_agent,
    tools=synthesis_agent_tools,
    verbose=True
)

# ----------------------------------------------------------------------------
# STEP 4: DEFINE THE DEMONSTRATION FUNCTION
# This function showcases the multi-agent system in action.
# ----------------------------------------------------------------------------

def generate_rare_disease_dossier(disease: str, gene: str):
    """
    Runs the full multi-agent RAG pipeline to generate a dossier.
    
    Args:
        disease (str): The name of the rare disease to investigate.
        gene (str): The gene symbol or variant ID to investigate.
    """
    print("="*60)
    print(f"🚀 Starting dossier generation for {disease} ({gene})...")
    print("="*60)
    
    query = f"Generate an evidence dossier for the rare disease '{disease}', focusing on the gene '{gene}'."
    
    result = synthesis_agent_executor.invoke({"input": query})
    
    print("\n" + "="*60)
    print("✅ Dossier Generation Complete. Final Report:")
    print("="*60)
    print(result['output'])
    print("="*60)

# ----------------------------------------------------------------------------
# STEP 5: RUN THE DEMONSTRATION
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    # Example Case: Marfan Syndrome, caused by mutations in the FBN1 gene.
    generate_rare_disease_dossier(
        disease="Marfan Syndrome",
        gene="FBN1"
    )
