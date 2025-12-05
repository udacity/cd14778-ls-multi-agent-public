# main_starter.py

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
# These functions are assumed to be in a file named `ls_action_space.py`
from ls_action_space import query_pubmed, query_clinvar

# ----------------------------------------------------------------------------
# STEP 1: INITIALIZE THE LLM
# This is provided for you. The model is configured for the course environment.
# ----------------------------------------------------------------------------

llm = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=os.getenv("OPENAI_API_KEY"),
    base_url="https://openai.vocareum.com/v1",
    temperature=0.0
)

# ----------------------------------------------------------------------------
# STEP 2: DEFINE THE LOW-LEVEL RETRIEVAL TOOLS
# These are the basic tools that interact with the external APIs.
# ----------------------------------------------------------------------------

@tool
def pubmed_search_tool(query: str) -> List[Dict[str, Any]]:
    """
    Searches PubMed for scientific articles related to a given query (e.g., a disease).
    """
    return query_pubmed(query, max_results=5)

@tool
def clinvar_search_tool(variant_or_gene: str) -> Dict[str, Any]:
    """
    Searches ClinVar for information about a specific genetic variant or gene.
    """
    return query_clinvar(variant_or_gene)

# ----------------------------------------------------------------------------
# TODO: IMPLEMENT THE AGENTIC SYSTEM
# Your task is to build the multi-agent system below.
# ----------------------------------------------------------------------------

# TODO: Part 1 - Create the PubMed Retriever Agent
# This agent's only job is to search PubMed.
# 1. Define a simple ChatPromptTemplate for its system message.
# 2. Create the agent using create_tool_calling_agent().
# 3. Create the AgentExecutor.
#
# pubmed_retriever_prompt = ...
# pubmed_retriever_agent = ...
# pubmed_agent_executor = ...


# TODO: Part 2 - Create the ClinVar Retriever Agent
# This agent's only job is to search ClinVar.
# 1. Define a simple ChatPromptTemplate for its system message.
# 2. Create the agent using create_tool_calling_agent().
# 3. Create the AgentExecutor.
#
# clinvar_retriever_prompt = ...
# clinvar_retriever_agent = ...
# clinvar_agent_executor = ...


# TODO: Part 3 - Create Tools for the Synthesis Agent
# These tools will *delegate* work to the retriever agents you just created.
# This is the key to the multi-agent design.
@tool
def run_pubmed_retrieval(disease_query: str) -> str:
    """
    Delegates the task of searching PubMed to the specialized PubMed Retriever Agent.
    """
    print(f"--- Synthesis Agent delegating to PubMed Retriever for: {disease_query} ---")
    # Hint: Use pubmed_agent_executor.invoke(...) and return the 'output'.
    # response = ...
    # return json.dumps(response['output'])
    pass # Remove this line after implementing

@tool
def run_clinvar_retrieval(variant_query: str) -> str:
    """
    Delegates the task of searching ClinVar to the specialized ClinVar Retriever Agent.
    """
    print(f"--- Synthesis Agent delegating to ClinVar Retriever for: {variant_query} ---")
    # Hint: Use clinvar_agent_executor.invoke(...) and return the 'output'.
    # response = ...
    # return json.dumps(response['output'])
    pass # Remove this line after implementing


# TODO: Part 4 - Create and Run the Synthesis Agent (Orchestrator)
# This agent will use the delegation tools to gather data and write the final report.
def generate_rare_disease_dossier(disease: str, gene: str):
    """
    Runs the full multi-agent RAG pipeline to generate a dossier.
    """
    print(f"🚀 Starting dossier generation for {disease} ({gene})...")

    # 1. Define the ChatPromptTemplate for the synthesis agent.
    #    Instruct it to call the retrieval tools and then format a final report
    #    with sections for "Disease Overview", "Literature Summary", and "Genetic Variant Report".
    #
    # synthesis_prompt = ...

    # 2. Create the synthesis agent. Its tools are the two delegation tools you made.
    #
    # synthesis_agent_tools = [run_pubmed_retrieval, run_clinvar_retrieval]
    # synthesis_agent = ...

    # 3. Create the AgentExecutor for the synthesis agent.
    #
    # synthesis_agent_executor = ...

    # 4. Invoke the synthesis agent with a clear, user-facing query.
    query = f"Generate an evidence dossier for the rare disease '{disease}', focusing on the gene '{gene}'."
    # result = synthesis_agent_executor.invoke({"input": query})

    # 5. Print the final output from the result.
    print("\n✅ Dossier Generation Complete. Final Report:")
    # print(result['output'])


# ----------------------------------------------------------------------------
# Main execution block
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    # Example Case: Marfan Syndrome, caused by mutations in the FBN1 gene.
    # After you implement the TODOs, running this file will execute your system.
    generate_rare_disease_dossier(
        disease="Marfan Syndrome",
        gene="FBN1"
    )
