import sqlite3
import pandas as pd
import threading
import time
from datetime import datetime
import os
from typing import List, Type
import atexit
import functools

# --- LangChain and OpenAI Imports ---
# Make sure to `pip install langchain langchain-openai python-dotenv`
from langchain_openai import ChatOpenAI
from langchain.agents import AgentExecutor, create_openai_tools_agent
from langchain.tools import Tool
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.pydantic_v1 import BaseModel, Field
from dotenv import load_dotenv

# --- Environment Setup ---
# Students should create a .env file with their Udacity OpenAI key:
# OPENAI_API_KEY="your_key_here"
load_dotenv("../../.env")
openai_api_key = os.getenv("OPENAI_API_KEY")

# --- Core State Management and Database Setup (Provided) ---

DB_NAME = "mimic3_demo.db"
SIMULATION_PATIENT_ID = 10006

class PatientTimeline:
    """
    A thread-safe, in-memory data structure to manage patient medication timelines.
    This class acts as the shared state between the logging and monitoring agents.
    """
    def __init__(self):
        self._timeline = {}
        # A lock is crucial for preventing race conditions when multiple threads access the timeline
        self._lock = threading.Lock()

    # TODO 1: Implement the `add_prescription` method.
    # It should acquire the lock before modifying the `_timeline` dictionary
    # to ensure that writes are thread-safe.
    def add_prescription(self, patient_id, prescription):
        # Your implementation here
        pass

    # TODO 2: Implement the `get_timeline` method.
    # It should also acquire the lock to safely read from the `_timeline` dictionary.
    # Remember to return a copy of the list to prevent external modification.
    def get_timeline(self, patient_id):
        # Your implementation here
        return []

def setup_database():
    """Initializes the SQLite database. (Provided)"""
    if os.path.exists(DB_NAME): return
    print(f"[Setup] Database '{DB_NAME}' not found. Creating...")
    conn = sqlite3.connect(DB_NAME, check_same_thread=False)
    df = pd.read_csv('../PRESCRIPTIONS.csv')
    df.columns = df.columns.str.strip()
    df.to_sql('prescriptions', conn, if_exists='replace', index=False)
    conn.close()
    print("[Setup] Database setup complete.")

# --- Agent Tools (To Be Implemented by Student) ---

# TODO 3: Implement the logic for the agent's tools.
# Each of these functions will be called by the LangChain agent.
# They should use the `shared_timeline` object to safely get the patient data and analyze it.

class GetTimelineInput(BaseModel): pass

def get_medication_timeline(shared_timeline: PatientTimeline, patient_id: int) -> str:
    """Retrieves a formatted string of the patient's medication timeline."""
    # Hint: Call the get_timeline method you implemented.
    # Format the data into a readable string.
    return "Tool 'get_medication_timeline' not implemented yet."

class CheckOverlapsInput(BaseModel): pass

def check_for_drug_overlaps(shared_timeline: PatientTimeline, patient_id: int) -> str:
    """Checks for and reports any overlapping prescriptions."""
    # Hint: Get the timeline, then loop through it to find prescriptions with overlapping dates.
    return "Tool 'check_for_drug_overlaps' not implemented yet."

class CheckGapsInput(BaseModel):
    min_gap_days: int = Field(description="The minimum number of days for a gap to be reported.")

def check_for_medication_gaps(shared_timeline: PatientTimeline, patient_id: int, min_gap_days: int) -> str:
    """Checks for medication gaps of a certain duration for the same drug."""
    # Hint: This is more advanced. Group prescriptions by drug, sort by date, and check the time between them.
    return "Tool 'check_for_medication_gaps' not implemented yet."

# --- Procedural Agent (Provided) ---

def logging_agent(patient_id, shared_timeline, db_connection, stop_event):
    """Agent that runs in the background, logging prescriptions to the shared timeline."""
    print(f"[Logger Agent] Started. Simulating real-time data feed...")
    cursor = db_connection.cursor()
    query = "SELECT startdate, enddate, drug FROM prescriptions WHERE subject_id = ? ORDER BY startdate"
    cursor.execute(query, (patient_id,))
    
    for row in cursor.fetchall():
        if stop_event.is_set(): break
        prescription = {
            'startdate': datetime.fromisoformat(row[0]) if row[0] else None,
            'enddate': datetime.fromisoformat(row[1]) if row[1] else None,
            'drug': row[2]
        }
        # This calls the `add_prescription` method you will implement
        shared_timeline.add_prescription(patient_id, prescription)
        time.sleep(0.3)
    print(f"\n[Logger Agent] Finished processing historical data.")

# --- Main Application Logic (Provided) ---

if __name__ == "__main__":
    if not openai_api_key:
        raise ValueError("OpenAI API key not found. Please set OPENAI_API_KEY in your .env file.")
    
    setup_database()
    
    patient_medication_timeline = PatientTimeline() # Your implementation will be used here
    db_conn = sqlite3.connect(DB_NAME, check_same_thread=False)
    stop_event = threading.Event()
    atexit.register(lambda: (db_conn.close(), print("\n[Main] Database connection closed.")))

    print(f"--- Starting EHR Agentic AI Simulation for Patient ID: {SIMULATION_PATIENT_ID} ---")
    
    logging_thread = threading.Thread(
        target=logging_agent,
        args=(SIMULATION_PATIENT_ID, patient_medication_timeline, db_conn, stop_event)
    )
    logging_thread.daemon = True
    logging_thread.start()
    
    # Initialize the LangChain Agent
    llm = ChatOpenAI(model="gpt-4o-mini", api_key=openai_api_key, base_url="https://openai.vocareum.com/v1")
    
    # The tools list will use the functions you implement above
    tools = [
        Tool(name="get_medication_timeline", func=functools.partial(get_medication_timeline, shared_timeline=patient_medication_timeline, patient_id=SIMULATION_PATIENT_ID), description="Get the full medication timeline for the current patient.", args_schema=GetTimelineInput),
        Tool(name="check_for_drug_overlaps", func=functools.partial(check_for_drug_overlaps, shared_timeline=patient_medication_timeline, patient_id=SIMULATION_PATIENT_ID), description="Check for any overlapping medications in the patient's timeline.", args_schema=CheckOverlapsInput),
        Tool(name="check_for_medication_gaps", func=functools.partial(check_for_medication_gaps, shared_timeline=patient_medication_timeline, patient_id=SIMULATION_PATIENT_ID), description="Check for gaps between prescriptions for the same drug.", args_schema=CheckGapsInput),
    ]

    prompt = ChatPromptTemplate.from_messages([
        ("system", "You are a helpful clinical assistant. Use your tools to answer questions about the patient's real-time medication timeline."),
        ("user", "{input}"),
        MessagesPlaceholder(variable_name="agent_scratchpad"),
    ])

    agent = create_openai_tools_agent(llm, tools, prompt)
    agent_executor = AgentExecutor(agent=agent, tools=tools, verbose=True)

    print("\n[Clinical Assistant Agent] Ready. How can I help? (Type 'exit' to quit)")
    
    try:
        while True:
            query = input("> ")
            if query.lower() == 'exit': break
            agent_executor.invoke({"input": query})
    except KeyboardInterrupt:
        print("\n[Main] Shutting down simulation...")
    finally:
        stop_event.set()
        print("[Main] Simulation stopped.")

