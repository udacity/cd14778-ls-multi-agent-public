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
from langchain.tools import StructuredTool
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.pydantic_v1 import (BaseModel, Field)
from dotenv import load_dotenv

# --- Environment Setup ---
# Load environment variables from a .env file.
# Students should create a .env file with their Udacity OpenAI key:
# OPENAI_API_KEY="your_key_here"
load_dotenv("../../.env")
openai_api_key = os.getenv("OPENAI_API_KEY")

# --- Core State Management and Database Setup ---

DB_NAME = "mimic3_demo.db"
CSV_FILES = {
    'patients': '../PATIENTS.csv',
    'prescriptions': '../PRESCRIPTIONS.csv'
}
SIMULATION_PATIENT_ID = 10006

class PatientTimeline:
    """
    A thread-safe, in-memory data structure to manage patient medication timelines.
    This class acts as the shared state between different agents, ensuring that
    all operations (reads/writes) are synchronized to prevent race conditions.
    """
    def __init__(self):
        self._timeline = {}
        self._lock = threading.Lock()

    def add_prescription(self, patient_id, prescription):
        with self._lock:
            if patient_id not in self._timeline:
                self._timeline[patient_id] = []
            self._timeline[patient_id].append(prescription)

    def get_timeline(self, patient_id):
        with self._lock:
            return self._timeline.get(patient_id, []).copy()

def setup_database():
    """Initializes the SQLite database from MIMIC-III CSV files."""
    if os.path.exists(DB_NAME):
        print(f"[Setup] Database '{DB_NAME}' already exists. Skipping.")
        return
    print(f"[Setup] Database '{DB_NAME}' not found. Creating...")
    try:
        conn = sqlite3.connect(DB_NAME, check_same_thread=False)
        for table_name, file_name in CSV_FILES.items():
            if not os.path.exists(file_name):
                raise FileNotFoundError(f"Error: Required CSV '{file_name}' not found.")
            df = pd.read_csv(file_name)
            df.columns = df.columns.str.strip()
            df.to_sql(table_name, conn, if_exists='replace', index=False)
        print("[Setup] Database setup complete.")
        conn.close()
    except Exception as e:
        print(f"An error occurred during database setup: {e}")
        if os.path.exists(DB_NAME): os.remove(DB_NAME)
        exit()

# --- Agent Tools (Interacting with Shared State) ---
# Note: We now define input schemas for LangChain tools using Pydantic

class GetTimelineInput(BaseModel):
    # This tool requires no arguments from the LLM, as patient_id is fixed for the simulation
    pass

def get_medication_timeline(shared_timeline: PatientTimeline, patient_id: int) -> str:
    """Retrieves a formatted string of the patient's medication timeline."""
    timeline = shared_timeline.get_timeline(patient_id)
    if not timeline:
        return f"No medication data found for patient {patient_id}."
    
    response = f"--- Medication Timeline for Patient {patient_id} ---\n"
    sorted_timeline = sorted(timeline, key=lambda x: x['startdate'] or datetime.min)
    for p in sorted_timeline:
        start = p['startdate'].date() if p['startdate'] else 'N/A'
        end = p['enddate'].date() if p['enddate'] else 'N/A'
        response += f"- {p['drug']}: From {start} to {end}\n"
    return response

class CheckOverlapsInput(BaseModel):
    pass # No arguments needed from LLM

def check_for_drug_overlaps(shared_timeline: PatientTimeline, patient_id: int) -> str:
    """Checks for and reports any overlapping prescriptions in the patient's timeline."""
    timeline = shared_timeline.get_timeline(patient_id)
    if len(timeline) < 2:
        return "Not enough prescriptions to check for overlaps."

    overlaps = []
    sorted_timeline = sorted(timeline, key=lambda x: x['startdate'] or datetime.min)
    
    for i in range(len(sorted_timeline)):
        for j in range(i + 1, len(sorted_timeline)):
            p1, p2 = sorted_timeline[i], sorted_timeline[j]
            if not all([p1['startdate'], p1['enddate'], p2['startdate'], p2['enddate']]): continue
            if p1['startdate'] < p2['enddate'] and p1['enddate'] > p2['startdate']:
                overlaps.append(
                    f"  - Overlap: '{p1['drug']}' ({p1['startdate'].date()} to {p1['enddate'].date()}) "
                    f"and '{p2['drug']}' ({p2['startdate'].date()} to {p2['enddate'].date()})"
                )
    
    return "--- [ALERT] Potential Overlaps ---\n" + "\n".join(overlaps) if overlaps else "No drug overlaps detected."

class CheckGapsInput(BaseModel):
    min_gap_days: int = Field(description="The minimum number of days for a gap to be reported.")

def check_for_medication_gaps(shared_timeline: PatientTimeline, patient_id: int, min_gap_days: int) -> str:
    """Checks for medication gaps of a certain duration for the same drug."""
    timeline = shared_timeline.get_timeline(patient_id)
    if len(timeline) < 2: return "Not enough prescriptions to check for gaps."

    gaps, drugs = [], {}
    for p in timeline:
        drugs.setdefault(p['drug'], []).append(p)
    
    for drug_name, prescriptions in drugs.items():
        if len(prescriptions) < 2: continue
        
        sorted_prescriptions = sorted(prescriptions, key=lambda x: x['startdate'] or datetime.min)
        for i in range(len(sorted_prescriptions) - 1):
            p1, p2 = sorted_prescriptions[i], sorted_prescriptions[i+1]
            if p1['enddate'] and p2['startdate']:
                gap = p2['startdate'] - p1['enddate']
                if gap.days > min_gap_days:
                    gaps.append(
                        f"  - '{drug_name}': {gap.days} day gap "
                        f"between {p1['enddate'].date()} and {p2['startdate'].date()}."
                    )
    
    return f"--- [ALERT] Gaps > {min_gap_days} Days ---\n" + "\n".join(gaps) if gaps else f"No gaps > {min_gap_days} days detected."

# --- Procedural Agent (Populates Shared State) ---

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
        shared_timeline.add_prescription(patient_id, prescription)
        time.sleep(0.3)
    print(f"\n[Logger Agent] Finished processing historical data.")

# --- Main Application Logic ---

if __name__ == "__main__":
    if not openai_api_key:
        raise ValueError("OpenAI API key not found. Please set OPENAI_API_KEY in your .env file.")
    
    setup_database()
    
    patient_medication_timeline = PatientTimeline()
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
    
    # --- Initialize the LangChain Clinical Assistant Agent ---
    llm = ChatOpenAI(
        model="gpt-4o-mini",
        api_key=openai_api_key,
        base_url="https://openai.vocareum.com/v1"
    )
    
    # Bind the shared timeline and patient ID to the tool functions
    # This allows LangChain to call them without needing these context arguments
    def _get_timeline() -> str:
        return get_medication_timeline(
            shared_timeline=patient_medication_timeline,
            patient_id=SIMULATION_PATIENT_ID
        )

    def _check_overlaps() -> str:
        return check_for_drug_overlaps(
            shared_timeline=patient_medication_timeline,
            patient_id=SIMULATION_PATIENT_ID
        )

    def _check_gaps(min_gap_days: int) -> str:
        return check_for_medication_gaps(
            shared_timeline=patient_medication_timeline,
            patient_id=SIMULATION_PATIENT_ID,
            min_gap_days=min_gap_days
        )

    tools = [
        StructuredTool.from_function(
            name="get_medication_timeline",
            description="Get the full medication timeline for the current patient.",
            func=_get_timeline,
            args_schema=GetTimelineInput,  # empty schema → no args
        ),
        StructuredTool.from_function(
            name="check_for_drug_overlaps",
            description="Check for any overlapping medications in the patient's timeline.",
            func=_check_overlaps,
            args_schema=CheckOverlapsInput,  # empty schema → no args
        ),
        StructuredTool.from_function(
            name="check_for_medication_gaps",
            description="Check for gaps between prescriptions for the same drug.",
            func=_check_gaps,
            args_schema=CheckGapsInput,  # has min_gap_days: int
        ),
    ]

    prompt = ChatPromptTemplate.from_messages([
        ("system", "You are a helpful clinical assistant. You have access to a patient's real-time medication timeline. Use your tools to answer questions about it."),
        ("user", "{input}"),
        MessagesPlaceholder(variable_name="agent_scratchpad"),
    ])

    agent = create_openai_tools_agent(llm, tools, prompt)
    agent_executor = AgentExecutor(agent=agent, tools=tools, verbose=True)

    print("""
Example questions: e.g.
- “Show the full medication timeline for the current patient.” (tests get_medication_timeline)
- “Are there any overlapping prescriptions right now? If yes, list the pairs.” (tests check_for_drug_overlaps)
- “Check for medication gaps longer than 14 days.” (tests check_for_medication_gaps with a practical threshold)
""")
    print()

    print("\n[Clinical Assistant Agent] Ready. How can I help? (Type 'exit' to quit)")
    
    try:
        while True:
            query = input("> ")
            if query.lower() == 'exit': break
            
            response = agent_executor.invoke({"input": query})
            print(f"\n[Agent] {response['output']}")

    except KeyboardInterrupt:
        print("\n[Main] Shutting down simulation...")
    finally:
        stop_event.set()
        print("[Main] Simulation stopped.")


