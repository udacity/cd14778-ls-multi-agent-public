# High-Throughput Screening Plate Allocator - Starter Code
# =========================================================
# OBJECTIVE: Implement state coordination mechanisms to handle concurrency
# and resolve conflicts in a multi-agent system.

import os
import re
import threading
from dataclasses import dataclass, field
from typing import Dict, List, Any

# --- LangChain Imports ---
from langchain_core.prompts import ChatPromptTemplate
from langchain.agents import AgentExecutor, create_tool_calling_agent
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI

# ==================================================
# SECTION 1: CONFIGURATION (PROVIDED)
# ==================================================

import dotenv
dotenv.load_dotenv("../../.env")
openai_api_key = os.getenv("OPENAI_API_KEY")

model = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=openai_api_key,
    # base_url="https://openai.vocareum.com/v1", # Uncomment if in Udacity env
    temperature=0.0
)

# ==================================================
# SECTION 2: SHARED STATE MANAGEMENT (PARTIALLY PROVIDED)
# ==================================================

@dataclass
class ScreeningJob:
    """Represents a single screening job request."""
    job_id: str
    compound_ids: List[str]
    plate_id: str
    wells: List[str] = field(default_factory=list)
    status: str = "pending"
    priority: int = 1  # 1=Normal, 2=High

@dataclass
class HTSLabState:
    """Manages the state of the HTS lab, including plates and well bookings."""
    plate_layouts: Dict[str, Dict[str, Any]] = field(default_factory=lambda: {
        "PLATE-96-A": {"rows": 8, "cols": 12},
        "PLATE-384-A": {"rows": 16, "cols": 24},
    })
    well_bookings: Dict[str, Dict[str, str]] = field(default_factory=dict)
    job_queue: List[ScreeningJob] = field(default_factory=list)
    job_counter: int = 0
    
    # TODO 1: Add the core state coordination mechanism.
    # A lock is needed to prevent multiple agents from trying to book the
    # same well at the exact same time (a "race condition").
    # Add a `threading.Lock` instance here.
    # lock: threading.Lock = field(default_factory=threading.Lock)


# Initialize the singleton shared state for the entire system
lab_state = HTSLabState()

# ==================================================
# SECTION 3: AGENT TOOLS (STUDENT TO IMPLEMENT)
# ==================================================
# Helper function (provided)
def _well_generator(rows: int, cols: int):
    for row in range(rows):
        for col in range(1, cols + 1):
            yield f"{chr(ord('A') + row)}{col}"

@tool
def find_available_wells(plate_id: str, num_wells_needed: int) -> Dict[str, Any]:
    """
    Finds a contiguous block of available wells on a specified plate.
    """
    # TODO 2: Implement the logic to find available wells.
    # 1. Check if plate_id exists in `lab_state.plate_layouts`. If not, return an error.
    # 2. Get the list of all wells for the plate using `_well_generator`.
    # 3. Get the set of currently booked wells for that plate from `lab_state.well_bookings`.
    # 4. Find the wells that are in the "all_wells" list but not in the "booked_wells" set.
    # 5. If you have enough available wells, return the first `num_wells_needed` of them.
    #    Otherwise, return an error message indicating insufficient space.
    pass


@tool
def allocate_wells(
    job_id: str,
    plate_id: str,
    compound_ids: List[str],
    priority: int
) -> Dict[str, Any]:
    """
    Allocates a set of wells on a plate to a screening job, handling conflicts.
    This is the most important function for the learning objective.
    """
    # TODO 3: Implement the well allocation and conflict resolution logic.
    # 1. Use the `find_available_wells` tool to see if there's space. If not, return a "Conflict" message.
    #
    # 2. **CRITICAL STEP**: Use a `with lab_state.lock:` block to ensure the
    #    following steps are "atomic" (they happen all at once without interruption).
    #
    # 3. Inside the lock, you MUST re-check for conflicts. Another agent might have
    #    booked the wells between your step 1 and step 2.
    #    - Get the current bookings for `plate_id`.
    #    - Loop through the wells you plan to book. If any is now in the current bookings,
    #      return a "Conflict" message. This is the core of conflict detection.
    #
    # 4. If there are still no conflicts, update the shared state:
    #    - Add the new bookings to `lab_state.well_bookings`.
    #    - Create a `ScreeningJob` instance.
    #    - Append the new job to `lab_state.job_queue`.
    #
    # 5. **PRIORITY QUEUING**: After adding the job, sort `lab_state.job_queue`
    #    by priority (highest first) to ensure "must-test" compounds are handled first.
    #
    # 6. Return a success message with the job_id and the wells that were assigned.
    pass

# Provided simple tools
@tool
def generate_job_id() -> str:
    """Generates a new, unique job ID for a screening request."""
    with lab_state.lock:
        lab_state.job_counter += 1
        return f"JOB-{lab_state.job_counter:03d}"

@tool
def list_available_plates() -> List[str]:
    """Returns a list of all plate IDs available in the lab."""
    return list(lab_state.plate_layouts.keys())


# ==================================================
# SECTION 4: AGENT & ORCHESTRATOR SETUP (PARTIALLY PROVIDED)
# ==================================================

# --- Agent Definitions (Provided) ---
request_parser_tools = [generate_job_id, list_available_plates]
request_parser_prompt = ChatPromptTemplate.from_messages([
    ("system", "You parse screening requests to identify compounds, priority, and generate a job ID."),
    ("user", "{input}"), ("placeholder", "{agent_scratchpad}"),
])
request_parser_agent = create_tool_calling_agent(model, request_parser_tools, request_parser_prompt)
request_parser_executor = AgentExecutor(agent=request_parser_agent, tools=request_parser_tools, verbose=False)

scheduler_tools = [find_available_wells, allocate_wells]
scheduler_prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a plate scheduling expert. You find available wells and use the `allocate_wells` tool to book them."),
    ("user", "{input}"), ("placeholder", "{agent_scratchpad}"),
])
scheduler_agent = create_tool_calling_agent(model, scheduler_tools, scheduler_prompt)
scheduler_executor = AgentExecutor(agent=scheduler_agent, tools=scheduler_tools, verbose=False)

# --- Orchestrator (Student to Implement) ---
class HTSOrchestrator:
    """Coordinates agents to process concurrent screening requests."""
    def extract_details(self, text: str) -> Dict[str, Any]:
        details = {"compound_ids": re.findall(r'C-\d{3}', text)}
        job_id_match = re.search(r"JOB-\d{3}", text)
        if job_id_match: details["job_id"] = job_id_match.group(0)
        return details

    def process_screening_request(self, request_text: str) -> str:
        """Processes a user request from start to finish."""
        print(f"\n🚀 Orchestrator received request: \"{request_text}\"")
        priority = 2 if "must-test" in request_text.lower() else 1

        # TODO 4: Implement the orchestration flow.
        # 1. Invoke the `request_parser_executor` to get a Job ID and compound list.
        #    - The input prompt should be something like:
        #      f"Parse this request and get a job ID: '{request_text}'"
        #
        # 2. Extract the `job_id` and `compound_ids` from the parser's response.
        #
        # 3. Invoke the `scheduler_executor` to book the wells.
        #    - The input prompt should include the job_id, compounds, priority, and a target plate
        #      (e.g., "PLATE-96-A").
        #    - It should explicitly tell the agent to use the `allocate_wells` tool.
        #
        # 4. Based on the scheduler's response (check for "Success" or "Conflict"),
        #    return a final, user-friendly status message.
        return "Orchestrator logic not yet implemented."

# ==================================================
# SECTION 5: DEMONSTRATION (PROVIDED)
# ==================================================

def run_demo():
    """Simulates requests to demonstrate state coordination."""
    orchestrator = HTSOrchestrator()
    print("=" * 50 + "\n🔬 Welcome to the HTS Plate Allocator System 🔬\n" + "=" * 50)

    requests = [
        "Please screen compounds C-001, C-002, C-003, C-004.",
        "This is a must-test batch for compounds C-101, C-102.",
        "Requesting screen for 5 compounds: C-201, C-202, C-203, C-204, C-205."
    ]

    for req in requests:
        response = orchestrator.process_screening_request(req)
        print(f"   - Final Response: {response}")
    
    print("\n" + "=" * 50 + "\n✅ Demo Complete. Final Lab State:\n" + "-" * 50)
    import json
    # A small fix to make the lock serializable for printing
    lab_state_dict = lab_state.__dict__.copy()
    lab_state_dict.pop('lock', None)
    lab_state_dict['job_queue'] = [j.__dict__ for j in lab_state_dict['job_queue']]
    print(json.dumps(lab_state_dict, indent=2))
    print("\nCheck that the 'must-test' job is at the top of the job_queue.")

if __name__ == "__main__":
    run_demo()
