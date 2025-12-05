# High-Throughput Screening Plate Allocator - Solution Code
# =========================================================
# This solution demonstrates state coordination for a multi-agent system
# managing High-Throughput Screening (HTS) plate resources in a lab.

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

# --- Configuration ---
# Load your OpenAI API key from an environment variable
# Ensure you have a .env file with OPENAI_API_KEY="your-key"
import dotenv
dotenv.load_dotenv("../../.env")
openai_api_key = os.getenv("OPENAI_API_KEY")

# If using the Udacity environment, the base URL is required
# If running locally with a standard OpenAI key, you can remove the 'base_url' parameter
model = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=openai_api_key,
    # base_url="https://openai.vocareum.com/v1", # Uncomment if in Udacity env
    temperature=0.0
)

# ==================================================
# SECTION 1: SHARED STATE MANAGEMENT
# ==================================================
# This section defines the shared resources that all agents will access and modify.
# The core challenge is to manage this state safely, preventing conflicts.

@dataclass
class ScreeningJob:
    """Represents a single screening job request."""
    job_id: str
    compound_ids: List[str]
    plate_id: str
    wells: List[str] = field(default_factory=list)
    status: str = "pending"  # pending, scheduled, rejected, completed
    priority: int = 1  # 1=Normal, 2=High

@dataclass
class HTSLabState:
    """
    Manages the state of the HTS lab, including plates and well bookings.
    This is the central, shared resource that requires careful coordination.
    """
    # Plate configurations available in the lab.
    plate_layouts: Dict[str, Dict[str, Any]] = field(default_factory=lambda: {
        "PLATE-96-A": {"rows": 8, "cols": 12},
        "PLATE-96-B": {"rows": 8, "cols": 12},
        "PLATE-384-A": {"rows": 16, "cols": 24},
    })
    
    # The critical shared resource: tracks which well is booked by which job.
    # Structure: {'plate_id': {'well_coordinate': 'job_id'}}
    # e.g., {'PLATE-96-A': {'A1': 'JOB-001', 'H12': 'JOB-002'}}
    well_bookings: Dict[str, Dict[str, str]] = field(default_factory=dict)
    
    # A queue of jobs to be processed.
    job_queue: List[ScreeningJob] = field(default_factory=list)
    
    # A simple counter to generate unique job IDs.
    job_counter: int = 0
    
    # A lock to prevent race conditions when multiple agents try to
    # access or modify well_bookings simultaneously. This is the key
    # mechanism for state coordination and concurrency control.
    lock: threading.Lock = field(default_factory=threading.Lock)

    def to_dict(self) -> Dict[str, Any]:
        """Serializes the state for easy printing."""
        return {
            "plates": list(self.plate_layouts.keys()),
            "well_bookings": self.well_bookings,
            "job_queue": [job.__dict__ for job in self.job_queue],
            "job_counter": self.job_counter,
        }

# --- Initialize the singleton shared state for the entire system ---
lab_state = HTSLabState()

# ==================================================
# SECTION 2: AGENT TOOLS (THE BUILDING BLOCKS)
# ==================================================
# These tools are functions that agents can call. They directly interact
# with the shared HTSLabState, using the lock to ensure thread safety.

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

def _well_generator(rows: int, cols: int):
    """Helper to generate well coordinates like A1, A2, ..., H12."""
    for row in range(rows):
        for col in range(1, cols + 1):
            yield f"{chr(ord('A') + row)}{col}"

@tool
def find_available_wells(plate_id: str, num_wells_needed: int) -> Dict[str, Any]:
    """
    Finds a contiguous block of available wells on a specified plate.
    
    Args:
        plate_id: The ID of the plate to search.
        num_wells_needed: The number of wells required.
        
    Returns:
        A dictionary containing the list of available wells or an error message.
    """
    if plate_id not in lab_state.plate_layouts:
        return {"success": False, "error": f"Plate '{plate_id}' not found."}

    layout = lab_state.plate_layouts[plate_id]
    all_wells = list(_well_generator(layout["rows"], layout["cols"]))
    
    with lab_state.lock:
        booked_wells = lab_state.well_bookings.get(plate_id, {}).keys()
        available_wells = [well for well in all_wells if well not in booked_wells]

    if len(available_wells) < num_wells_needed:
        return {
            "success": False,
            "error": "Insufficient available wells",
            "wells_available": len(available_wells)
        }
    
    return {"success": True, "wells": available_wells[:num_wells_needed]}

@tool
def allocate_wells(
    job_id: str,
    plate_id: str,
    compound_ids: List[str],
    priority: int
) -> Dict[str, Any]:
    """
    Allocates a set of wells on a plate to a screening job. This is the
    core conflict resolution function. It checks for availability and
    books the wells in a single, atomic operation.
    
    Args:
        job_id: The ID for the new job.
        plate_id: The plate where compounds will be screened.
        compound_ids: A list of compound identifiers to be screened.
        priority: The job priority (1=Normal, 2=High).
        
    Returns:
        A dictionary confirming the booking or detailing the conflict.
    """
    num_wells_needed = len(compound_ids)
    
    # Step 1: Find available wells
    well_search_result = find_available_wells.invoke(
        {"plate_id": plate_id, "num_wells_needed": num_wells_needed}
    )
    if not well_search_result["success"]:
        return {
            "success": False,
            "status": "rejected",
            "reason": "Conflict detected. Not enough available wells on the plate.",
            "details": well_search_result.get("error")
        }
    
    wells_to_book = well_search_result["wells"]

    # Step 2: Create the job object
    new_job = ScreeningJob(
        job_id=job_id,
        compound_ids=compound_ids,
        plate_id=plate_id,
        wells=wells_to_book,
        status="scheduled",
        priority=priority,
    )

    # Step 3: Book the wells within a lock to ensure atomicity
    with lab_state.lock:
        # Re-check for conflicts inside the lock to prevent race conditions
        current_bookings = lab_state.well_bookings.get(plate_id, {})
        for well in wells_to_book:
            if well in current_bookings:
                return {
                    "success": False,
                    "status": "rejected",
                    "reason": f"Conflict detected. Well {well} on {plate_id} was booked by another process."
                }

        # If no conflicts, update the state
        if plate_id not in lab_state.well_bookings:
            lab_state.well_bookings[plate_id] = {}
        for well in wells_to_book:
            lab_state.well_bookings[plate_id][well] = job_id
        
        lab_state.job_queue.append(new_job)
        # Sort queue by priority (desc) and then by job ID (asc) to ensure high-priority jobs go first
        lab_state.job_queue.sort(key=lambda j: (-j.priority, j.job_id))

    return {
        "success": True,
        "status": "scheduled",
        "job_id": job_id,
        "plate_id": plate_id,
        "wells_assigned": wells_to_book,
    }


# ==================================================
# SECTION 3: AGENT DEFINITIONS
# ==================================================
# We define specialized agents. In a real system, these might run as
# separate services. Here, they are logical roles handled by a single orchestrator.

# --- Agent 1: Request Parsing Agent ---
# Role: To understand the user's natural language request and extract key parameters.
request_parser_tools = [generate_job_id, list_available_plates]
request_parser_prompt = ChatPromptTemplate.from_messages([
    ("system", "You are an assistant that parses screening requests. Your goal is to identify the compounds, determine if it's high priority, and generate a job ID."),
    ("user", "{input}"),
    ("placeholder", "{agent_scratchpad}"),
])
request_parser_agent = create_tool_calling_agent(model, request_parser_tools, request_parser_prompt)
request_parser_executor = AgentExecutor(agent=request_parser_agent, tools=request_parser_tools, verbose=False)


# --- Agent 2: Plate Scheduler Agent ---
# Role: To find available space on plates and formally allocate the wells.
scheduler_tools = [find_available_wells, allocate_wells]
scheduler_prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a plate scheduling expert. Your job is to find available wells and allocate them to a job. You must use the `allocate_wells` tool to finalize any booking."),
    ("user", "{input}"),
    ("placeholder", "{agent_scratchpad}"),
])
scheduler_agent = create_tool_calling_agent(model, scheduler_tools, scheduler_prompt)
scheduler_executor = AgentExecutor(agent=scheduler_agent, tools=scheduler_tools, verbose=False)


# ==================================================
# SECTION 4: THE ORCHESTRATOR
# ==================================================
# The orchestrator coordinates the agents to fulfill a user request.
# It manages the flow of information and handles the final response.

class HTSOrchestrator:
    """Coordinates agents to process concurrent screening requests."""
    
    def extract_details(self, text: str) -> Dict[str, Any]:
        """A simple helper to extract structured data from agent responses."""
        details = {}
        # Extract Job ID
        job_id_match = re.search(r"JOB-\d{3}", text)
        if job_id_match:
            details["job_id"] = job_id_match.group(0)
        # Extract Compound IDs (assuming format C-XXX)
        details["compound_ids"] = re.findall(r'C-\d{3}', text)
        return details

    def process_screening_request(self, request_text: str) -> str:
        """
        Processes a user request from start to finish.
        This method demonstrates the core orchestration logic.
        """
        print(f"\n🚀 Orchestrator received request: \"{request_text}\"")

        # Step 1: Determine Priority (Mirroring real-world priority queuing)
        is_high_priority = "high-priority" in request_text.lower() or "must-test" in request_text.lower()
        priority = 2 if is_high_priority else 1
        print(f"   - Priority detected: {'High' if is_high_priority else 'Normal'}")

        # Step 2: Use the Request Parser Agent to understand the request
        parser_input = f"Parse the following request, determine the compound IDs, and generate a job ID for it. Request: '{request_text}'"
        parser_response = request_parser_executor.invoke({"input": parser_input})
        
        parsed_details = self.extract_details(parser_response['output'])
        job_id = parsed_details.get("job_id")
        compound_ids = parsed_details.get("compound_ids")

        if not job_id or not compound_ids:
            return f"Error: Could not parse the request to get a Job ID and Compound IDs. Parser output: {parser_response['output']}"
        
        print(f"   - Request Parser Agent identified {len(compound_ids)} compounds and assigned Job ID: {job_id}")

        # Step 3: Use the Scheduler Agent to find space and allocate wells
        # For simplicity, we'll try the first available 96-well plate.
        target_plate = "PLATE-96-A" 
        scheduler_input = (
            f"Please schedule job '{job_id}' for compounds {compound_ids} with priority {priority} "
            f"on plate '{target_plate}'. Use the 'allocate_wells' tool to finalize the booking."
        )
        
        print(f"   - Handing off to Scheduler Agent to book wells on {target_plate}...")
        scheduler_response = scheduler_executor.invoke({"input": scheduler_input})
        
        # Step 4: Interpret the result and formulate the final response
        # The tool output is often buried in the 'tool_output' part of the response
        final_output = scheduler_response.get('output', '')
        
        if "scheduled" in final_output.lower():
            wells_match = re.search(r"wells_assigned': \[(.*?)\]", final_output, re.DOTALL)
            wells_str = wells_match.group(1) if wells_match else "N/A"
            wells_str = wells_str.replace(' ', '').replace('\'', '')
            return f"✅ SUCCESS: Job {job_id} was successfully scheduled on plate {target_plate} in wells: {wells_str}."
        elif "conflict" in final_output.lower() or "rejected" in final_output.lower():
            return f"❌ CONFLICT: Job {job_id} could not be scheduled. The system reported a booking conflict or lack of space on plate {target_plate}."
        else:
            return f"⚠️ UNKNOWN: The scheduler returned an unexpected response for Job {job_id}: {final_output}"


# ==================================================
# SECTION 5: DEMONSTRATION
# ==================================================

def run_demo():
    """Simulates concurrent requests to demonstrate state coordination."""
    orchestrator = HTSOrchestrator()
    
    print("=" * 50)
    print("🔬 Welcome to the HTS Plate Allocator System 🔬")
    print("=" * 50)
    print("Initial Lab State:")
    print(lab_state.to_dict())
    print("-" * 50)

    # --- Simulate a series of concurrent requests ---
    requests = [
        "Please screen compounds C-001, C-002, C-003, C-004.",
        "This is a must-test batch, please run a high-priority screen for compounds C-101, C-102, C-103.",
        "We have a new set of 6 compounds to test: C-005 through C-010.",
        # This request will likely fail due to lack of space on PLATE-96-A
        "Requesting screen for 5 compounds: C-201, C-202, C-203, C-204, C-205."
    ]

    for req in requests:
        response = orchestrator.process_screening_request(req)
        print(f"   - Final Response: {response}")
    
    print("\n" + "=" * 50)
    print("✅ Demo Complete. Final Lab State:")
    print("-" * 50)
    import json
    print(json.dumps(lab_state.to_dict(), indent=2))
    print("\nNotice that the high-priority JOB-002 is at the top of the job_queue, even though it was the second request.")


if __name__ == "__main__":
    run_demo()
