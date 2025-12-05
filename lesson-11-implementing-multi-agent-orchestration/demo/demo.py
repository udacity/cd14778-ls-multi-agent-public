# qPCR Machine Booker demo
#
# - Centralized shared state owned by a scheduler agent
# - Atomic "check-and-reserve" with a lock (prevents double booking)
# - Tool-like functions for a pluggable feel
# - Orchestrator pipeline: parse → validate → reserve → queue → ETA
# - Natural-language parsing (lightweight) + graceful error handling
# - Concurrent scenarios: conflict vs. no-conflict
#

from __future__ import annotations

import os
import re
import json
import threading
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Iterable, List, Optional, Set, Tuple


# -----------------------------
# Shared plate state & helpers
# -----------------------------

def all_wells_96() -> Set[str]:
    """Return all wells for a 96-well plate as labels A1..H12."""
    rows = "ABCDEFGH"
    cols = range(1, 13)
    return {f"{r}{c}" for r in rows for c in cols}


def wells_in_row(row: str, start: int, end: int) -> Set[str]:
    """Build a contiguous range within a row (e.g., B1..B6)."""
    if row not in "ABCDEFGH":
        raise ValueError("Row must be A..H")
    if not (1 <= start <= 12 and 1 <= end <= 12 and start <= end):
        raise ValueError("Columns must be 1..12 and start <= end")
    return {f"{row}{c}" for c in range(start, end + 1)}


def parse_well_token(tok: str) -> Optional[Tuple[str, int]]:
    """
    Parse a well token like 'B6' into ('B', 6). Returns None if invalid.
    """
    m = re.fullmatch(r"([A-Ha-h])\s*([1-9]|1[0-2])", tok.strip())
    if not m:
        return None
    row = m.group(1).upper()
    col = int(m.group(2))
    return (row, col)


def expand_range_token(tok: str) -> Optional[Set[str]]:
    """
    Parse tokens like 'B1-B6' or 'C3..C8' into a set of wells.
    """
    tok = tok.strip().upper().replace("..", "-")
    m = re.fullmatch(r"([A-H])\s*(\d{1,2})\s*-\s*([A-H])\s*(\d{1,2})", tok)
    if not m:
        return None
    r1, c1, r2, c2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
    if r1 != r2:
        return None  # only row-contiguous ranges for simplicity
    return wells_in_row(r1, min(c1, c2), max(c1, c2))


def parse_well_list_from_text(text: str) -> Set[str]:
    """
    Extract wells from free text. Supports:
    - Single wells: 'B6', 'c12'
    - Ranges: 'B1-B6', 'C3..C8'
    - Commas/and/whitespace separated lists
    """
    # Normalize list separators
    norm = re.sub(r"\band\b", ",", text, flags=re.I)
    tokens = re.split(r"[,\n;/]+", norm)

    wells: Set[str] = set()
    for raw in tokens:
        tok = raw.strip()
        if not tok:
            continue

        # Try range first
        rng = expand_range_token(tok)
        if rng:
            wells |= rng
            continue

        # Single well?
        single = parse_well_token(tok)
        if single:
            wells.add(f"{single[0]}{single[1]}")
            continue

        # Soft match patterns like 'row B 1-6'
        m = re.fullmatch(
            r"(row\s+)?([A-Ha-h])\s*(\d{1,2})\s*-\s*(\d{1,2})", tok
        )
        if m:
            row = m.group(2).upper()
            c1, c2 = int(m.group(3)), int(m.group(4))
            wells |= wells_in_row(row, min(c1, c2), max(c1, c2))

    return wells


# -----------------------------
# Data model
# -----------------------------

@dataclass(frozen=True)
class BookingRequest:
    experiment_id: str
    assay_name: str
    wells: Set[str]


@dataclass
class ProductionItem:
    booking: BookingRequest
    eta: datetime
    status: str = "queued"
    created_at: datetime = field(default_factory=datetime.utcnow)


# -----------------------------
# Tool-like functions
# -----------------------------

def tool_calculate_eta(requested_wells: int,
                       queue_depth: int) -> datetime:
    """
    Simple ETA heuristic:
    - Base runtime: 45 minutes per 96-well run (qPCR cycles)
    - Setup buffer: 10 minutes
    - Scale with fraction of plate used
    - Add small queue penalty per existing item
    """
    base_minutes = 45
    setup_minutes = 10
    fraction = max(0.15, requested_wells / 96.0)
    run_minutes = int(base_minutes * fraction) + setup_minutes
    queue_penalty = queue_depth * 5
    return datetime.utcnow() + timedelta(
        minutes=run_minutes + queue_penalty
    )


def tool_snapshot_state(scheduler: "QPCRSchedulerAgent",
                        queue: List[ProductionItem]) -> str:
    return json.dumps({
        "available_wells": scheduler.available_count(),
        "production_queue": [{
            "exp_id": it.booking.experiment_id,
            "assay": it.booking.assay_name,
            "wells": len(it.booking.wells),
            "eta_utc": it.eta.isoformat(timespec="minutes"),
            "status": it.status,
        } for it in queue]
    }, indent=2)

# -----------------------------
# Agents
# -----------------------------

class QPCRSchedulerAgent:
    """
    Centralized scheduler that owns the plate's availability.
    Invariant: 'check' and 'reserve' happen atomically
    inside one lock-guarded critical section.
    """

    def __init__(self) -> None:
        self._available: Set[str] = all_wells_96()
        self._lock = threading.Lock()

    def check_and_reserve(self, req: BookingRequest) -> bool:
        """
        Atomically:
        1) Check that all requested wells are available.
        2) If so, reserve them; else, deny.
        """
        with self._lock:
            if req.wells.issubset(self._available):
                self._available -= req.wells
                print(
                    f"[Scheduler] RESERVED {len(req.wells)} wells for "
                    f"{req.experiment_id} ({req.assay_name})"
                )
                return True
            conflicted = sorted(req.wells - self._available)
            print(
                f"[Scheduler] DENIED {req.experiment_id} "
                f"(conflicts: {conflicted[:6]}"
                f"{'...' if len(conflicted) > 6 else ''})"
            )
            return False

    def release(self, wells: Iterable[str]) -> None:
        """Release wells back to availability (e.g., after run completes)."""
        with self._lock:
            self._available |= set(wells)

    def available_count(self) -> int:
        with self._lock:
            return len(self._available)


class ExperimentAgent:
    """
    An experiment/assay agent that requests wells for its run.
    """

    def __init__(self, experiment_id: str, assay_name: str) -> None:
        self.experiment_id = experiment_id
        self.assay_name = assay_name

    def request_booking(self, scheduler: QPCRSchedulerAgent,
                        wells: Set[str]) -> bool:
        preview = sorted(list(wells))[:6]
        print(
            f"[{self.experiment_id}] Requesting {len(wells)} wells "
            f"(sample: {preview}{'...' if len(wells) > 6 else ''})"
        )
        req = BookingRequest(
            experiment_id=self.experiment_id,
            assay_name=self.assay_name,
            wells=wells,
        )
        ok = scheduler.check_and_reserve(req)
        print(
            f"[{self.experiment_id}] Booking "
            f"{'GRANTED' if ok else 'DENIED'}"
        )
        return ok


# -----------------------------
# Orchestrator
# -----------------------------

class LabOrchestrator:
    """
    Orchestrates the pipeline:
    parse → validate → reserve → queue → ETA
    """

    def __init__(self) -> None:
        self.scheduler = QPCRSchedulerAgent()
        self.production_queue: List[ProductionItem] = []

    # ---- Parsing & validation ----

    def parse_order(self, user_text: str) -> Tuple[str, str, Set[str]]:
        """
        Extract (experiment_id, assay_name, wells) from a free-text request.
        Expected hints in text, but robust to missing pieces.
        """
        # Experiment ID: e.g., EXP-123 or fallback
        m_id = re.search(r"\bEXP-[A-Za-z0-9_-]+\b", user_text)
        experiment_id = m_id.group(0) if m_id else "EXP-AUTO"

        # Assay name: take a simple phrase after 'assay' or 'for'
        m_assay = re.search(
            r"(assay|for)\s*[:\-]?\s*([A-Za-z0-9 _\-]+)",
            user_text, flags=re.I
        )
        assay_name = (m_assay.group(2).strip()
                      if m_assay else "qPCR run")

        # Wells: parse flexible list; require at least 1
        wells = parse_well_list_from_text(user_text)
        if not wells:
            # Provide a helpful message in the exception
            raise ValueError(
                "No wells found. Try formats like "
                "'B1-B6, C1-C6' or 'A1, A2, A3'."
            )
        return experiment_id, assay_name, wells

    def validate_wells(self, wells: Set[str]) -> None:
        """Ensure all wells are valid 96-well identifiers."""
        valid = all_wells_96()
        invalid = sorted(wells - valid)
        if invalid:
            raise ValueError(
                f"Invalid wells detected: {invalid[:6]}"
                f"{'...' if len(invalid) > 6 else ''}"
            )

    # ---- Reservation & queueing ----

    def reserve_and_enqueue(self, exp_id: str,
                            assay: str, wells: Set[str]) -> ProductionItem:
        """Reserve atomically; on success, enqueue with ETA."""
        req = BookingRequest(exp_id, assay, wells)
        if not self.scheduler.check_and_reserve(req):
            raise RuntimeError("Reservation conflict")

        eta = tool_calculate_eta(len(wells), len(self.production_queue))
        item = ProductionItem(booking=req, eta=eta)
        self.production_queue.append(item)
        return item

    # ---- Public API ----

    def handle_request_text(self, text: str) -> str:
        """
        End-to-end: parse, validate, reserve, enqueue, respond.
        Gracefully returns user-facing messages on failure.
        """
        print("\n--- Orchestrator: Processing Request ---")
        print(f"User: {text}")

        try:
            exp_id, assay, wells = self.parse_order(text)
            self.validate_wells(wells)
        except Exception as e:
            return (
                "Sorry, I couldn't understand your request.\n"
                f"Reason: {e}\n"
                "Example: 'Please book EXP-001 for assay GeneA: "
                "B1-B6, C1-B6'."
            )

        try:
            item = self.reserve_and_enqueue(exp_id, assay, wells)
        except RuntimeError:
            return (
                f"Conflict: Some requested wells are already reserved. "
                "Please pick a different set."
            )
        except Exception as e:
            return f"Unexpected error during reservation: {e}"

        eta_str = item.eta.isoformat(timespec="minutes") + "Z"
        preview = sorted(list(wells))[:8]
        response = (
            "✅ Booking confirmed!\n"
            f"- Experiment: {exp_id}\n"
            f"- Assay: {assay}\n"
            f"- Wells: {len(wells)} (sample: {preview}"
            f"{'...' if len(wells) > 8 else ''})\n"
            f"- ETA (UTC): {eta_str}\n"
        )

        print("\n--- State Snapshot (After Booking) ---")
        print(tool_snapshot_state(self.scheduler, self.production_queue))
        return response

    def state_snapshot(self) -> str:
        return tool_snapshot_state(self.scheduler, self.production_queue)


# -----------------------------
# Demo scenarios
# -----------------------------

def demo_conflict_threads() -> None:
    """
    Two agents concurrently request overlapping wells.
    Atomic reservation ensures only one succeeds.
    """
    print("\n=== Demo: Conflict (concurrent overlapping requests) ===")
    orch = LabOrchestrator()

    exp1 = ExperimentAgent("EXP-001", "GeneA quant")
    exp2 = ExperimentAgent("EXP-002", "GeneB quant")

    wells_exp1 = wells_in_row("B", 1, 6) | wells_in_row("C", 1, 6)
    wells_exp2 = wells_in_row("B", 1, 6) | wells_in_row("D", 1, 6)

    t1 = threading.Thread(
        target=lambda: exp1.request_booking(orch.scheduler, wells_exp1),
        daemon=True,
    )
    t2 = threading.Thread(
        target=lambda: exp2.request_booking(orch.scheduler, wells_exp2),
        daemon=True,
    )

    print("\n--- State (initial) ---")
    print(orch.state_snapshot())

    t1.start()
    t2.start()
    t1.join()
    t2.join()

    print("\n--- State (after threads) ---")
    print(orch.state_snapshot())


def demo_orchestrator_end_to_end() -> None:
    """
    Natural language → parsed wells → validation → reserve → queue → ETA.
    Shows: success, invalid wells, and conflict via sequential requests.
    """
    print("\n=== Demo: Orchestrator end-to-end ===")
    orch = LabOrchestrator()

    print("\n--- Initial State ---")
    print(orch.state_snapshot())

    # Success
    req1 = (
        "Please book EXP-101 for assay GeneX: B1-B6, C1-C6. "
        "Thanks!"
    )
    print("\n[REQ1] ", req1)
    print(orch.handle_request_text(req1))

    # Invalid wells
    req2 = "Book EXP-102 for GeneY: Z1-Z3, A1."
    print("\n[REQ2] ", req2)
    print(orch.handle_request_text(req2))

    # Conflict: overlaps with B1-B6 requested in REQ1
    req3 = "Reserve EXP-103 for assay GeneB: B1-B6, D1-D6."
    print("\n[REQ3] ", req3)
    print(orch.handle_request_text(req3))

    print("\n--- Final State ---")
    print(orch.state_snapshot())


def demo_no_conflict_threads() -> None:
    """
    Two agents request disjoint wells concurrently; both succeed.
    """
    print("\n=== Demo: No Conflict (disjoint concurrent requests) ===")
    orch = LabOrchestrator()

    exp1 = ExperimentAgent("EXP-201", "GeneX quant")
    exp2 = ExperimentAgent("EXP-202", "GeneY quant")

    wells_exp1 = wells_in_row("E", 1, 6) | wells_in_row("F", 1, 6)
    wells_exp2 = wells_in_row("G", 7, 12) | wells_in_row("H", 7, 12)

    print("\n--- State (initial) ---")
    print(orch.state_snapshot())

    t1 = threading.Thread(
        target=lambda: exp1.request_booking(orch.scheduler, wells_exp1),
        daemon=True,
    )
    t2 = threading.Thread(
        target=lambda: exp2.request_booking(orch.scheduler, wells_exp2),
        daemon=True,
    )

    t1.start()
    t2.start()
    t1.join()
    t2.join()

    # Enqueue via orchestrator for ETA to demonstrate full pipeline
    print("\n--- Orchestrator queuing for ETA ---")
    try:
        item1 = orch.reserve_and_enqueue("EXP-201-QUEUE",
                                         "GeneX queued", wells_exp1)
    except RuntimeError:
        pass  # already reserved by agent; this is fine for demo

    try:
        item2 = orch.reserve_and_enqueue("EXP-202-QUEUE",
                                         "GeneY queued", wells_exp2)
    except RuntimeError:
        pass  # already reserved by agent; this is fine for demo

    print("\n--- State (after) ---")
    print(orch.state_snapshot())


# -----------------------------
# Entry point
# -----------------------------

if __name__ == "__main__":
    # 1) show concurrent conflict (threads + atomic reservation)
    demo_conflict_threads()

    # 2) show orchestrator full flow incl. ETA + error handling
    demo_orchestrator_end_to_end()

    # 3) show concurrent no-conflict and queueing
    demo_no_conflict_threads()

