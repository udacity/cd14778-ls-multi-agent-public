# logging_utils.py

import json
from datetime import datetime, timezone
from typing import Dict, Any

class JSONLLogger:
    """A simple, append-only JSONL logger for auditable agent actions."""
    def __init__(self, filepath: str, seed: int):
        self.filepath = filepath
        self.seed = seed

    def _digest_value(self, value: Any) -> Any:
        """Creates a short, readable summary of a value for logging."""
        if isinstance(value, str):
            return f"{value[:50]}..." if len(value) > 53 else value
        if isinstance(value, list):
            return f"[...{len(value)} items]" if len(value) > 5 else value
        if isinstance(value, dict):
             return {k: self._digest_value(v) for i, (k, v) in enumerate(value.items()) if i < 3}
        return value

    def _digest_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Creates a digest of a dictionary for logging."""
        if not isinstance(data, dict):
            return {"value": self._digest_value(data)}
        return {key: self._digest_value(value) for key, value in data.items()}

    def log_action(self, agent: str, tool: str, target: str,
                   inputs: Dict[str, Any], outputs: Dict[str, Any],
                   latency_ms: int):
        """Logs a single agent action to the JSONL file."""
        log_entry = {
            "ts": datetime.now(timezone.utc).isoformat(),
            "agent": agent,
            "tool": tool,
            "target": target,
            "inputs_digest": self._digest_data(inputs),
            "outputs_digest": self._digest_data(outputs),
            "latency_ms": latency_ms,
            "seed": self.seed
        }
        with open(self.filepath, 'a') as f:
            f.write(json.dumps(log_entry) + '\n')

    def log_orchestrator_start(self, run_id: str, input_counts: Dict[str, int]):
        """Logs the start of the main workflow."""
        log_entry = {
            "ts": datetime.now(timezone.utc).isoformat(), "agent": "Orchestrator",
            "tool": "start", "target": run_id, "inputs_digest": input_counts,
            "outputs_digest": {}, "latency_ms": 0, "seed": self.seed
        }
        with open(self.filepath, 'a') as f:
            f.write(json.dumps(log_entry) + '\n')

    def log_orchestrator_finish(self, run_id: str, status: str):
        """Logs the completion of the main workflow."""
        log_entry = {
            "ts": datetime.now(timezone.utc).isoformat(), "agent": "Orchestrator",
            "tool": "finish", "target": run_id, "inputs_digest": {},
            "outputs_digest": {"status": status}, "latency_ms": 0, "seed": self.seed
        }
        with open(self.filepath, 'a') as f:
            f.write(json.dumps(log_entry) + '\n')