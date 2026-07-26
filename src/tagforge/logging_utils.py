from __future__ import annotations

import logging
import time
from logging.handlers import RotatingFileHandler
from pathlib import Path


def sample_logger(sample: str, path: Path):
    logger = logging.getLogger(f"tagforge.{sample}")
    logger.setLevel(logging.INFO); logger.propagate = False
    if not logger.handlers:
        path.parent.mkdir(parents=True, exist_ok=True)
        handler = RotatingFileHandler(
            path, maxBytes=5 * 1024 * 1024, backupCount=1, encoding="utf-8"
        )
        handler.setFormatter(logging.Formatter(f"%(asctime)s\t%(levelname)s\t{sample}\t%(message)s"))
        logger.addHandler(handler)
    return logger



class ProgressThrottle:
    """Rate-limit transient progress records while always preserving milestones."""
    def __init__(self, interval_seconds: float = 5.0):
        self.interval_seconds = interval_seconds
        self.last_emitted = float("-inf")

    def allow(self, status: str) -> bool:
        if status in {"completed", "failed", "aggregated", "aggregate_complete", "sort_complete"}:
            self.last_emitted = time.monotonic()
            return True
        now = time.monotonic()
        if now - self.last_emitted >= self.interval_seconds:
            self.last_emitted = now
            return True
        return False
