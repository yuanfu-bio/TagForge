from __future__ import annotations

import logging
from pathlib import Path


def sample_logger(sample: str, path: Path):
    logger = logging.getLogger(f"tagforge.{sample}")
    logger.setLevel(logging.INFO); logger.propagate = False
    if not logger.handlers:
        path.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(path, encoding="utf-8")
        handler.setFormatter(logging.Formatter(f"%(asctime)s\t%(levelname)s\t{sample}\t%(message)s"))
        logger.addHandler(handler)
    return logger

