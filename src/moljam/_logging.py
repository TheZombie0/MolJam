import logging
import sys


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(f"moljam.{name}")


_root = logging.getLogger("moljam")
if not _root.handlers:
    _handler = logging.StreamHandler(sys.stdout)
    _handler.setFormatter(logging.Formatter("%(message)s"))
    _root.addHandler(_handler)
    _root.setLevel(logging.INFO)
