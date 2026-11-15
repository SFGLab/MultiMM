from .bridge import SimulationEngine as SimulationEngine


def main(*args, **kwargs):
    """Command-line entry point (imported lazily so `python -m multimm.run` stays warning-free)."""
    from .run import main as _main

    return _main(*args, **kwargs)


__version__ = "2.1.1"
