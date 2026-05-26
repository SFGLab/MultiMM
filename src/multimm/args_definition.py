import logging

from .config import SimulationConfig

logger = logging.getLogger(__name__)

# Expose a default instance of SimulationConfig for compatibility
args = SimulationConfig()
