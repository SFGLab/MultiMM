import logging
import os
import subprocess
import sys
from typing import Any, Dict

import openmm

from .config import SimulationConfig
from .model import MultiMM
from .run import write_config

logger = logging.getLogger(__name__)


class SimulationEngine:
    """Provides parameter validation, schema export, in-process and subprocess
    execution."""

    @classmethod
    def get_schema(cls) -> Dict[str, Any]:
        """Returns the JSON schema of the simulation configuration."""
        return SimulationConfig.model_json_schema()

    @classmethod
    def validate_params(cls, params: Dict[str, Any]) -> Dict[str, Any]:
        """Validates configuration parameters and returns a dictionary."""
        return SimulationConfig(**params).model_dump()

    @classmethod
    def run_in_process(cls, config_params: Dict[str, Any], fallback_to_cpu: bool = True) -> str:
        """Runs the simulation synchronously inside the current process.

        Handles type coercion, default values, optional fields,
        ensemble iteration.
        Logs are streamed to both standard logger and a file.

        Args:
            config_params: Dictionary of configuration parameters.
            fallback_to_cpu: If True, falls back to CPU if OpenCL/CUDA fails.

        Returns:
            str: The path to the auto-generated configuration file.
        """
        config = SimulationConfig(**config_params)

        os.makedirs(config.OUT_PATH, exist_ok=True)
        metadata_dir = os.path.join(config.OUT_PATH, "metadata")
        os.makedirs(metadata_dir, exist_ok=True)
        log_path = os.path.join(metadata_dir, "output.log")

        file_handler = logging.FileHandler(log_path, mode="w")
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        sim_logger = logging.getLogger("multimm")
        sim_logger.addHandler(file_handler)
        old_level = sim_logger.level
        if old_level == logging.NOTSET or old_level > logging.INFO:
            sim_logger.setLevel(logging.INFO)

        def attempt_run(cfg: SimulationConfig) -> bool:
            try:
                md = MultiMM(cfg)
                md.run()
                return True
            except openmm.OpenMMException as e:
                err_msg = str(e)
                current_platform = cfg.PLATFORM
                logger.error(f"Simulation failed on platform {current_platform} with error: {err_msg}")

                platform_errors = [
                    "Error initializing context",
                    "clGetPlatformInfo",
                    "Error launching CUDA compiler",
                    "CUDA error",
                ]
                is_platform_error = any(pe in err_msg for pe in platform_errors)

                if is_platform_error and fallback_to_cpu and current_platform in ("OpenCL", "CUDA"):
                    logger.warning(f"Platform {current_platform} failed. Falling back to CPU.")
                    cfg.PLATFORM = "CPU"
                    md_fallback = MultiMM(cfg)
                    md_fallback.run()
                    return True
                raise e
            except ValueError as e:
                if "Given point must have three values" in str(e):
                    logger.warning(f"Simulation finished but plotting failed: {e}")
                    return True
                raise e

        try:
            base_out_path = config.OUT_PATH

            write_config(config)

            if config.GENERATE_ENSEMBLE and config.N_ENSEMBLE is not None:
                start_seed = config.SHUFFLING_SEED
                for i in range(config.N_ENSEMBLE):
                    config.SHUFFLING_SEED = start_seed + i
                    config.OUT_PATH = f"{base_out_path}_{i+1}"

                    for attempt in range(3):
                        try:
                            if attempt_run(config):
                                break
                        except Exception as ex:
                            if attempt == 2:
                                raise ex
                            logger.warning(f"Ensemble {i+1} attempt {attempt + 1} failed, retrying... Error: {ex}")
            else:
                for attempt in range(3):
                    try:
                        if attempt_run(config):
                            break
                    except Exception as ex:
                        if attempt == 2:
                            raise ex
                        logger.warning(f"Attempt {attempt + 1} failed, retrying... Error: {ex}")

        finally:
            sim_logger.removeHandler(file_handler)
            file_handler.close()
            sim_logger.setLevel(old_level)

        return os.path.join(metadata_dir, "config_auto.ini")

    @classmethod
    def run_subprocess(cls, config_params: Dict[str, Any]) -> str:
        """Runs the simulation in a separate subprocess.

        Args:
            config_params: Dictionary of configuration parameters.

        Returns:
            str: The path to the auto-generated configuration file.
        """
        config = SimulationConfig(**config_params)

        os.makedirs(config.OUT_PATH, exist_ok=True)
        metadata_dir = os.path.join(config.OUT_PATH, "metadata")
        os.makedirs(metadata_dir, exist_ok=True)
        config_path = os.path.join(metadata_dir, "config_auto.ini")

        write_config(config)

        cmd = [
            sys.executable,
            "-m",
            "multimm.run",
            "-c",
            config_path,
        ]

        log_path = os.path.join(metadata_dir, "output.log")
        with open(log_path, "w") as log_file:
            subprocess.run(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                text=True,
                check=True,
            )

        return config_path
