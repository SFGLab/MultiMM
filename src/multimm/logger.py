import logging
import sys

class _ColorFormatter(logging.Formatter):
    """Pretty colored formatter for terminal logs."""

    COLORS = {
        "DEBUG": "\033[1;34m",    # blue
        "INFO": "\033[1;32m",     # green
        "WARNING": "\033[1;33m",  # yellow
        "ERROR": "\033[1;31m",    # red
        "CRITICAL": "\033[1;41m", # red background
    }

    RESET = "\033[0m"
    BOLD = "\033[1m"

    def format(self, record):
        level_color = self.COLORS.get(record.levelname, "")
        
        record.levelname = f"{level_color}{record.levelname:<8}{self.RESET}"
        record.name = f"\033[1;35m{record.name}{self.RESET}"
        record.msg = str(record.msg)

        return super().format(record)


def setup_logger(level=logging.INFO, debug=False):
    """
    Clean, colored logger for simulation pipelines.
    """

    root = logging.getLogger()

    # Avoid duplicate handlers
    if root.handlers:
        return

    handler = logging.StreamHandler(sys.stdout)

    formatter = _ColorFormatter(
        "\033[1;36m[%(asctime)s]\033[0m "
        "%(levelname)s "
        "%(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    handler.setFormatter(formatter)

    root.setLevel(logging.DEBUG if debug else level)
    root.addHandler(handler)