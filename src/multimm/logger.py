import logging
import sys
import os
import time
from typing import Optional


# ── ANSI primitives ───────────────────────────────────────────────────────────

RESET  = "\033[0m"
BOLD   = "\033[1m"
DIM    = "\033[2m"
ITALIC = "\033[3m"

def fg(n: int) -> str: return f"\033[38;5;{n}m"
def bg(n: int) -> str: return f"\033[48;5;{n}m"


# ── Per-level styles  (colour-string, icon, 3-char tag) ───────────────────────

LEVEL_FMT = {
    "DEBUG":    (fg(63)  + DIM,              "·", "DBG"),
    "INFO":     (fg(114),                    "✦", "INF"),
    "WARNING":  (fg(214) + BOLD,             "▲", "WRN"),   # 256-color orange
    "ERROR":    (fg(196) + BOLD,             "✖", "ERR"),   # 256-color bright red
    "CRITICAL": (fg(231) + bg(88) + BOLD,    "☠", "CRT"),
}

# Fixed visual widths (all badges are 5 visible chars: "✦ INF")
_BADGE_W  = 5   # icon(1) + space(1) + tag(3)
_MODULE_W = 22  # truncated module name column


# ── Colour detection ──────────────────────────────────────────────────────────

def _tty_colors() -> bool:
    """True when the terminal is colour-capable.

    Respects the standard NO_COLOR / FORCE_COLOR / COLORTERM conventions,
    and falls back to checking stderr when stdout has been redirected.
    """
    if os.environ.get("NO_COLOR"):          # https://no-color.org
        return False
    if os.environ.get("FORCE_COLOR") or os.environ.get("COLORTERM"):
        return True                          # kitty, alacritty, truecolor terms
    if os.environ.get("TERM") in ("dumb", ""):
        return False
    # Check the *original* stderr — it's rarely redirected and reliable.
    # We use _STDERR_IS_TTY which is cached at import time (before any Tee
    # wrapping replaces sys.stderr/sys.stdout).
    return _STDERR_IS_TTY or _STDOUT_IS_TTY


# Cache TTY state at import time, BEFORE any Tee-wrapping of sys.stdout.
_STDOUT_IS_TTY: bool = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()
_STDERR_IS_TTY: bool = hasattr(sys.stderr, "isatty") and sys.stderr.isatty()

# Single cached flag used everywhere in this module.
_COLOR_ENABLED: bool = _tty_colors()


# ── Formatter ─────────────────────────────────────────────────────────────────

class CozyFormatter(logging.Formatter):
    """
    One structured line per record:

        [HH:MM:SS]  ✦ INF  module.sub            ›  message text

    Colour is detected once at formatter-creation time (or forced via
    force_color=True).  CRITICAL also emits a bordered banner.
    """

    def __init__(self, force_color: bool = False) -> None:
        super().__init__()
        self._color: bool = force_color or _COLOR_ENABLED

    # Colourize plain text; padding must happen BEFORE this call so that
    # ANSI codes do not distort f-string alignment specifiers.
    def _c(self, plain: str, *ansi: str) -> str:
        if not self._color:
            return plain
        return "".join(ansi) + plain + RESET

    def format(self, record: logging.LogRecord) -> str:  # noqa: A003
        lvl_color, icon, tag = LEVEL_FMT.get(record.levelname, ("", "•", "???"))

        # ── timestamp ─────────────────────────────────────────────────────────
        ts = time.strftime("%H:%M:%S", time.localtime(record.created))
        ts_ = self._c(f"[{ts}]", fg(242))

        # ── level badge: pad plain text FIRST, then colourize ─────────────────
        badge_plain = f"{icon} {tag}"          # e.g. "✦ INF"  (5 visible chars)
        badge = self._c(badge_plain, lvl_color)

        # ── module name: last two dotted segments, fixed-width ─────────────────
        parts  = record.name.split(".")
        mod    = ".".join(parts[-2:]) if len(parts) > 2 else record.name
        mod_   = self._c(f"{mod:{_MODULE_W}}", fg(111))   # pad plain, then colour

        # ── separator ─────────────────────────────────────────────────────────
        sep = self._c("›", fg(238))

        # ── message ───────────────────────────────────────────────────────────
        msg = record.getMessage()
        if record.exc_info:
            msg += "\n" + self.formatException(record.exc_info)

        if record.levelname == "DEBUG":
            msg_ = self._c(msg, fg(240) + ITALIC)
        elif record.levelname == "WARNING":
            msg_ = self._c(msg, fg(214))          # orange message text
        elif record.levelname in ("ERROR", "CRITICAL"):
            msg_ = self._c(msg, lvl_color)        # red / critical message text
        else:
            msg_ = msg                             # INFO: pass ANSI in msg as-is

        line = f"{ts_}  {badge}  {mod_} {sep}  {msg_}"

        # ── CRITICAL: bordered banner ──────────────────────────────────────────
        if record.levelname == "CRITICAL":
            w   = max(len(msg) + 8, 64)
            bar = self._c("━" * w, fg(196) + BOLD)
            lbl = self._c(f"  ☠  CRITICAL  ☠  {msg}  ", fg(231) + bg(88) + BOLD)
            line = f"\n{bar}\n{lbl}\n{bar}\n"

        return line


# ── Public API ────────────────────────────────────────────────────────────────

# Third-party libraries that log routine/advisory INFO-DEBUG chatter which
# isn't ours — e.g. matplotlib.category's "Using categorical units..." notice
# whenever a plot call gets numeric-looking string data. Since setup_logger
# attaches its formatter to the ROOT logger, every library that propagates up
# to root would otherwise get printed as if it were a MultiMM log line. Capped
# at WARNING so a real problem from these libraries still surfaces.
_NOISY_THIRD_PARTY_LOGGERS = (
    "matplotlib",
    "PIL",
    "numba",
    "h5py",
    "fontTools",
)


def _quiet_third_party_loggers() -> None:
    for name in _NOISY_THIRD_PARTY_LOGGERS:
        logging.getLogger(name).setLevel(logging.WARNING)


def setup_logger(
    level: int = logging.INFO,
    debug: bool = False,
    force_color: bool = False,
    log_file: Optional[str] = None,
) -> None:
    """
    Attach CozyFormatter to the root logger, replacing any existing handlers.

    Safe to call multiple times.  If *log_file* is given and the logger
    already has a CozyFormatter console handler, only the file handler is
    added (if not already present), so subsequent calls with different
    *log_file* paths will add a new file sink without re-initialising the
    console handler.

    Parameters
    ----------
    level       Minimum log level when *debug* is False (default INFO).
    debug       If True, sets level to DEBUG regardless of *level*.
    force_color Emit ANSI codes even when stdout is not a tty (e.g. piped).
    log_file    Optional path to a plain-text log file written alongside the
                coloured console output.  The file uses a simple timestamped
                format (no ANSI codes).  The parent directory is created
                automatically.  Recommended: ``<OUT_PATH>/multimm.log``.
    """
    import os

    root = logging.getLogger()
    effective_level = logging.DEBUG if debug else level

    # keep third-party library chatter (matplotlib, etc.) out of our output —
    # see _NOISY_THIRD_PARTY_LOGGERS
    _quiet_third_party_loggers()

    # ── Console handler (CozyFormatter) ──────────────────────────────────────
    has_cozy = any(
        isinstance(getattr(h, "formatter", None), CozyFormatter)
        for h in root.handlers
    )
    if not has_cozy:
        # Remove any handlers added by third-party libraries before us.
        for h in list(root.handlers):
            root.removeHandler(h)
            h.close()

        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(CozyFormatter(force_color=force_color or _COLOR_ENABLED))
        root.setLevel(effective_level)
        root.addHandler(stream_handler)

    # ── File handler (plain text, no ANSI) ───────────────────────────────────
    if log_file is not None:
        # Don't add a duplicate file handler pointing at the same path.
        existing_paths = {
            getattr(h, "baseFilename", None)
            for h in root.handlers
            if isinstance(h, logging.FileHandler)
        }
        abs_log = os.path.abspath(log_file)
        if abs_log not in existing_paths:
            os.makedirs(os.path.dirname(abs_log) or ".", exist_ok=True)
            file_handler = logging.FileHandler(abs_log, mode="a", encoding="utf-8")
            file_handler.setLevel(effective_level)
            file_handler.setFormatter(
                logging.Formatter(
                    fmt="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                )
            )
            root.addHandler(file_handler)
            root.info("Log file: %s", abs_log)


# ── Progress bar ──────────────────────────────────────────────────────────────

def _format_eta(seconds: float) -> str:
    """Format a duration in seconds as H:MM:SS (or MM:SS under an hour)."""
    import math
    if not math.isfinite(seconds) or seconds < 0:
        return "—"
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


class ProgressLogger:
    """A single in-place progress line — deliberately not tqdm, since tqdm's
    bare '\\r' redraw interleaves badly with this project's newline-per-record
    logger. On a TTY, writes the bar directly to the console with '\\r'
    (the logger is used only for a start line and the final 100% line);
    otherwise falls back to throttled logger lines (every ``every_pct``%).

    Shared across the project: ``from .logger import ProgressLogger``.
    """

    def __init__(self, n: int, log, label: str, every_pct: float = 5.0, bar_width: int = 24):
        self.n = max(1, int(n))
        self.log = log
        self.label = label
        self.every_pct = every_pct
        self.bar_width = bar_width
        self.start = time.time()
        self._last_pct = -1
        self._last_line_len = 0
        # Prefer stdout (where the console logger writes), fall back to
        # stderr, so the bar and the logger's own lines never fight over
        # which stream "owns" the cursor.
        self._stream = (
            sys.stdout if getattr(sys.stdout, "isatty", lambda: False)()
            else (sys.stderr if getattr(sys.stderr, "isatty", lambda: False)() else None)
        )
        self._started = False

    def _format_line(self, done: int, pct: float) -> str:
        elapsed = time.time() - self.start
        rate = done / elapsed if elapsed > 0 else 0.0
        eta = (self.n - done) / rate if rate > 0 else float("nan")
        filled = int(round(self.bar_width * pct / 100.0))
        bar = "█" * filled + "░" * (self.bar_width - filled)
        return (
            f"{self.label:<20} │{bar}│ {pct:3.0f}%  ({done}/{self.n})  "
            f"{rate:.1f} it/s  eta {_format_eta(eta)}"
        )

    def update(self, i: int) -> None:
        """Call once per iteration with the 0-indexed loop counter."""
        done = i + 1
        pct = 100.0 * done / self.n
        is_last = done == self.n

        if self._stream is not None:
            # In-place redraw: always refresh so the bar/ETA feel live, but
            # the *logger* (which would print a brand-new timestamped line
            # every call) is only ever touched once, at the very end.
            if not self._started:
                self.log.info("  %s: starting (%d total)…", self.label, self.n)
                self._started = True
            text = self._format_line(done, pct)
            pad = max(self._last_line_len - len(text), 0)
            self._stream.write("\r  " + text + (" " * pad))
            self._stream.flush()
            self._last_line_len = len(text)
            if is_last:
                self._stream.write("\n")
                self._stream.flush()
                self.log.info("  %s", text)
            return

        # Non-interactive fallback (piped/captured output): in-place redraw
        # can't work without a terminal, so fall back to throttled, one-
        # line-per-update logger output, capped at most once per whole
        # percentage point and never more often than every ``every_pct``%.
        step = max(1, round(self.n * self.every_pct / 100.0))
        if not is_last and (done % step != 0 or int(pct) == self._last_pct):
            return
        self._last_pct = int(pct)
        self.log.info("  %s", self._format_line(done, pct))


# ── Table helper ─────────────────────────────────────────────────────────────

def log_table(
    rows,
    title: str = "",
    log_fn=None,
    width: int = 62,
) -> None:
    """Emit a compact, aligned key-value table through *log_fn*.

    The table auto-expands beyond *width* when any row's content would
    otherwise overflow the right border — the caller's *width* is treated
    as a *minimum*, not a cap.

    Parameters
    ----------
    rows    : iterable of (label, value) pairs or bare strings (section
              separators / headers).
    title   : optional title embedded in the top border.
    log_fn  : callable (e.g. ``logger.info``).  Defaults to ``print``.
    width   : minimum total visual width including border characters.

    Example
    -------
    >>> log_table(
    ...     [("Mode", "svd"), ("K", 3), ("k_scale", "200.0 kJ/mol")],
    ...     title="Hi-C Force", log_fn=logger.info,
    ... )
    """
    if log_fn is None:
        log_fn = print

    # Materialise once so we can scan twice (once for widths, once for output)
    rows = list(rows)

    tuple_rows = [r for r in rows if isinstance(r, tuple) and len(r) == 2]
    str_rows   = [r for r in rows if isinstance(r, str)]

    # Label column: natural width, hard-capped at 40 to avoid eating the line
    label_w = max((len(str(r[0])) for r in tuple_rows), default=1)
    label_w = min(label_w, 40)

    # ── Auto-expand: compute the minimum width that makes every row fit ────────
    # Each tuple row needs:  │  <label:<label_w>  <1-space gap><val>│
    #   = 1 + 2 + label_w + 2 + 1 + len(val) + 1  = label_w + len(val) + 7
    if tuple_rows:
        max_val_len  = max(len(str(r[1])) for r in tuple_rows)
        min_w_data   = label_w + max_val_len + 7
    else:
        min_w_data = 0

    # Each string row needs:  │  <str>│  = 1 + 2 + len(str) + 1 = len(str) + 4
    min_w_str = max((len(s) + 4 for s in str_rows), default=0)

    # Title row needs at least: ┌─ <title> ─┐  (2 border + 2 dash + 1 space each side)
    min_w_title = len(title) + 6 if title else 0

    width = max(width, min_w_data, min_w_str, min_w_title)
    inner = width - 2  # usable columns between │ characters

    # Re-apply label_w cap now that inner is final
    label_w = min(label_w, inner - 8)

    # ── top border ────────────────────────────────────────────────────────────
    if title:
        head = f"─ {title} "
        pad  = max(0, inner - len(head))
        top  = "┌" + head + "─" * pad + "┐"
    else:
        top = "┌" + "─" * inner + "┐"
    log_fn(top)

    # ── body ──────────────────────────────────────────────────────────────────
    for row in rows:
        if isinstance(row, str):
            # section header / separator — left-aligned, padded to fill inner
            log_fn(f"│  {row:<{inner - 2}}│")
            continue
        label, val = row
        val_str = str(val)
        gap     = inner - 2 - label_w - 2 - len(val_str)
        gap     = max(gap, 1)
        log_fn(f"│  {str(label):<{label_w}}  {' ' * gap}{val_str}│")

    # ── bottom border ─────────────────────────────────────────────────────────
    log_fn("└" + "─" * inner + "┘")


# ── Section / success helpers ─────────────────────────────────────────────────

_SECTION_W = 70   # visual width of section dividers

# ANSI strings resolved once using the cached flag so they work even after
# sys.stdout has been replaced by a Tee object.
_BLUE  = f"\033[38;5;75m\033[1m" if _COLOR_ENABLED else ""
_GREEN = f"\033[38;5;82m\033[1m" if _COLOR_ENABLED else ""
_RST   = "\033[0m"               if _COLOR_ENABLED else ""


# Cycle of 256-colour codes spanning the visible spectrum, for rainbow(). Picked
# for even hue spacing and readability on both light and dark terminals.
_RAINBOW_CYCLE = (196, 208, 220, 82, 51, 33, 93, 201)


def rainbow(text: str) -> str:
    """Bold text with each non-space character coloured a different hue,
    cycling through `_RAINBOW_CYCLE`. For one-off festive log lines (e.g. a
    banner message) — not for routine logging. Respects the same colour
    detection as the rest of this module (NO_COLOR/FORCE_COLOR/tty check),
    so it degrades to plain text when colour is off.
    """
    if not _COLOR_ENABLED:
        return text
    out = []
    i = 0
    for ch in text:
        if ch.isspace():
            out.append(ch)       # don't spend a colour code on whitespace
            continue
        out.append(fg(_RAINBOW_CYCLE[i % len(_RAINBOW_CYCLE)]) + BOLD + ch)
        i += 1
    out.append(RESET)
    return "".join(out)


def log_section(name: str) -> None:
    """Print a blue section-divider directly to stdout (bypasses formatter).

    Prints a blank line, then a full-width ruled header, then another blank
    line so the section stands out clearly from the surrounding log output.

    Example output (with colour):
        ─────────────────────────────────────────────────────────────────────
          ▸  Energy Minimization
        ─────────────────────────────────────────────────────────────────────
    """
    rule = "─" * _SECTION_W
    print(f"\n{_BLUE}{rule}{_RST}")
    print(f"{_BLUE}  ▸  {name}{_RST}")
    print(f"{_BLUE}{rule}{_RST}")


def log_success(name: str, logger: logging.Logger | None = None) -> None:
    """Log a green success message through the logger (or stdout as fallback).

    Call this immediately after a section completes without raising an
    exception.  Emits a blank line after the message to visually close
    the section.

    Example output (with colour):
        [HH:MM:SS]  ✦ INF  multimm.model         ›  ✅  Energy Minimization — completed successfully
    """
    GREEN = f"\033[38;5;82m\033[1m" if _COLOR_ENABLED else ""
    RST   = "\033[0m"               if _COLOR_ENABLED else ""
    msg   = f"{GREEN}  ✅  {name} — completed successfully{RST}"
    if logger is not None:
        logger.info(msg)
    else:
        # Fallback: emit via root logger so it always goes through the formatter.
        logging.getLogger("multimm").info(msg)
    # Blank line to visually close the section.
    print("")


# ── Demo ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    setup_logger(debug=True, force_color=True)
    log = logging.getLogger("src.multimm.model")

    log.debug("Initialising force-field parameters  EV_POWER=6  sigma=0.10 nm")
    log.info("Hi-C matrix loaded — shape=(2000, 2000), non-zero entries: 2.6%")
    log.info("SVD decomposed: top-10 components explain 11.5% of |λ| variance")
    log.warning("KR normalisation unavailable — falling back to NONE")
    log.error("read_hic_matrix() got an unexpected keyword argument 'start'")
    log.critical("Simulation diverged — NaN detected in particle positions")

    log_section("Energy Minimization")
    log.info("Running L-BFGS minimizer …")
    log_success("Energy Minimization", log)
