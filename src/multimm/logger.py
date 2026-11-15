"""MultiMM logging: brand colours, four verbosity levels, pipeline stages with tqdm.

Verbosity (set with ``-v`` flags or ``set_verbosity``):
    0  (no flag)  banner, pipeline stages and progress bars only
    1  (-v)       + the main log messages
    2  (-vv)      + tables with statistics of the input and output data
    3  (-vvv)     everything, including debug messages and module names
Warnings and errors are always shown. The log file always gets everything.
"""
import logging
import math
import os
import re
import shutil
import sys
import textwrap
import threading
import time
from typing import Optional

from tqdm.auto import tqdm
from tqdm import tqdm as _tqdm_std

# ── ANSI primitives ───────────────────────────────────────────────────────────

RESET  = "\033[0m"
BOLD   = "\033[1m"
DIM    = "\033[2m"
ITALIC = "\033[3m"


def fg(n: int) -> str: return f"\033[38;5;{n}m"
def bg(n: int) -> str: return f"\033[48;5;{n}m"


# ── Colour detection ──────────────────────────────────────────────────────────

# TTY state is cached at import time, before stdout gets wrapped by a Tee.
_STDOUT_IS_TTY: bool = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()
_STDERR_IS_TTY: bool = hasattr(sys.stderr, "isatty") and sys.stderr.isatty()


def _tty_colors() -> bool:
    """True when colours are wanted (honours NO_COLOR / FORCE_COLOR / COLORTERM)."""
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR") or os.environ.get("COLORTERM"):
        return True
    if os.environ.get("TERM") in ("dumb", ""):
        return False
    return _STDERR_IS_TTY or _STDOUT_IS_TTY


_COLOR_ENABLED: bool = _tty_colors()
_TRUECOLOR: bool = os.environ.get("COLORTERM", "").lower() in ("truecolor", "24bit")


# ── Brand palette (from the MultiMM logo, lifted a bit for dark terminals) ────

BRAND = {
    "indigo":  (91, 87, 217),
    "violet":  (138, 82, 208),
    "crimson": (217, 85, 143),
    "gold":    (242, 176, 58),
}
_STOPS = [BRAND[k] for k in ("indigo", "violet", "crimson", "gold")]


def _to256(r: int, g: int, b: int) -> int:
    q = lambda v: round(v / 255 * 5)  # noqa: E731
    return 16 + 36 * q(r) + 6 * q(g) + q(b)


def rgb(c) -> str:
    """Foreground colour code for an (r, g, b) tuple or a BRAND key."""
    if not _COLOR_ENABLED:
        return ""
    r, g, b = BRAND[c] if isinstance(c, str) else c
    return f"\033[38;2;{r};{g};{b}m" if _TRUECOLOR else fg(_to256(r, g, b))


def _lerp(t: float):
    """Colour at position t in [0, 1] along indigo → violet → crimson → gold."""
    t = min(max(t, 0.0), 1.0) * (len(_STOPS) - 1)
    i = min(int(t), len(_STOPS) - 2)
    f = t - i
    a, b = _STOPS[i], _STOPS[i + 1]
    return tuple(round(a[k] + (b[k] - a[k]) * f) for k in range(3))


def gradient(text: str, t0: float = 0.0, t1: float = 1.0) -> str:
    """Colour text left to right along the logo gradient."""
    if not _COLOR_ENABLED:
        return text
    n = max(len(text) - 1, 1)
    out = [
        ch if ch == " " else rgb(_lerp(t0 + (t1 - t0) * i / n)) + ch
        for i, ch in enumerate(text)
    ]
    return BOLD + "".join(out) + RESET


def _tint(text: str, t0: float = 0.0, t1: float = 1.0, italic: bool = False) -> str:
    """Like gradient() but light (not bold), optionally italic: for body text."""
    if not _COLOR_ENABLED:
        return text
    n = max(len(text) - 1, 1)
    out = [ch if ch == " " else rgb(_lerp(t0 + (t1 - t0) * i / n)) + ch for i, ch in enumerate(text)]
    return (ITALIC if italic else "") + "".join(out) + RESET


_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def strip_ansi(s: str) -> str:
    return _ANSI_RE.sub("", s)


# ── Levels and verbosity ──────────────────────────────────────────────────────

TABLE = 15   # statistics tables      (shown from -vv)
STAGE = 25   # pipeline stage events  (always shown)
logging.addLevelName(TABLE, "TABLE")
logging.addLevelName(STAGE, "STAGE")

_MIN_LEVEL = {0: STAGE, 1: logging.INFO, 2: TABLE, 3: logging.DEBUG}
_VERBOSITY = 0


def set_verbosity(n: int) -> int:
    """Set the console verbosity (0-3) and return the value used."""
    global _VERBOSITY
    _VERBOSITY = max(0, min(3, int(n)))
    return _VERBOSITY


def get_verbosity() -> int:
    return _VERBOSITY


try:
    set_verbosity(os.environ.get("MULTIMM_VERBOSITY", 0))
except ValueError:
    pass


# ── Console formatter ─────────────────────────────────────────────────────────

_ICONS = {
    "DEBUG":    ("·", "violet",  True),
    "INFO":     ("●", "indigo",  False),
    "WARNING":  ("▲", "gold",    False),
    "ERROR":    ("✖", "crimson", False),
    "CRITICAL": ("☠", "crimson", False),
}
_PREFIX_W = 13  # width of "HH:MM:SS  ●  " for continuation lines

# Kept for backward compatibility with older imports.
LEVEL_FMT = {
    "DEBUG":    (fg(63) + DIM, "·", "DBG"),
    "INFO":     (fg(114),      "✦", "INF"),
    "WARNING":  (fg(214) + BOLD, "▲", "WRN"),
    "ERROR":    (fg(196) + BOLD, "✖", "ERR"),
    "CRITICAL": (fg(231) + bg(88) + BOLD, "☠", "CRT"),
}


class CozyFormatter(logging.Formatter):
    """Brand-coloured console lines.

        12:03:51  ●  message                      (INFO)
        12:03:51  ▲  warning in amber
          ┌─ table lines are printed bare ─┐
         ▸ [3/8]  Forces                         (stage events)
    """

    def __init__(self, force_color: bool = False) -> None:
        super().__init__()
        self._color = force_color or _COLOR_ENABLED

    def _c(self, text: str, *codes: str) -> str:
        return "".join(codes) + text + RESET if self._color and any(codes) else text

    def format(self, record: logging.LogRecord) -> str:  # noqa: A003
        msg = record.getMessage()
        if record.exc_info:
            msg += "\n" + self.formatException(record.exc_info)

        if record.levelno == STAGE:
            return self._format_stage(record, msg)
        if record.levelno == TABLE:
            return self._c("  " + msg, rgb("violet"))

        icon, colour, dim = _ICONS.get(record.levelname, ("•", "indigo", False))
        ts = self._c(time.strftime("%H:%M:%S", time.localtime(record.created)), DIM, rgb("violet"))
        badge = self._c(icon, BOLD, rgb(colour))

        mod = ""
        if _VERBOSITY >= 3:
            parts = record.name.split(".")
            mod = self._c(f"{'.'.join(parts[-2:])[:20]:<20} ", DIM, rgb("indigo"))

        if record.levelno >= logging.ERROR:
            msg = self._c(msg, BOLD, rgb("crimson"))
        elif record.levelno == logging.WARNING:
            msg = self._c(msg, rgb("gold"))
        elif dim:
            msg = self._c(msg, DIM, ITALIC)

        pad = " " * (_PREFIX_W + (21 if mod else 0))
        body = ("\n" + pad).join(msg.split("\n"))
        return f"{ts}  {badge}  {mod}{body}"

    def _format_stage(self, record: logging.LogRecord, msg: str) -> str:
        kind = getattr(record, "stage_kind", "start")
        if kind == "start":
            return "\n " + self._c("▸ ", BOLD, rgb("gold")) + self._c(msg, BOLD)
        if kind == "done":
            return "   " + self._c("✓ ", BOLD, rgb("gold")) + self._c(msg, DIM)
        if kind == "skip":
            return "   " + self._c("– " + msg, DIM)
        return "   " + self._c("✖ " + msg, BOLD, rgb("crimson"))


class PlainFormatter(logging.Formatter):
    """Log-file lines without colour codes."""

    def format(self, record: logging.LogRecord) -> str:  # noqa: A003
        return strip_ansi(super().format(record))


# ── Handler: print through tqdm so progress bars are never torn apart ─────────

class VerbosityFilter(logging.Filter):
    """Console filter: applies the current verbosity to each record."""

    def filter(self, record: logging.LogRecord) -> bool:  # noqa: A003
        if record.levelno >= logging.WARNING:
            return True
        if record.levelno < _MIN_LEVEL[_VERBOSITY]:
            return False
        if _VERBOSITY < 3 and not record.name.startswith(("multimm", "src.multimm", "__main__")):
            return False
        return True


class TqdmHandler(logging.Handler):
    def emit(self, record: logging.LogRecord) -> None:
        try:
            _tqdm_std.write(self.format(record), file=sys.stdout)
        except Exception:
            self.handleError(record)


_NOISY_THIRD_PARTY_LOGGERS = (
    "matplotlib", "PIL", "numba", "h5py", "fontTools", "urllib3", "asyncio", "trame", "vtk", "pyvista",
)


def _quiet_third_party_loggers() -> None:
    for name in _NOISY_THIRD_PARTY_LOGGERS:
        logging.getLogger(name).setLevel(logging.WARNING)


def setup_logger(
    level: int = logging.INFO,
    debug: bool = False,
    force_color: bool = False,
    log_file: Optional[str] = None,
    verbosity: Optional[int] = None,
) -> None:
    """Attach the console handler (and optionally a log file) to the root logger.

    Safe to call repeatedly. *verbosity* (0-3) overrides the current level;
    ``debug=True`` means verbosity 3. The log file always receives every
    record, whatever the console verbosity.
    """
    if debug or level <= logging.DEBUG:
        set_verbosity(3)
    if verbosity is not None:
        set_verbosity(verbosity)

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)  # filtering happens per handler
    _quiet_third_party_loggers()

    has_console = any(isinstance(getattr(h, "formatter", None), CozyFormatter) for h in root.handlers)
    if not has_console:
        for h in list(root.handlers):
            root.removeHandler(h)
            h.close()
        console = TqdmHandler()
        console.setFormatter(CozyFormatter(force_color=force_color or _COLOR_ENABLED))
        console.addFilter(VerbosityFilter())
        root.addHandler(console)

    if log_file is not None:
        abs_log = os.path.abspath(log_file)
        existing = {getattr(h, "baseFilename", None) for h in root.handlers if isinstance(h, logging.FileHandler)}
        if abs_log not in existing:
            os.makedirs(os.path.dirname(abs_log) or ".", exist_ok=True)
            fh = logging.FileHandler(abs_log, mode="a", encoding="utf-8")
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(PlainFormatter(
                fmt="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
            root.addHandler(fh)
            logging.getLogger("multimm").debug("Log file: %s", abs_log)


# ── Output helpers ────────────────────────────────────────────────────────────

def _emit(text: str = "") -> None:
    """Print a line without breaking progress bars."""
    _tqdm_std.write(text, file=sys.stdout)


def _is_tty() -> bool:
    return hasattr(sys.stdout, "isatty") and sys.stdout.isatty()


def _format_eta(seconds: float) -> str:
    if not math.isfinite(seconds) or seconds < 0:
        return "—"
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


def _fmt_seconds(s: float) -> str:
    return f"{s:.1f} s" if s < 60 else _format_eta(s)


# ── Pipeline stages ───────────────────────────────────────────────────────────

STAGES = (
    "Input data", "Initial structure", "Forces", "Energy minimization",
    "MD simulation", "Analysis", "Plotting", "Outputs",
)
_stage_log = logging.getLogger("multimm.pipeline")
_ACTIVE: "Optional[Stage]" = None
_SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"


def _stage_event(msg: str, kind: str) -> None:
    _stage_log.log(STAGE, msg, extra={"stage_kind": kind})


class Stage:
    """One step of the pipeline: a header, a tqdm bar and a closing line.

    ``total`` items  -> a progress bar. Call ``next(label)`` when an item starts
                        (it completes the previous one) and ``finish()`` at the end.
    ``total=None``   -> a spinner with the elapsed time (for opaque steps).
    Usable as a context manager.
    """

    def __init__(self, name: str, total: Optional[int] = None):
        self.name = name
        self.total = total
        self.bar = None
        self._t0 = time.perf_counter()
        self._thread = None
        self._stop = threading.Event()
        self._started = False
        self._open_item = False
        self.idx = STAGES.index(name) + 1 if name in STAGES else None

    # -- lifecycle ------------------------------------------------------------
    def start(self) -> "Stage":
        global _ACTIVE
        _ACTIVE = self
        self._started = True
        self._t0 = time.perf_counter()
        tag = f"[{self.idx}/{len(STAGES)}]  " if self.idx else ""
        _stage_event(f"{tag}{self.name}", "start")
        if _is_tty():
            if self.total:
                self.bar = tqdm(
                    total=self.total, file=sys.stdout, leave=True, mininterval=0.1,
                    bar_format="   {desc:<24} {bar:26} {n_fmt}/{total_fmt}  {elapsed}",
                    colour=self._hex(), desc=self.name,
                )
            else:
                self.bar = tqdm(total=0, file=sys.stdout, leave=False,
                                bar_format="   {desc}  {elapsed}", desc=self.name)
                self._thread = threading.Thread(target=self._spin, daemon=True)
                self._thread.start()
        return self

    def _hex(self) -> Optional[str]:
        if not _COLOR_ENABLED:
            return None
        r, g, b = _lerp(((self.idx or 1) - 1) / (len(STAGES) - 1))
        return f"#{r:02x}{g:02x}{b:02x}"

    def _spin(self) -> None:
        i = 0
        rst = RESET if _COLOR_ENABLED else ""
        while not self._stop.wait(0.15):
            if self.bar is not None:
                self.bar.set_description_str(f"{rgb('gold')}{_SPINNER[i % len(_SPINNER)]}{rst} {self.name}")
                self.bar.refresh()
            i += 1

    # -- progress -------------------------------------------------------------
    def next(self, label: str) -> None:
        """Complete the previous item (if any) and start *label*."""
        self._complete_item()
        self._open_item = True
        if self.bar is not None:
            self.bar.set_description_str(label)
            self.bar.refresh()
        else:
            _stage_log.info("%s …", label)

    def _complete_item(self) -> None:
        if self._open_item:
            self._open_item = False
            if self.bar is not None and self.total:
                self.bar.update(1)

    def tick(self, n: int = 1) -> None:
        """Advance the bar by n items (for plain loops)."""
        if self.bar is not None and self.total:
            self.bar.update(n)

    def detail(self, text: str) -> None:
        """Short note beside the bar (e.g. a file name)."""
        if self.bar is not None and self.total:
            self.bar.set_postfix_str(text, refresh=True)

    def finish(self, ok: bool = True) -> None:
        global _ACTIVE
        if not self._started:
            return
        self._started = False
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1)
        if ok:
            self._complete_item()
            if self.bar is not None and self.total and self.bar.n < self.total:
                self.bar.update(self.total - self.bar.n)
        if self.bar is not None:
            self.bar.close()
            self.bar = None
        dt = _fmt_seconds(time.perf_counter() - self._t0)
        if ok:
            _stage_event(f"{self.name}  ·  {dt}", "done")
        else:
            _stage_event(f"{self.name} failed after {dt}", "fail")
        if _ACTIVE is self:
            _ACTIVE = None

    def __enter__(self) -> "Stage":
        return self.start()

    def __exit__(self, exc_type, exc, tb) -> bool:
        self.finish(ok=exc_type is None)
        return False


def stage(name: str, total: Optional[int] = None) -> Stage:
    """Create a pipeline stage (``with stage(...) as st:`` or ``.start()`` / ``.finish()``)."""
    return Stage(name, total)


def skip_stage(name: str, reason: str = "") -> None:
    """Show a stage that is not run."""
    idx = STAGES.index(name) + 1 if name in STAGES else None
    tag = f"[{idx}/{len(STAGES)}]  " if idx else ""
    _stage_event(f"{tag}{name} skipped" + (f" ({reason})" if reason else ""), "skip")


def close_active_stage(ok: bool = False) -> None:
    """Close a dangling stage (used when the run fails)."""
    if _ACTIVE is not None:
        _ACTIVE.finish(ok=ok)


def vtqdm(iterable=None, **kw):
    """tqdm for inner loops: shown from -v on a terminal, never at the default verbosity."""
    kw.setdefault("leave", False)
    kw.setdefault("file", sys.stdout)
    kw.setdefault("disable", get_verbosity() < 1 or not _is_tty())
    return tqdm(iterable, **kw)


# ── Inner progress (used by validation) ───────────────────────────────────────

class ProgressLogger:
    """Inner progress bar: a tqdm bar from -v on a terminal, throttled log lines otherwise.

    ``update(i)`` takes the 0-indexed loop counter.
    """

    def __init__(self, n: int, log, label: str, every_pct: float = 5.0, bar_width: int = 24):
        self.n = max(1, int(n))
        self.log = log
        self.label = label
        self.every_pct = every_pct
        self._last = 0
        self._last_pct = -1
        self.bar = None
        if get_verbosity() >= 1 and _is_tty():
            self.bar = tqdm(total=self.n, desc=f"  {label}", file=sys.stdout, leave=False,
                            colour="#8a52d0" if _COLOR_ENABLED else None)

    def update(self, i: int) -> None:
        done = i + 1
        if self.bar is not None:
            self.bar.update(done - self._last)
            self._last = done
            if done >= self.n:
                self.bar.close()
                self.bar = None
            return
        pct = 100.0 * done / self.n
        step = max(1, round(self.n * self.every_pct / 100.0))
        if done != self.n and (done % step != 0 or int(pct) == self._last_pct):
            return
        self._last_pct = int(pct)
        self.log.info("  %s: %3.0f%%  (%d/%d)", self.label, pct, done, self.n)


# ── Tables ────────────────────────────────────────────────────────────────────

def log_table(rows, title: str = "", log_fn=None, width: int = 62) -> None:
    """Emit an aligned key-value table (shown from -vv).

    *rows* are ``(label, value)`` pairs or bare strings (section headers).
    The table widens to fit its content; *width* is only a minimum.
    If *log_fn* is a logger method, the lines are sent at the TABLE level.
    """
    owner = getattr(log_fn, "__self__", None)
    if isinstance(owner, logging.Logger):
        emit = lambda s: owner.log(TABLE, s)  # noqa: E731
    else:
        emit = log_fn or print

    rows = list(rows)
    tuple_rows = [r for r in rows if isinstance(r, tuple) and len(r) == 2]
    str_rows = [r for r in rows if isinstance(r, str)]

    label_w = min(max((len(str(r[0])) for r in tuple_rows), default=1), 40)
    min_w_data = label_w + max((len(str(r[1])) for r in tuple_rows), default=0) + 7 if tuple_rows else 0
    min_w_str = max((len(s) + 4 for s in str_rows), default=0)
    min_w_title = len(title) + 6 if title else 0
    width = max(width, min_w_data, min_w_str, min_w_title)
    inner = width - 2
    label_w = min(label_w, inner - 8)

    if title:
        head = f"─ {title} "
        emit("┌" + head + "─" * max(0, inner - len(head)) + "┐")
    else:
        emit("┌" + "─" * inner + "┐")
    for row in rows:
        if isinstance(row, str):
            emit(f"│  {row:<{inner - 2}}│")
            continue
        label, val = row
        val_str = str(val)
        gap = max(inner - 2 - label_w - 2 - len(val_str), 1)
        emit(f"│  {str(label):<{label_w}}  {' ' * gap}{val_str}│")
    emit("└" + "─" * inner + "┘")


# ── Section / success helpers (shown from -v) ─────────────────────────────────

_SECTION_W = 70


def log_section(name: str) -> None:
    """A ruled sub-heading inside a stage (from -v)."""
    if get_verbosity() < 1:
        return
    rule = "─" * _SECTION_W
    c, r = rgb("violet"), RESET if _COLOR_ENABLED else ""
    _emit(f"\n{c}{rule}{r}")
    _emit(f"{c}  ▸  {name}{r}")
    _emit(f"{c}{rule}{r}")


def log_success(name: str, logger: Optional[logging.Logger] = None) -> None:
    """A closing line for a sub-section (from -v)."""
    if get_verbosity() < 1:
        return
    (logger or logging.getLogger("multimm")).info(f"{name} — done")


def rainbow(text: str) -> str:
    """Text coloured along the logo gradient (for one-off festive lines)."""
    return gradient(text) if _COLOR_ENABLED else text


# ── Banner ────────────────────────────────────────────────────────────────────

_LETTERS = (
    "███╗   ███╗██╗   ██╗██╗  ████████╗██╗███╗   ███╗███╗   ███╗",
    "████╗ ████║██║   ██║██║  ╚══██╔══╝██║████╗ ████║████╗ ████║",
    "██╔████╔██║██║   ██║██║     ██║   ██║██╔████╔██║██╔████╔██║",
    "██║╚██╔╝██║██║   ██║██║     ██║   ██║██║╚██╔╝██║██║╚██╔╝██║",
    "██║ ╚═╝ ██║╚██████╔╝███████╗██║   ██║██║ ╚═╝ ██║██║ ╚═╝ ██║",
    "╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝╚═╝     ╚═╝╚═╝     ╚═╝",
)
_SCALES = ("nucleosomes", "loops", "TADs", "compartments", "chromosomes", "nucleus")


def _hilbert_d2xy(order: int, d: int):
    """Point d of a 2D Hilbert curve on a 2**order grid."""
    n = 2 ** order
    x = y = 0
    t = d
    s = 1
    while s < n:
        rx = 1 & (t // 2)
        ry = 1 & (t ^ rx)
        if ry == 0:
            if rx == 1:
                x, y = s - 1 - x, s - 1 - y
            x, y = y, x
        x += s * rx
        y += s * ry
        t //= 4
        s *= 2
    return x, y


def _hilbert_art(order: int = 2):
    """The Hilbert curve (the logo's symbol) drawn with box characters.

    Returns rows of cells; a cell is ``(char, t)`` with t the position along
    the curve in [0, 1], or ``None`` for blank.
    """
    n = 2 ** order
    pts = [_hilbert_d2xy(order, d) for d in range(n * n)]
    H, W = 2 * n - 1, 3 * (n - 1) + 1
    grid = [[None] * W for _ in range(H)]
    dirs = {p: set() for p in pts}
    total = len(pts) - 1
    for k in range(total):
        (x0, y0), (x1, y1) = pts[k], pts[k + 1]
        t = (k + 0.5) / total
        if y0 == y1:   # horizontal link
            dirs[(x0, y0)].add("R" if x1 > x0 else "L")
            dirs[(x1, y1)].add("L" if x1 > x0 else "R")
            for c in range(min(x0, x1) * 3 + 1, max(x0, x1) * 3):
                grid[2 * y0][c] = ("─", t)
        else:          # vertical link
            dirs[(x0, y0)].add("D" if y1 > y0 else "U")
            dirs[(x1, y1)].add("U" if y1 > y0 else "D")
            grid[min(y0, y1) * 2 + 1][3 * x0] = ("│", t)
    glyph = {frozenset("LR"): "─", frozenset("UD"): "│", frozenset("RD"): "┌",
             frozenset("LD"): "┐", frozenset("RU"): "└", frozenset("LU"): "┘"}
    for k, (x, y) in enumerate(pts):
        d = frozenset(dirs[(x, y)])
        grid[2 * y][3 * x] = ("●" if len(d) == 1 else glyph.get(d, "●"), k / total)
    return grid


def print_banner(version: Optional[str] = None) -> None:
    """Print the startup banner: the Hilbert-curve symbol beside the logo text."""
    if version is None:
        try:
            from . import __version__ as version
        except Exception:
            version = ""

    cols = shutil.get_terminal_size((80, 24)).columns
    sym = _hilbert_art(2)
    sym_w = len(sym[0])
    show_sym = cols >= len(_LETTERS[0]) + sym_w + 6
    text_w = len(_LETTERS[0])
    rst = RESET if _COLOR_ENABLED else ""
    d = DIM if _COLOR_ENABLED else ""

    n_rows = max(len(sym), len(_LETTERS) + 1)
    top_off = (n_rows - len(_LETTERS) + 1) // 2   # vertically centre the letters
    lines = [""]
    for r in range(n_rows):
        left = ""
        if show_sym:
            left = "".join(
                " " if cell is None else gradient(cell[0], cell[1], cell[1])
                for cell in (sym[r] if r < len(sym) else [None] * sym_w)
            )
        li = r - top_off
        right = ""
        if 0 <= li < len(_LETTERS):
            right = "".join(
                " " if ch == " " else rgb(_lerp((c + 3 * li) / (text_w + 18))) + BOLD + ch
                for c, ch in enumerate(_LETTERS[li])
            ) + rst
        lines.append("  " + (left + "   " if show_sym else "") + right)

    pad = "  " + (" " * (sym_w + 3) if show_sym else "")
    beads = "".join(gradient("●─" if i < 15 else "●", i / 15, i / 15) for i in range(16))
    lines += [
        "",
        pad + gradient(" ┈ ".join(_SCALES)),
        pad + f"{d}multiscale modelling of chromatin{rst}" + (f"   {gradient('v' + version)}" if version else ""),
        pad + beads,
        "",
        f"  {d}Creator{rst}      Sebastian Korsak  (Warsaw, Poland)",
        f"  {d}Nucleosomes{rst}  Krzysztof Banecki    {d}Web server{rst}  Patryk Prusak",
        f"  {d}Contact{rst}      s.korsak@datascience.edu.pl · d.plewczynski@datascience.edu.pl",
        "",
    ]
    # one random chromatin fact: gradient rule, italic gradient text, closing mark
    from .facts import random_fact

    width = max(40, min(cols, 100) - 10)
    wrapped = textwrap.wrap(random_fact(), width)
    rule = "─" * min(width + 2, 60)
    lines.append("  " + _tint("◆ did you know?  " + rule[:max(0, len(rule) - 17)], italic=True))
    for k, ln in enumerate(wrapped):
        t0, t1 = k / len(wrapped), (k + 1) / len(wrapped)
        mark = "❝ " if k == 0 else "  "
        end = " ❞" if k == len(wrapped) - 1 else ""
        lines.append("  " + _tint(mark + ln + end, t0, t1, italic=True))
    lines.append("")
    for line in lines:
        _emit(line)


# ── Demo ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import random

    v = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    setup_logger(verbosity=v)
    log = logging.getLogger("multimm.demo")
    print_banner("2.1.0")
    with stage("Input data", 3) as s:
        for name in ("Loops", "Compartments", "Hi-C"):
            s.next(name)
            log.info("Loaded %s", name)
            log_table([("File", f"{name}.bed"), ("Rows", random.randint(100, 9999))], title=name, log_fn=log.info)
            time.sleep(0.4)
    with stage("Energy minimization"):
        log.debug("L-BFGS running")
        time.sleep(1.2)
    log.warning("A warning looks like this")
    skip_stage("MD simulation", "SIM_RUN_MD=False")
