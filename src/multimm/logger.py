import logging
import sys
import os
import time


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
    "WARNING":  (fg(214) + BOLD,             "▲", "WRN"),
    "ERROR":    (fg(203) + BOLD,             "✖", "ERR"),
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
    # stdout first; if it's been redirected, check stderr as fallback
    if hasattr(sys.stdout, "isatty") and sys.stdout.isatty():
        return True
    if hasattr(sys.stderr, "isatty") and sys.stderr.isatty():
        return True
    return False


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
        self._color: bool = force_color or _tty_colors()

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
        # All badges are the same visible width so no extra padding is needed,
        # but we assemble plain text before wrapping in ANSI codes.
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
        elif record.levelname in ("ERROR", "CRITICAL"):
            msg_ = self._c(msg, lvl_color)
        else:
            msg_ = msg

        line = f"{ts_}  {badge}  {mod_} {sep}  {msg_}"

        # ── CRITICAL: bordered banner ──────────────────────────────────────────
        if record.levelname == "CRITICAL":
            w   = max(len(msg) + 8, 64)
            bar = self._c("━" * w, fg(196) + BOLD)
            lbl = self._c(f"  ☠  CRITICAL  ☠  {msg}  ", fg(231) + bg(88) + BOLD)
            line = f"\n{bar}\n{lbl}\n{bar}\n"

        return line


# ── Public API ────────────────────────────────────────────────────────────────

def setup_logger(
    level: int = logging.INFO,
    debug: bool = False,
    force_color: bool = False,
) -> None:
    """
    Attach CozyFormatter to the root logger.  Safe to call multiple times.

    Parameters
    ----------
    level       Minimum log level when *debug* is False (default INFO).
    debug       If True, sets level to DEBUG regardless of *level*.
    force_color Emit ANSI codes even when stdout is not a tty (e.g. piped).
    """
    root = logging.getLogger()
    if root.handlers:
        return

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(CozyFormatter(force_color=force_color))

    root.setLevel(logging.DEBUG if debug else level)
    root.addHandler(handler)


# ── Table helper ─────────────────────────────────────────────────────────────

def log_table(
    rows,
    title: str = "",
    log_fn=None,
    width: int = 62,
) -> None:
    """Emit a compact, aligned key-value table through *log_fn*.

    Parameters
    ----------
    rows    : iterable of (label, value) pairs or bare strings (section
              separators / headers).
    title   : optional title embedded in the top border.
    log_fn  : callable (e.g. ``logger.info``).  Defaults to ``print``.
    width   : total visual width of the table including border characters.

    Example
    -------
    >>> log_table(
    ...     [("Mode", "svd"), ("K", 3), ("k_scale", "200.0 kJ/mol")],
    ...     title="Hi-C Force", log_fn=logger.info,
    ... )
    """
    if log_fn is None:
        log_fn = print

    inner = width - 2  # usable columns between │ characters

    # ── top border ────────────────────────────────────────────────────────────
    if title:
        head = f"─ {title} "
        pad  = max(0, inner - len(head))
        top  = "┌" + head + "─" * pad + "┐"
    else:
        top = "┌" + "─" * inner + "┐"
    log_fn(top)

    # ── body ──────────────────────────────────────────────────────────────────
    tuple_rows = [r for r in rows if isinstance(r, tuple) and len(r) == 2]
    label_w    = max((len(str(r[0])) for r in tuple_rows), default=1)
    label_w    = min(label_w, inner - 8)  # always leave room for values

    for row in rows:
        if isinstance(row, str):
            # section header / separator
            log_fn(f"│  {row:<{inner - 2}}│")
            continue
        label, val = row
        val_str = str(val)
        gap     = inner - 2 - label_w - 2 - len(val_str)
        gap     = max(gap, 1)
        log_fn(f"│  {str(label):<{label_w}}  {' ' * gap}{val_str}│")

    # ── bottom border ─────────────────────────────────────────────────────────
    log_fn("└" + "─" * inner + "┘")


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
