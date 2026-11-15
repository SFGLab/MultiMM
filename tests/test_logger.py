"""Quick checks of the logger: verbosity filtering, tables, stages, banner."""
import logging

import pytest

from multimm import logger as L


@pytest.fixture(autouse=True)
def _reset():
    v = L.get_verbosity()
    yield
    L.set_verbosity(v)


def _rec(level, name="multimm.x"):
    return logging.LogRecord(name, level, "f", 1, "m", None, None)


@pytest.mark.parametrize("v,shown", [
    (0, {L.STAGE}), (1, {L.STAGE, logging.INFO}),
    (2, {L.STAGE, logging.INFO, L.TABLE}), (3, {L.STAGE, logging.INFO, L.TABLE, logging.DEBUG}),
])
def test_verbosity_levels(v, shown):
    L.set_verbosity(v)
    f = L.VerbosityFilter()
    for lvl in (logging.DEBUG, logging.INFO, L.TABLE, L.STAGE):
        assert f.filter(_rec(lvl)) == (lvl in shown)
    assert f.filter(_rec(logging.WARNING))  # warnings always shown


def test_verbosity_is_clamped():
    assert L.set_verbosity(9) == 3 and L.set_verbosity(-2) == 0


def test_log_table_lines():
    out = []
    L.log_table([("a", 1), "header", ("bb", "x")], title="T", log_fn=out.append)
    assert out[0].startswith("┌") and out[-1].startswith("└")
    assert len({len(x) for x in out}) == 1  # aligned box


def test_stage_runs_without_tty(capsys):
    with L.stage("Input data", 2) as st:
        st.next("a"); st.next("b")
    assert L._ACTIVE is None


def test_close_active_stage():
    st = L.stage("Forces", 1).start()
    L.close_active_stage()
    assert L._ACTIVE is None and st.bar is None


def test_banner_and_gradient(capsys):
    L.print_banner("1.2.3")
    out = L.strip_ansi(capsys.readouterr().out)
    assert "1.2.3" in out and "Sebastian Korsak" in out


def test_hilbert_art_visits_all_cells():
    art = L._hilbert_art(2)
    assert sum(c is not None for row in art for c in row) > 16


def test_progress_logger_updates(caplog):
    log = logging.getLogger("multimm.test")
    L.set_verbosity(0)
    p = L.ProgressLogger(10, log, "x", every_pct=50)
    with caplog.at_level(logging.INFO, logger="multimm.test"):
        for i in range(10):
            p.update(i)
    assert any("100%" in r.getMessage() for r in caplog.records)


def test_facts_are_unique_and_forty():
    from multimm.facts import FACTS, random_fact
    assert len(FACTS) == 40 and len(set(FACTS)) == 40
    assert random_fact() in FACTS


def test_banner_shows_a_fact(capsys):
    from multimm.facts import FACTS
    L.print_banner("1.0")
    out = " ".join(L.strip_ansi(capsys.readouterr().out).split())
    assert any(f.split()[0] in out and f.split()[-1] in out for f in FACTS)
