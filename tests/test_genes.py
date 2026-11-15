"""Gene lookup, multi-gene regions and region annotation."""
import numpy as np
import pytest

from multimm import utils as u
from multimm.config import SimulationConfig

TSV = SimulationConfig().GENE_TSV


def test_single_gene_region():
    chrom, region, gene = u.get_gene_region(TSV, gene_name="MYC", window_size=1000)
    assert chrom == "chr8" and region == [gene[0] - 1000, gene[1] + 1000]


def test_multi_gene_region_is_union():
    _, _, g1 = u.get_gene_region(TSV, gene_name="MYC")
    _, _, g2 = u.get_gene_region(TSV, gene_name="PVT1")
    _, _, both = u.get_gene_region(TSV, gene_name="MYC, PVT1")
    assert both == [min(g1[0], g2[0]), max(g1[1], g2[1])]


def test_gene_errors():
    with pytest.raises(ValueError):
        u.get_gene_region(TSV, gene_name="NOT_A_GENE")
    with pytest.raises(ValueError):
        u.get_gene_region(TSV, gene_name="MYC,TP53")  # chr8 vs chr17
    with pytest.raises(ValueError):
        u.get_gene_region(TSV)


def test_genes_in_region_beads():
    _, region, _ = u.get_gene_region(TSV, gene_name="MYC", window_size=200000)
    g = u.genes_in_region(TSV, "8", region, 1000, targets=["MYC"])
    assert len(g) > 1 and g["is_target"].sum() == 1
    assert g[["bead_start", "bead_end"]].min().min() >= 0 and g["bead_end"].max() <= 999
    assert (g["bead_end"] >= g["bead_start"]).all()


def test_gene_geometry_columns():
    _, region, _ = u.get_gene_region(TSV, gene_name="MYC", window_size=100000)
    g = u.genes_in_region(TSV, "chr8", region, 200, targets=["MYC"])
    V = np.random.default_rng(0).random((200, 3))
    out, cen = u.annotate_gene_geometry(g, V)
    assert cen.shape == (len(g), 3) and {"rg", "dist_to_target"} <= set(out.columns)
    assert out.loc[out["is_target"], "dist_to_target"].min() < out["dist_to_target"].max()


def test_gene_plots(tmp_path):
    from multimm import plots as p
    _, region, _ = u.get_gene_region(TSV, gene_name="MYC", window_size=100000)
    g = u.genes_in_region(TSV, "chr8", region, 200, targets=["MYC"])
    V = np.cumsum(np.random.default_rng(1).normal(size=(200, 3)), 0) * 0.1
    g, cen = u.annotate_gene_geometry(g, V)
    p.plot_gene_track(g, 200, cen, save_path=str(tmp_path / "t.png"))
    p.save_chimera_genes_cmd(g, str(tmp_path / "g.cmd"))
    p.viz_gene_annotation(V, g, save_path=str(tmp_path / "s.png"))
    assert all((tmp_path / f).exists() for f in ("t.png", "g.cmd", "s.png"))


def test_gene_level_model_maps_all_genes(tmp_path):
    from multimm.model import MultiMM
    cfg = SimulationConfig(LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe", OUT_PATH=str(tmp_path),
                           MODELLING_LEVEL="GENE", GENE_NAME="MYC,PVT1", GENE_WINDOW=300000, SAVE_PLOTS=False)
    md = MultiMM(cfg)
    g = md.region_genes
    assert set(g.loc[g["is_target"], "gene_name"]) >= {"MYC", "PVT1"}
    assert 0 <= md.gene_start <= md.gene_end <= cfg.N_BEADS - 1
