"""Generate docs/reference/arguments.md from src/multimm/config.py (AST only)."""
import ast, re, sys
from pathlib import Path

CFG = Path(sys.argv[1] if len(sys.argv) > 1 else "src/multimm/config.py")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/reference/arguments.md")

# section title -> (intro, keys)
SECTIONS = [
    ("Run and hardware", "Where and how the computation runs.",
     "PLATFORM CPU_THREADS DEVICE"),
    ("Modelling setup", "What is modelled and where results go.",
     "MODELLING_LEVEL N_BEADS CHROM LOC_START LOC_END OUT_PATH SAVE_PLOTS FORCEFIELD_PATH"),
    ("Input data", "Optional experimental inputs ([details](input.md)).",
     "LOOPS_PATH COMPARTMENT_PATH HIC_PATH ATACSEQ_PATH DOWNSAMPLING_PROB COMPARTMENT_FLIP_PROB COMPARTMENT_NOISE_STD"),
    ("Initial structure", "Starting conformation ([model](../model/initial-structure.md)).",
     "BUILD_INITIAL_STRUCTURE INITIAL_STRUCTURE_TYPE INITIAL_STRUCTURE_PATH"),
    ("Ensembles and randomness", "Many structures, shuffling and seeds.",
     "GENERATE_ENSEMBLE N_ENSEMBLE SHUFFLE_CHROMS SHUFFLING_SEED"),
    ("Backbone", "Bonds, angles and excluded volume ([model](../model/backbone.md)).",
     "POL_USE_HARMONIC_BOND POL_HARMONIC_BOND_R0 POL_HARMONIC_BOND_K POL_USE_HARMONIC_ANGLE "
     "POL_HARMONIC_ANGLE_R0 POL_HARMONIC_ANGLE_CONSTANT_K EV_USE_EXCLUDED_VOLUME EV_FORCE_TYPE "
     "EV_EPSILON EV_R_SMALL EV_POWER"),
    ("Loops", "Loop bonds from a `.bedpe` file ([model](../model/loops.md)).",
     "LE_USE_HARMONIC_BOND LE_LOOP_FORCE_TYPE LE_FIXED_DISTANCES LE_HARMONIC_BOND_R0 LE_HARMONIC_BOND_K"),
    ("Compartments and chromosomes", "A/B, subcompartment and chromosome-territory forces ([model](../model/compartments.md)).",
     "COB_USE_COMPARTMENT_BLOCKS COB_FORCE_TYPE COB_DISTANCE COB_EA COB_EB COB_CUTOFF "
     "SCB_USE_SUBCOMPARTMENT_BLOCKS SCB_FORCE_TYPE SCB_DISTANCE SCB_EA1 SCB_EA2 SCB_EB1 SCB_EB2 "
     "CHB_USE_CHROMOSOMAL_BLOCKS CHB_FORCE_TYPE CHB_KC CHB_DE"),
    ("Hi-C force", "Restraints from a Hi-C map ([model](../model/hic-force.md)).",
     "HIC_USE_FORCE HIC_NORMALIZATION HIC_K_SCALE HIC_RC HIC_AUTO_SCALE HIC_BOLTZMANN_KERNEL "
     "HIC_BOLTZMANN_ALPHA HIC_FORCE_OE HIC_NOISE_INTENSITY HIC_MAX_GAP HIC_INSULATION_WINDOW "
     "HIC_BLOCK_COPOLYMER HIC_BLOCK_COPOLYMER_MIN_BP"),
    ("Nuclear organisation", "Container, lamina and nucleolus ([model](../model/lamina.md)).",
     "SC_USE_SPHERICAL_CONTAINER SC_SCALE SC_RADIUS1 SC_RADIUS2 IBL_USE_B_LAMINA_INTERACTION "
     "BLAMINA_FORCE_TYPE IBL_SCALE CF_USE_CENTRAL_FORCE CENTRAL_FORCE_TYPE CF_STRENGTH"),
    ("Gene mode", "Modelling around a gene ([advanced](../model/advanced.md)).",
     "GENE_TSV GENE_NAME GENE_ID GENE_WINDOW"),
    ("Nucleosome interpolation", "Post-processing with nucleosomes ([advanced](../model/advanced.md)).",
     "NUC_DO_INTERPOLATION MAX_NUCS_PER_BEAD NUC_RADIUS POINTS_PER_NUC PHI_NORM"),
    ("Molecular dynamics", "Optional MD after minimisation.",
     "SIM_RUN_MD SIM_N_STEPS TRJ_FRAMES SIM_SAMPLING_STEP SIM_INTEGRATOR_TYPE SIM_INTEGRATOR_STEP "
     "SIM_FRICTION_COEFF SIM_TEMPERATURE SIM_SET_INITIAL_VELOCITIES SIM_ERROR_TOLERANCE SIM_AMD_ALPHA SIM_AMD_E"),
]
# fields accepted by the config but never read elsewhere in src/
UNUSED = {"DEVICE", "SC_RADIUS1", "SC_RADIUS2", "COB_DISTANCE", "SCB_DISTANCE"}

src = CFG.read_text()
tree = ast.parse(src)
cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "SimulationConfig")
F = {}
for n in cls.body:
    if not isinstance(n, ast.AnnAssign):
        continue
    d = desc = None
    for k in getattr(n.value, "keywords", []):
        if k.arg == "default":
            d = ast.get_source_segment(src, k.value)
        elif k.arg == "description":
            desc = ast.literal_eval(k.value) if isinstance(k.value, (ast.Constant, ast.JoinedStr)) is False or True else None
    F[n.target.id] = (ast.get_source_segment(src, n.annotation), d, desc)

def esc(s):
    return " ".join(str(s).split()).replace("|", "\\|")

def short(desc, limit=110):
    """One-line summary: first sentence, long parentheticals and asides dropped, options kept."""
    d = " ".join(str(desc).split())
    m = re.search(r"Options?:(.*?)(?:\.\s|\.$|$)", d)
    opts = ""
    if m:
        names = m.group(1)
        while re.search(r"\([^()]*\)", names):  # drop (nested) parentheses
            names = re.sub(r"\s*\([^()]*\)", "", names)
        names = re.sub(r"\s*[—-]+\s.*?(?=,|$)", "", names).strip().rstrip(".")
        opts = f" Options: {names}."
        d = d[:m.start()].strip()
    d = re.sub(r"\s*\([^()]{26,}\)", "", d)               # long parentheticals
    d = re.split(r"(?<=[.!?])\s+(?=[A-Z])|\s+[—;]\s+|(?<=[a-z-]):\s", d)[0]  # first sentence / aside
    d = d[:1].upper() + d[1:]
    if len(d) > limit:
        d = d[:limit].rsplit(" ", 1)[0].rstrip(",;:") + " …"
    d = d.rstrip()
    if d and d[-1] not in ".!?…":
        d += "."
    return esc(d + opts)


def fmt_default(d):
    if d is None:
        return "—"
    d = d.strip()
    if d in ("None", '""', "''"):
        return "*none*"
    if d.startswith("InitialStructureType."):
        return f"`{d.split('.')[-1].lower()}`"
    if d in ("default_xml_path", "default_gene_path"):
        return "*bundled file*"
    return "`" + d.strip('"').strip("'") + "`"

def fmt_type(t):
    t = t.replace("Optional[", "").replace("]", "").replace("OpenMMQuantity", "quantity").replace("Boolean", "bool").replace("ChromStr", "str").replace("InitialStructureType", "enum")
    return t

used = set()
lines = ["# All arguments", "",
         "*This page is generated from `src/multimm/config.py` by `docs_tools/gen_arguments.py`, so names, types and defaults always match the code. Descriptions are shortened here; `MultiMM --help` and the model pages give the full text.*", "",
         "Types: `quantity` is a number with a unit written as text, for example `0.1 nanometer` or `310 kelvin` ([syntax](../model/scales.md)); `bool` accepts `True/False`, `yes/no`, `1/0`. "
         "\"*none*\" means the argument is empty by default. Arguments marked **unused** are accepted but not read anywhere in the current code.", ""]
for title, intro, keys in SECTIONS:
    lines += [f"## {title}", "", intro, "", "| Argument | Type | Default | Description |", "|---|---|---|---|"]
    for k in keys.split():
        t, d, desc = F[k]
        used.add(k)
        flag = " **(unused)**" if k in UNUSED else ""
        lines.append(f"| `{k}` | {fmt_type(t)} | {fmt_default(d)} | {short(desc)}{flag} |")
    lines.append("")
missing = set(F) - used
assert not missing, f"fields without a section: {missing}"
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("\n".join(lines))
print(f"{len(used)} arguments written to {OUT}")
