# Installation

Getting started takes a few minutes and a few commands. If you have never installed a Python package, do not worry: just copy the lines below into a terminal one by one.

## Requirements

- **Python 3.10 or newer.**
- **Linux** is the primary, tested platform (Ubuntu, Debian, Red Hat). macOS also works, but without CUDA. Windows is not recommended (the `pyBigWig` dependency does not support it).
- A **GPU is optional but strongly recommended** for anything larger than a few thousand beads. MultiMM uses OpenMM, which supports NVIDIA GPUs through CUDA and many GPUs through OpenCL. Without a GPU the CPU platform is used.

## Install from PyPI

```bash
pip install MultiMM
```

Reading `.hic` Hi-C files needs one optional package, `hic-straw`, which has to compile some C++ code and can be fussy to install. It is therefore not installed by default; if you use `.hic` files, add it with

```bash
pip install "MultiMM[hic]"
```

`.cool` and `.mcool` files need nothing extra.

Using a virtual environment keeps the many scientific dependencies separate from your system Python:

```bash
python -m venv multimm-env
source multimm-env/bin/activate
pip install MultiMM
```

Check that the command-line tool is available:

```bash
MultiMM --help
```

## Install from source

Clone the repository and install it in editable mode if you want to modify the code:

```bash
git clone https://github.com/SFGLab/MultiMM.git
cd MultiMM
pip install -e .
```

The project also ships a `uv.lock` file, so [uv](https://docs.astral.sh/uv/) users can instead run `uv sync` inside the cloned folder to get exactly the dependency versions the project was developed with.

## Choosing a compute platform

OpenMM selects the hardware through the `PLATFORM` setting in the config file (see [Arguments](../reference/arguments.md)):

| `PLATFORM` | Hardware | When to use |
|---|---|---|
| `CUDA` | NVIDIA GPU | Fastest option if you have an NVIDIA card and drivers |
| `OpenCL` | Most GPUs (and CPUs) | Good general choice on GPUs without CUDA |
| `CPU` | CPU cores | Small systems; set `CPU_THREADS` to control parallelism |
| `Reference` | Single CPU thread | Debugging only, very slow |

If the requested GPU platform is not found, MultiMM prints a message and falls back to the CPU. To see which platforms your OpenMM installation provides, run:

```bash
python -m openmm.testInstallation
```

## Example data

A small example dataset (GM12878 Hi-C from Rao et al., CALDER subcompartments and ENCODE ATAC-seq) is available [on Google Drive](https://drive.google.com/drive/folders/1nFAPE4pCaHpeL5nw6nq0VvfUFoc24aXm?usp=sharing). The `examples/` folder of the repository holds ready-to-edit configuration files for the common scenarios (genome-wide, a single region, single-cell data, and the Hi-C force).

## Installing the documentation tools

To build this documentation locally:

```bash
pip install mkdocs-material
mkdocs serve        # live preview at http://127.0.0.1:8000
```
