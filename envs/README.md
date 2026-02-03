# Conda Environment Files

This directory contains conda environment YAML files required for the mCRC manuscript analysis. These YAML files are the **source of truth** for reproducible environment setup and contain exact package versions used in the analysis.

For detailed package documentation, see:
- `../PACKAGE_DOCUMENTATION.md` - **Start here**: Guide to all package documentation files
- `../Software_packages_list.md` - Complete package documentation with descriptions
- `../Software_packages_concise.md` - Concise version for methods section

## Available Environments

### 1. `seurat5_env.yml`
**Purpose**: R analysis environment for Seurat v5 single-cell and spatial transcriptomics analysis

**Usage**:
```bash
conda env create -f seurat5_env.yml
conda activate seurat5_env
```

**Build time**: ~30-60 minutes (depending on network speed and system)

**Original location**: `/diskmnt/Users2/epeng/tools/conda_envs/seurat5_env`

**Key packages**:
- R 4.4.3
- Seurat 5.x
- Bioconductor packages
- JupyterLab and R kernel
- Python 3.13.3 with scanpy, pandas, numpy

### 2. `morph_env.yml`
**Purpose**: Morphological annotation environment

**Usage**:
```bash
conda env create -f morph_env.yml
conda activate morph_env
# Install morph from GitHub
pip install git+https://github.com/ding-lab/morph.git
```

**Build time**: ~5-10 minutes

**Original location**: `/diskmnt/Users2/epeng/tools/conda_envs/morph_env`

**Key packages**:
- Python 3.10.19
- Minimal dependencies for morphological annotation workflows
- **Morph** - Must be installed separately from GitHub: https://github.com/ding-lab/morph

### 3. `3d-analysis_env.yml`
**Purpose**: 3D reconstruction and spatial analysis environment

**Usage**:
```bash
conda env create -f 3d-analysis_env.yml
conda activate 3d-analysis
```

**Build time**: ~15-30 minutes

**Original location**: `/diskmnt/Users2/epeng/tools/conda_envs/3d-analysis`

**Key packages**:
- Python 3.13.7
- scanpy 1.11.4
- squidpy 1.2.2
- scikit-image 0.25.2
- matplotlib, seaborn
- JupyterLab with extensions
- Image processing libraries (imageio, tifffile, pillow)

### 4. `banksy_env.yml`
**Purpose**: Banksy clustering for Xenium spatial transcriptomics

**Usage**:
```bash
conda env create -f banksy_env.yml
conda activate banksy
```

**Build time**: ~20-40 minutes

**Original location**: `/diskmnt/Users2/simonmo/Software/miniforge3/envs/banksy`

**Key packages**:
- Python 3.11.6
- scanpy 1.9.5
- squidpy 1.6.5
- spatialdata 0.4.0
- harmonypy 0.0.10
- secuer 1.1 (via pip)
- leidenalg 0.10.1
- JupyterLab and notebooks

## Package Versions

The exact versions of all packages (including dependencies) are specified in each YAML file. Key package versions:

**seurat5_env.yml:**
- R 4.4.3
- Seurat 5.3.0
- Signac 1.14.0
- Python 3.13.3
- See the YAML file for complete dependency list

**morph_env.yml:**
- Python 3.10.19
- Minimal dependencies for morphological annotation
- **Note**: Morph package must be installed separately after environment creation:
  ```bash
  pip install git+https://github.com/ding-lab/morph.git
  ```

**3d-analysis_env.yml:**
- Python 3.13.7
- scanpy 1.11.4
- squidpy 1.2.2
- scikit-image 0.25.2
- See the YAML file for complete dependency list

**banksy_env.yml:**
- Python 3.11.6
- scanpy 1.9.5
- squidpy 1.6.5
- spatialdata 0.4.0
- harmonypy 0.0.10
- secuer 1.1 (via pip)
- See the YAML file for complete dependency list

## Notes

- All environment files were exported with `--no-builds` flag to ensure portability across different platforms
- To recreate an environment from a YAML file, use: `conda env create -f <env_file>.yml`
- To update an existing environment: `conda env update -f <env_file>.yml --prune`
- The prefix paths in the YAML files point to the original installation locations but will be ignored when creating new environments
- For a human-readable list of packages, see `../Software_packages_list.md`
