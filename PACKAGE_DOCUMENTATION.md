# Package Documentation Guide

This document explains the structure and purpose of package documentation files in this repository.

## Documentation Hierarchy

### 1. **Conda Environment Files** (`envs/` directory) - **Source of Truth**
   - **Purpose**: Complete, reproducible conda environments with exact package versions
   - **Files**:
     - `seurat5_env.yml` - R analysis environment (Seurat v5, Bioconductor, Python packages)
     - `morph_env.yml` - Morphological annotation environment
     - `3d-analysis_env.yml` - 3D reconstruction and spatial analysis environment
     - `banksy_env.yml` - Banksy clustering for Xenium spatial transcriptomics
   - **Use when**: Setting up environments, ensuring reproducibility, or checking exact versions
   - **Format**: YAML (conda environment specification)

### 2. **Detailed Documentation** (`Software_packages_list.md`)
   - **Purpose**: Comprehensive documentation with package descriptions, organized by category
   - **Use when**: Understanding what each package does, finding packages by function
   - **Format**: Markdown with detailed descriptions

### 3. **Concise Documentation** (`Software_packages_concise.md`)
   - **Purpose**: Brief package list organized by category, suitable for methods section
   - **Use when**: Writing methods section, quick reference
   - **Format**: Markdown with concise lists

### 4. **Version Tables** (Deprecated)
   - **Status**: Removed - information is now available directly from conda environment YAML files
   - **Alternative**: Use `conda list -n <env_name>` or check the YAML files in `envs/` for exact versions

## Quick Reference

| Need | File to Use |
|------|-------------|
| Set up analysis environment | `envs/*.yml` files |
| Understand what packages do | `Software_packages_list.md` |
| Write methods section | `Software_packages_concise.md` |
| Quick version lookup | Check YAML files in `envs/` or use `conda list -n <env_name>` |
| Complete environment spec | `envs/*.yml` files |

## File Relationships

```
envs/
├── seurat5_env.yml ──────────┐
├── morph_env.yml ────────────┤
├── 3d-analysis_env.yml ─────┼─── Source of truth (complete environments)
├── banksy_env.yml ───────────┘
└── README.md ──────────────────── Environment setup instructions

Software_packages_list.md ─────── Detailed documentation (human-readable)
Software_packages_concise.md ───── Concise version (for methods)
```

## Recommendations

1. **For reproducibility**: Always use conda environment files (`envs/*.yml`)
2. **For understanding**: Read `Software_packages_list.md`
3. **For methods**: Use `Software_packages_concise.md`
4. **For quick lookup**: Check YAML files in `envs/` or use `conda list -n <env_name>`

## Maintenance

- When packages are updated, regenerate conda environment files:
  ```bash
  conda env export -p /path/to/env --no-builds > envs/env_name.yml
  ```
- Update `Software_packages_list.md` and `Software_packages_concise.md` when adding new packages or environments
- Package versions are always sourced from the conda environment YAML files (source of truth)
