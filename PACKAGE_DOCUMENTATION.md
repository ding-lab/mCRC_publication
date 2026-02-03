# Package Documentation Guide

This document explains the structure and purpose of package documentation files in this repository.

## Documentation Hierarchy

### 1. **Conda Environment Files** (`envs/` directory) - **Source of Truth**
   - **Purpose**: Complete, reproducible conda environments with exact package versions
   - **Files**:
     - `seurat5_env.yml` - R analysis environment (Seurat v5, Bioconductor, Python packages)
     - `morph_env.yml` - Morphological annotation environment
     - `3d-analysis_env.yml` - 3D reconstruction and spatial analysis environment
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

### 4. **Version Tables** (Machine-readable)
   - **Purpose**: Tabular format for quick version lookup
   - **Files**:
     - `packages_versions_table.md` - R packages with versions (markdown table)
     - `python_packages_versions_table.md` - Python packages with versions (markdown table)
   - **Use when**: Quick version lookup, programmatic access
   - **Format**: Markdown tables

## Quick Reference

| Need | File to Use |
|------|-------------|
| Set up analysis environment | `envs/seurat5_env.yml` |
| Understand what packages do | `Software_packages_list.md` |
| Write methods section | `Software_packages_concise.md` |
| Quick version lookup | `packages_versions_table.md` or `python_packages_versions_table.md` |
| Complete environment spec | `envs/*.yml` files |

## File Relationships

```
envs/
├── seurat5_env.yml ──────────┐
├── morph_env.yml ────────────┼─── Source of truth (complete environments)
├── 3d-analysis_env.yml ──────┘
└── README.md ──────────────────── Environment setup instructions

Software_packages_list.md ─────── Detailed documentation (human-readable)
Software_packages_concise.md ───── Concise version (for methods)

packages_versions_table.md ─────── R packages version table
python_packages_versions_table.md ─ Python packages version table
```

## Recommendations

1. **For reproducibility**: Always use conda environment files (`envs/*.yml`)
2. **For understanding**: Read `Software_packages_list.md`
3. **For methods**: Use `Software_packages_concise.md`
4. **For quick lookup**: Use version tables

## Maintenance

- When packages are updated, regenerate conda environment files:
  ```bash
  conda env export -p /path/to/env --no-builds > envs/env_name.yml
  ```
- Update `Software_packages_list.md` and `Software_packages_concise.md` when adding new packages
- Version tables can be regenerated from the conda environments if needed
