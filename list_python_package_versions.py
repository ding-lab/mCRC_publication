#!/usr/bin/env python3
"""
Script to extract and list all Python packages used in this project and their versions.
"""

import sys
import subprocess

# List of Python packages used in the project
python_packages = [
    'scanpy',
    'anndata',
    'squidpy',
    'spatialdata',
    'spatialdata_io',
    'banksy',
    'harmonypy',
    'secuer',
    'umap',
    'sklearn',
    'scikit-learn',
    'numpy',
    'pandas',
    'matplotlib',
    'scipy',
    'skimage',
    'scikit-image',
    'Morph'
]

def get_package_version(pkg_name):
    """Get version of a Python package."""
    try:
        # Handle package name variations
        import_name = pkg_name
        if pkg_name == 'sklearn':
            import_name = 'sklearn'
        elif pkg_name == 'scikit-learn':
            import_name = 'sklearn'
        elif pkg_name == 'skimage':
            import_name = 'skimage'
        elif pkg_name == 'scikit-image':
            import_name = 'skimage'
        
        module = __import__(import_name)
        if hasattr(module, '__version__'):
            return module.__version__
        else:
            return "Version not found"
    except ImportError:
        return "Not installed"
    except Exception as e:
        return f"Error: {str(e)}"

print("Checking Python package versions...")
print(f"Python version: {sys.version}\n")

package_versions = {}
for pkg in python_packages:
    version = get_package_version(pkg)
    package_versions[pkg] = version
    print(f"  {pkg:30s} {version}")

# Save to file
with open('python_packages_with_versions.txt', 'w') as f:
    f.write("Package\tVersion\n")
    for pkg, version in sorted(package_versions.items()):
        f.write(f"{pkg}\t{version}\n")

print("\n\nPython package versions saved to: python_packages_with_versions.txt")

# Create markdown table
with open('python_packages_versions_table.md', 'w') as f:
    f.write("| Package | Version |\n")
    f.write("|---------|----------|\n")
    for pkg, version in sorted(package_versions.items()):
        f.write(f"| {pkg} | {version} |\n")

print("Markdown table saved to: python_packages_versions_table.md")
