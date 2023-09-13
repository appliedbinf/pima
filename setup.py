#!/usr/bin/env python
from setuptools import setup
from pathlib    import Path

base_dir     = Path(__file__).parent.resolve()
readme_file  = base_dir / "README.md"
version_file = base_dir / "VERSION"

# Get the long description from the README file
with readme_file.open(encoding = "utf-8") as f:
    long_description = f.read()

with open(version_file, "r") as version_fp:
    version = version_fp.read().strip()

setup(
    name="pima",
    version=version,
    python_requires='>=3.8',
    scripts = [
        # 'pima.py',
        # 'pima_circ.py',
        # 'sam2psl.py',
        'pima_install.sh',
        'pima_doc.py',
        'InstallScript.sh'
    ],
    install_requires=[
        "colorama; platform_system == 'Linux'",
        "importlib-metadata; python_version <= '3.8'",
    ],
    package_dir={"pima":"./"},
    package_data={'pima': ['VERSION', 'README.md', 'data/*.*', 'src/pChunks.R']},
    # py_modules = ['src.building_pycircos_figures', 'src.MarkdownReport'],
    # data_files = [('data', ['data/*.*']),('src', ['src/pChunks.R'])],
    entry_points = {
        "console_scripts": [
            "pima = pima.pima:main"
        ],
    },
    long_description = long_description,
    long_description_content_type = "text/markdown",
    license = "MIT",
    project_urls = {
        "Source":      "https://github.com/appliedbinf/MergedPima",
    },
)
