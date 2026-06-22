"""Compatibility shim for editable installs with pip versions before PEP 660."""

from setuptools import find_packages, setup


setup(
    name="tagforge",
    version="0.1.3",
    description="Streaming barcode and UMI processing for paired-end FASTQ libraries",
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    install_requires=["cutadapt>=4.6", "umi-tools>=1.1.5"],
    entry_points={"console_scripts": ["tagforge=tagforge.cli:main"]},
    extras_require={"yaml": ["PyYAML>=6"], "dev": ["pytest>=7", "pytest-cov>=4"]},
)
