# Tiny example

Generate the gzipped paired FASTQ files, validate the configuration, and run:

```bash
conda env create -f environment.yml
conda activate tagforge
python -m pip install -e . --no-deps --no-build-isolation
python examples/generate_example_data.py
tagforge validate-config --config configs/config.example.yaml
tagforge run --config configs/config.example.yaml
```

The example intentionally contains exact barcodes, a one-base mismatch, a
one-base forward shift, a feature mismatch, duplicate reads, and two features.
