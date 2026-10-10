# Tiny synthetic example

For an installed wheel, create the same standalone public example:

```bash
tagforge init-example --out tagforge-demo
tagforge run --config tagforge-demo/00_config/config.yaml
```

It yields nine valid reads and six deduplicated molecules. Installation and
migration details are in [README](../README.md) and [MIGRATION](MIGRATION.md).

Generate the gzipped paired FASTQ files, validate the configuration, and run:

```bash
conda env create -f environment.yml
conda activate tagforge
python -m pip install -e . --no-deps --no-build-isolation
python tests/fixtures/small_fastq/generate_example_data.py
tagforge validate-config --config configs/config.example.yaml
tagforge run --config configs/config.example.yaml
```

The example intentionally contains exact barcodes, a one-base mismatch, a
one-base forward shift, a feature mismatch, duplicate reads, and two features.

For the four paired PI-seq development samples, see the
[PI-seq fixture README](pi_seq_dev.md).
