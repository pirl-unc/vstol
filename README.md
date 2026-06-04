# VSTOL

**V**ariant **S**tandardization, **T**abulation, and **O**perations **L**ibrary

VSTOL standardizes small and structural DNA variants from VCF files into the
[Occam's Variant Grammar](https://pirlblog.substack.com/p/variant-standardization)
TSV format, then runs list operations — *annotate*, *collapse*, *diff*,
*intersect*, *merge*, and *overlap* — on the standardized variant lists.

[![CI](https://github.com/pirl-unc/vstol/actions/workflows/ci.yml/badge.svg)](https://github.com/pirl-unc/vstol/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**Documentation**: [https://pirl-unc.github.io/vstol/](https://pirl-unc.github.io/vstol/)

## 01. Installation

`vstol` is available on [PyPI](https://pypi.org/project/vstol/):

```bash
pip install vstol
```

## 02. Dependencies

* Python (>=3.10)
* pandas (>=2.0.3)
* numpy (>=1.22.3)
* networkx (>=3.4.2)
* intervaltree (>=3.1.0)

## 03. Usage

### View all available subcommands
```bash
vstol --help
```

### View a subcommand's parameters
```bash
vstol <subcommand> --help
```

### Available subcommands

| Subcommand  | Description                                                               |
|-------------|--------------------------------------------------------------------------|
| `annotate`  | Annotate variants with gene-level context.                               |
| `collapse`  | Collapse a merged list into one row per group, applying somatic filters.  |
| `diff`      | Output the variants specific to a target list.                           |
| `intersect` | Identify variants shared across variant lists.                           |
| `merge`     | Merge multiple variant lists into one.                                    |
| `overlap`   | Identify variants overlapping genomic regions.                           |
| `vcf2tsv`   | Convert a VCF file to the Occam's Variant Grammar TSV.                    |

See the [Commands documentation](https://pirl-unc.github.io/vstol/cli/) for full
parameter documentation.

## 04. Supported Variant Callers

Most workflows start with `vstol vcf2tsv`, which converts a variant caller's VCF
into the Occam's Variant Grammar TSV. Select the caller with `--method` (e.g.
`vstol vcf2tsv --method savana ...`).

See the
[Supported Variant Callers](https://pirl-unc.github.io/vstol/variant-callers.html)
documentation for the full list of callers and the versions tested, or run
`vstol vcf2tsv --help` for the authoritative list of accepted `--method` values.

## 05. Output Format (Occam's Variant Grammar)

`vstol vcf2tsv` writes each callset as a tab-separated (TSV) file with the
following columns:

* variant_id
* source_id
* sample_id
* platform
* method
* chromosome_1
* position_1
* strand_1
* operation_1
* chromosome_2
* position_2
* strand_2
* operation_2
* sequence

Additional columns from the source VCF are preserved.

## 06. License

Licensed under the [Apache License, Version 2.0](https://opensource.org/licenses/Apache-2.0).
