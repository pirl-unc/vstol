# VSTOL
VSTOL (**V**ariant **S**tandardization, **T**abulation, and **O**perations **L**ibrary) is a toolkit that 
standardizes smalll and structural DNA variants from VCF files into the [Occam's Variant Grammar](https://pirlblog.substack.com/p/variant-standardization) 
TSV files. It also supports list operations (*annotate*, *diff*, *intersect*, *merge*, *overlap*) on variant lists.

[![build](https://github.com/pirl-unc/vstol/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/vstol/actions/workflows/main.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

## 01. Installation

```
pip install . --verbose
```

## 02. Dependencies
- python>=3.10
- pandas>=2.0.3
- numpy>=1.22.3
- networkx>=3.4.2
- intervaltree>=3.1.0

## 03. Usage

```
vstol [-h] [--version] {annotate,diff,intersect,merge,overlap,vcf2tsv}
```

## 04. Available Commands

| Command   | Description                                                                                           |
|-----------|-------------------------------------------------------------------------------------------------------|
| annotate  | Annotate variants.                                                                                    |
| diff      | Diff target variants from query variants.                                                             |
| intersect | Identify intersecting variants.                                                                       |
| merge     | Merge multiple variant lists into one list.                                                           |
| overlap   | Identify variants that overlap with a list of genomic regions.                                        |
| vcf2tsv   | Convert a VCF file (see below for supported variant callers) to the Occam's Variant Grammar TSV file. |

## 05. TSV File Headers

`vstol vcf2tsv` outputs each variant callset as a tab-separated values (TSV) file with the following columns:

* id
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

Additional columns from the VCF file will be included.

## 06. Supported Variant Callers

To use `VSTOL`, we recommend that you first convert a VCF file to a TSV file 
using the `vcf2tsv` command in `VSTOL`. The following variant callers are currently supported:

- [ClairS](https://github.com/HKU-BAL/ClairS)
- [CuteSV](https://github.com/tjiangHIT/cuteSV)
- [DeepVariant](https://github.com/google/deepvariant)
- [Delly2](https://github.com/dellytools/delly)
- [GATK4 Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
- [Lumpy](https://github.com/arq5x/lumpy-sv)
- [Manta](https://github.com/Illumina/manta)
- [PBSV](https://github.com/PacificBiosciences/pbsv)
- [Savana](https://github.com/cortes-ciriano-lab/savana)
- [Severus](https://github.com/KolmogorovLab/Severus)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- [Strelka2](https://github.com/Illumina/strelka)
- [Svim](https://github.com/eldariont/svim)
- [SVision-pro](https://github.com/songbowang125/SVision-pro)
