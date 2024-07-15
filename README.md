# VSTOL
VSTOL (**V**ariant **S**tandardization, **T**abulation, and **O**perations **L**ibrary) 
converts VCF files of small and structural variants to TSV files and facilitates 
list operations (anntoate, diff, filter, intersect, merge, overlap etc) on a 
variants list.

[![build](https://github.com/pirl-unc/vstol/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/vstol/actions/workflows/main.yml)

## 01. Installation

```
pip install . --verbose
```

## 02. Dependencies
- python>= 3.10
- pandas>=2.0.3
- numpy>=1.22.3
- maturin>=0.14,<0.15
- pysam
- rust

## 03. Usage

```
vstol [-h] [--version] {annotate,diff,filter,intersect,merge,overlap,vcf2tsv}
```

## 04. Available Commands

| Command   | Description                                                                                                                 |
|-----------|-----------------------------------------------------------------------------------------------------------------------------|
| annotate  | Annotate variant calls using [pyensembl](https://github.com/openvax/pyensembl) or [gencode](https://www.gencodegenes.org/). |
| collapse  | Collapse a variants list into unique variants.                                                                              |
| diff      | Identify variant calls specific to a list.                                                                                  |
| filter    | Filter variant calls (can be used to identify somatic variants).                                                            |
| intersect | Identify intersecting variant calls.                                                                                        |
| merge     | Merge variant calls from various variant callers.                                                                           |
| overlap   | Identify variants that overlap with a list of genomic ranges.                                                               |
| score     | Calculates average alignment score for each breakpoint.                                                                     | 
| vcf2tsv   | Convert a VCF file (see below for supported variant callers) to a TSV file.                                                 |

## 05. Supported Variant Callers

To use `VSTOL`, we recommend that you first convert a VCF file to a TSV file 
using the `vcf2tsv` command in `VSTOL`. 
The following variant callers are currently supported:

- [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [DeepVariant](https://github.com/google/deepvariant)
- [Delly2](https://github.com/dellytools/delly)
- [GATK4 Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
- [LUMPY](https://github.com/arq5x/lumpy-sv)
- [pbsv](https://github.com/PacificBiosciences/pbsv)
- [Savana](https://github.com/cortes-ciriano-lab/savana)
- [Severus](https://github.com/KolmogorovLab/Severus)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- [Strelka2](https://github.com/Illumina/strelka)
- [SVIM](https://github.com/eldariont/svim)
- [SVision-pro](https://github.com/songbowang125/SVision-pro)

