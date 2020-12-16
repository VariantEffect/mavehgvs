[![Build Status](https://travis-ci.com/VariantEffect/mavehgvs.svg?branch=main)](https://travis-ci.com/VariantEffect/mavehgvs)
[![Coverage Status](https://coveralls.io/repos/github/VariantEffect/mavehgvs/badge.svg?branch=main)](https://coveralls.io/github/VariantEffect/mavehgvs?branch=main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# mavehgvs
mavehgvs is the Python reference implementation of the MAVE-HGVS variant representation standard,
a strict subset of [HGVS](http://varnomen.hgvs.org/), used primarily for clinical genomics.

## The MAVE-HGVS Standard
MAVE-HGVS is a strict subset of the [HGVS Sequence Variant Nomenclature](https://varnomen.hgvs.org/), version 20.05.
HGVS nomenclature is comprehensive and very expressive and consequently includes a lot of syntax that is not needed to
represent variants from Multiplexed Assay of Variant Effect (MAVE) data and makes the variant strings more challenging 
to parse.

While packages exist for parsing HGVS (most notably the
[biocommons hgvs package](https://github.com/biocommons/hgvs/), they are intended for use in human genetics and
rely on sequence databases and reference sequence (called "target sequence" for MAVE-HGVS), which are not always
available for or relevant for multiplexed assays.

MAVE-HGVS is an attempt to define an easy-to-parse subset of the HGVS nomenclature that captures those variants that
occur in MAVE datasets, while excluding many variant types that are unlikely to be found. Importantly, the
mavehgvs implementation does not rely on external sequence databases or identifiers.

## Supported Variants
MAVE-HGVS supports DNA, RNA, and protein variants.
MAVE-HGVS supports a subset of HGVS variants including:

* substitutions
* deletions
* duplications
* insertions
* frame shifts

Many HGVS variants are unsupported including:

* inversions
* conversions
* extensions
* changes in methylation state
* RNA fusion transcripts
* mosaicism
* chimerism
* variants with uncertain consequence
* variants in trans or unknown phase
* complex variants (e.g. translocations)

For further details, including example variants, see the specification in the package documentation.

# Installation
Install mavehgvs from pip using:

```bash
pip3 install mavehgvs
```

# Feedback
To report a problem or request a new feature with either the mavehgvs package or the MAVE-HGVS standard,
please use the GitHub issue tracker.
