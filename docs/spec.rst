.. _spec-docs:

MAVE-HGVS specification
=======================

MAVE-HGVS is a strict subset of the `HGVS Sequence Variant Nomenclature <https://varnomen.hgvs.org/>`_, version 20.05.
HGVS nomenclature is comprehensive and very expressive and consequently includes a lot of syntax that is not needed to
represent variants from Multiplexed Assay of Variant Effect (MAVE) data and makes the variant strings challenging to
parse in general.

While packages exist for parsing HGVS (most notably the
`biocommons hgvs package <https://github.com/biocommons/hgvs/>`_), they are intended for use in human genetics and
rely on sequence databases and reference sequence, which are not always available for or relevant to multiplexed assay
targets.

MAVE-HGVS is an attempt to define an easy-to-parse subset of the HGVS nomenclature that captures those variants that
occur in MAVE datasets, while excluding many variant types that are unlikely to be found. Importantly, the
:ref:`corresponding implementation <api-docs>` of MAVE-HGVS does not rely on external sequence databases or identifiers.

Key differences between HGVS and MAVE-HGVS
------------------------------------------

Unlike standard HGVS strings that begin with an identifier (e.g. "NM_001130145.3:c.832C>T"), MAVE-HGVS strings are
usually encountered without this extra information (e.g. "c.832C>T"). Because all variants in a single dataset are
typically variants of the same target sequence, making this information redundant. Additionally some MAVE datasets
(such as those interrogating synthetic proteins or codon-optimized sequences) do not have a suitable identifier in a
reference database.

MAVE-HGVS supports a subset of variant types. Supported variants include:

* substitutions
* deletions
* duplications
* insertions
* frame shifts
* extensions

Many other variant types are not supported. Unsupported variants include:

* inversions
* conversions
* changes in methylation state
* RNA fusion transcripts
* mosaicism
* chimerism
* most variants describing uncertainty
* variants in trans or unknown phase
* complex variants (e.g. translocations)

Substitution
------------

# TODO

Deletion
--------

# TODO

Duplication
-----------

# TODO

Insertion
---------

MAVE-HGVS only supports simple insertions of a specified nucleotide or amino acid sequence.
Insertions of a number of unknown bases or amino acids are not supported.
Insertions using ambiguity characters (e.g. N or Xaa) are not supported.
Insertions must be specified by listing the full inserted sequence. Referring to the sequence that is inserted based on
its position in the reference sequence is not supported.

Examples of valid insertions include:

* g.234_235insT
* c.84_85insCTG
* c.99+6_99+7insA
* p.His7_Gln8insSer
* p.Ala12_Pro13insGlyProCys

Deletion-Insertion
------------------

# TODO

Example nucleotide variants
---------------------------

Example protein variants
------------------------

Unlike HGVS, MAVE-HGVS does not distinguish between protein variants that have been observed experimentally and those
that are predicted based on DNA level data.
