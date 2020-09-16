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

Standard HGVS strings have the format :code:`reference:variant` (e.g. :code:`NM_001130145.3:c.832C>T`).
MAVE-HGVS strings typically include the variant portion only and the reference portion is inferred based on the MAVE
design.

Reference identifiers in MAVE-HGVS are optional, and would typically be used in cases where a mix of MAVE datasets are
being analyzed jointly or for experimental designs that contain multiple target sequences.
Reference identifiers in MAVE-HGVS can contain any word characters, numbers, or the underscore.

MAVE-HGVS does not distinguish between variants that have been observed experimentally and the predicted consequence of
observed variants.
Therefore, variants that contain :code:`()` to denote predicted consequences are considered invalid.

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

Substitution
------------

MAVE-HGVS supports substitutions of a single nucleotide or amino acid.

Unlike in HGVS, variants that describe identity to the reference at a single position (e.g. :code:`c.44=`) are not
valid.

Variants describing identity to the full reference (e.g. :code:`c.=`) are valid and are the intended way to specify
identity to the target (wild-type) sequence.

MAVE-HGVS does not support extension variants, which extend an amino acid sequence to the N- or C- terminal end
(e.g. :code:`p.Met1ext-4` for gain of an upstream start or :code:`p.Ter345Lysext5` for a new downstream termination
codon). Variants that result in an N-terminal extension should use `Insertion`_ syntax and variants that remove a
termination codon should be written as standard substitution variants.

Substitutions of more than one base at a time are covered under `Deletion-Insertion`_.

Examples of valid substitutions include:

* g.48C>A
* c.=
* c.122-6T>A
* p.Glu27Trp
* p.Ter345Lys

Deletion
--------

MAVE-HGVS supports deletions of specified nucleotides or amino acids.

Deletions of an unknown number of bases or amino acids are not supported.
For example, deletions where the breakpoint is not known or where the deletion extends past the end of the target
cannot be represented with uncertainty.
To represent a deletion of a sequence including the start or end of the target, specify the deletion exactly as if it
extended to the first or last position.

Examples of valid deletions include:

* g.44del
* c.78+5_78+10del
* c.1_95del
* p.Gly18del
* p.Gln7_Asn19del

Duplication
-----------

MAVE-HGVS supports duplications of one or more nucleotides or amino acids.
The syntax is the same as HGVS.

Examples of valid insertions include:

* g.22_24dup
* c.77dup
* c.101+1_101+7dup
* p.Pro12_Gly18dup
* p.Cys5dup

Insertion
---------

MAVE-HGVS supports insertions of a specified nucleotide or amino acid sequence.

Insertions of a number of unspecified bases or amino acids or insertions using ambiguity characters (e.g. N or Xaa)
are not supported.

Insertions must be specified by listing the complete inserted sequence.
Referring to the sequence that is inserted based on its position in the reference sequence is not considered valid for
MAVE-HGVS.

Examples of valid insertions include:

* g.234_235insT
* c.84_85insCTG
* c.99+6_99+7insA
* p.His7_Gln8insSer
* p.Ala12_Pro13insGlyProCys

Deletion-Insertion
------------------

MAVE-HGVS supports deletion-insertions of a specified nucleotide or amino acid sequence.

Deletion-insertions of a number of unspecified bases or amino acids or insertions using ambiguity characters
(e.g. N or Xaa) are not supported. This includes deletion-insertions with uncertain breakpoints.

Examples of valid insertions include:

* g.22delinsAACG
* c.83_85delinsT
* c.43-6_595+12delinsCTT
* p.Ile71_Cys80delinsSer
* p.His44delinsValProGlyGlu
