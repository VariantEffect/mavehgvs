from fqfa.constants.iupac.rna import RNA_BASES
from mavehgvs.patterns.util import combine_patterns, remove_named_groups
from mavehgvs.patterns.dna import pos_cnr

rna_nt: str = rf"[{''.join(RNA_BASES).lower()}]"
"""str: Pattern matching any lowercase RNA base.

This does not include IUPAC ambiguity characters.
"""

rna_sub: str = rf"(?P<rna_sub>(?:(?P<position>{pos_cnr})(?P<ref>{rna_nt})>(?P<alt>{rna_nt}))|(?P<equal>=))"
"""str: Pattern matching a RNA substitution with numeric or relative-to-transcript positions.
"""

rna_del: str = rf"(?P<rna_del>(?:(?P<start>{pos_cnr})_(?P<end>{pos_cnr})del)|(?:(?P<pos>{pos_cnr})del))"
"""str: Pattern matching a RNA deletion with numeric or relative-to-transcript positions.
"""

rna_dup: str = rf"(?P<rna_dup>(?:(?P<start>{pos_cnr})_(?P<end>{pos_cnr})dup)|(?:(?P<pos>{pos_cnr})dup))"
"""str: Pattern matching a RNA duplication with numeric or relative-to-transcript positions.
"""

rna_ins: str = rf"(?P<rna_ins>(?P<start>{pos_cnr})_(?P<end>{pos_cnr})ins(?P<seq>{rna_nt}+))"
"""str: Pattern matching a RNA insertion with numeric or relative-to-transcript positions.
"""

rna_delins: str = rf"(?P<rna_delins>(?:(?P<start>{pos_cnr})_(?P<end>{pos_cnr})delins)|(?:(?P<pos>{pos_cnr})delins)(?P<seq>{rna_nt}+))"
"""str: Pattern matching a RNA deletion-insertion with numeric or relative-to-transcript positions.
"""

rna_variant: str = combine_patterns([rna_sub, rna_del, rna_dup, rna_ins, rna_delins], None)
"""str: Pattern matching any single RNA variant event.
"""

rna_single_variant: str = rf"(?P<rna>r\.{rna_variant}))"
"""str: Pattern matching any complete RNA variant, including the prefix character.
"""

rna_multi_variant: str = rf"(?P<rna_multi>r\.\[{remove_named_groups(rna_variant)}(?:;{remove_named_groups(rna_variant)}){{1,}}\])"
"""str: Pattern matching any complete RNA multi-variant, including the prefix character.

Named capture groups have been removed from the variant patterns because of non-uniqueness.
Another applications of single-variant regular expressions is needed to recover the named groups from each individual
variant in the multi-variant.
"""
