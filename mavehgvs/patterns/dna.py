from fqfa.constants import DNA_BASES
from mavehgvs.patterns.util import combine_patterns, remove_named_groups
from mavehgvs.patterns.position import pos, pos_intron, pos_intron_utr

dna_nt: str = rf"[{''.join(DNA_BASES)}]"
"""str: Pattern matching any uppercase DNA base.

This does not include IUPAC ambiguity characters.
"""

dna_sub_c: str = rf"(?P<dna_sub_c>(?:(?P<position>{pos_intron_utr})(?P<ref>{dna_nt})>(?P<new>{dna_nt}))|(?P<equal>=))"
"""str: Pattern matching a DNA substitution with numeric, intronic, or UTR positions.
"""

dna_del_c: str = rf"(?P<dna_del_c>(?:(?:(?P<start>{pos_intron_utr})_(?P<end>{pos_intron_utr}))|(?P<pos>{pos_intron_utr}))del)"
"""str: Pattern matching a DNA deletion with numeric, intronic, or UTR positions.
"""

dna_dup_c: str = rf"(?P<dna_dup_c>(?:(?:(?P<start>{pos_intron_utr})_(?P<end>{pos_intron_utr}))|(?P<pos>{pos_intron_utr}))dup)"
"""str: Pattern matching a DNA duplication with numeric, intronic, or UTR positions.
"""

dna_ins_c: str = rf"(?P<dna_ins_c>(?P<start>{pos_intron_utr})_(?P<end>{pos_intron_utr})ins(?P<seq>{dna_nt}+))"
"""str: Pattern matching a DNA insertion with numeric, intronic, or UTR positions.
"""

dna_delins_c: str = rf"(?P<dna_delins_c>(?:(?:(?P<start>{pos_intron_utr})_(?P<end>{pos_intron_utr}))|(?P<pos>{pos_intron_utr}))delins(?P<seq>{dna_nt}+))"
"""str: Pattern matching a DNA deletion-insertion with numeric, intronic, or UTR positions.
"""

dna_sub_n: str = dna_sub_c.replace(pos_intron_utr, pos_intron).replace(
    "(?P<dna_sub_c>", "(?P<dna_sub_n>"
)
"""str: Pattern matching a DNA substitution with numeric or intron positions for non-coding variants.
"""

dna_del_n: str = dna_del_c.replace(pos_intron_utr, pos_intron).replace(
    "(?P<dna_del_c>", "(?P<dna_del_n>"
)
"""str: Pattern matching a DNA deletion with numeric or intron positions for non-coding variants.
"""

dna_dup_n: str = dna_dup_c.replace(pos_intron_utr, pos_intron).replace(
    "(?P<dna_dup_c>", "(?P<dna_dup_n>"
)
"""str: Pattern matching a DNA duplication with numeric or intron positions for non-coding variants.
"""

dna_ins_n: str = dna_ins_c.replace(pos_intron_utr, pos_intron).replace(
    "(?P<dna_ins_c>", "(?P<dna_ins_n>"
)
"""str: Pattern matching a DNA insertion with numeric or intron positions for non-coding variants.
"""

dna_delins_n: str = dna_delins_c.replace(pos_intron_utr, pos_intron).replace(
    "(?P<dna_delins_c>", "(?P<dna_delins_n>"
)
"""str: Pattern matching a DNA deletion-insertion with numeric or intron positions for non-coding variants.
"""

dna_sub_gmo: str = dna_sub_c.replace(pos_intron_utr, pos).replace(
    "(?P<dna_sub_c>", "(?P<dna_sub_gmo>"
)
"""str: Pattern matching a DNA substitution with only numeric positions for genomic-style variants.
"""

dna_del_gmo: str = dna_del_c.replace(pos_intron_utr, pos).replace(
    "(?P<dna_del_c>", "(?P<dna_del_gmo>"
)
"""str: Pattern matching a DNA deletion with only numeric positions for genomic-style variants.
"""

dna_dup_gmo: str = dna_dup_c.replace(pos_intron_utr, pos).replace(
    "(?P<dna_dup_c>", "(?P<dna_dup_gmo>"
)
"""str: Pattern matching a DNA duplication with only numeric positions for genomic-style variants.
"""

dna_ins_gmo: str = dna_ins_c.replace(pos_intron_utr, pos).replace(
    "(?P<dna_ins_c>", "(?P<dna_ins_gmo>"
)
"""str: Pattern matching a DNA insertion with only numeric positions for genomic-style variants.
"""

dna_delins_gmo: str = dna_delins_c.replace(pos_intron_utr, pos).replace(
    "(?P<dna_delins_c>", "(?P<dna_delins_gmo>"
)
"""str: Pattern matching a DNA deletion-insertion with only numeric positions for genomic-style variants.
"""

dna_variant_c: str = combine_patterns(
    [dna_sub_c, dna_del_c, dna_dup_c, dna_ins_c, dna_delins_c], None
)
"""str: Pattern matching any of the coding DNA variants.
"""

dna_variant_n: str = combine_patterns(
    [dna_sub_n, dna_del_n, dna_dup_n, dna_ins_n, dna_delins_n], None
)
"""str: Pattern matching any of the non-coding DNA variants.
"""

dna_variant_gmo: str = combine_patterns(
    [dna_sub_gmo, dna_del_gmo, dna_dup_gmo, dna_ins_gmo, dna_delins_gmo], None
)
"""str: Pattern matching any of the genomic-style DNA variants.
"""

dna_single_variant: str = rf"(?P<dna_c>c\.{dna_variant_c})|(?P<dna_n>n\.{dna_variant_n})|(?P<dna_gmo>[gmo]\.{dna_variant_gmo})"
"""str: Pattern matching any complete single DNA variant, including the prefix character.
"""

dna_multi_variant: str = rf"(?P<dna_c_multi>c\.\[{remove_named_groups(dna_variant_c)}(?:;{remove_named_groups(dna_variant_c)}){{1,}}\])|(?P<dna_n_multi>n\.\[{remove_named_groups(dna_variant_n)}(?:;{remove_named_groups(dna_variant_n)}){{1,}}\])|(?P<dna_gmo_multi>[gmo]\.\[{remove_named_groups(dna_variant_gmo)}(?:;{remove_named_groups(dna_variant_gmo)}){{1,}})\]"
"""str: Pattern matching any complete DNA multi-variant, including the prefix character.

Named capture groups have been removed from the variant patterns because of non-uniqueness.
Another applications of single-variant regular expressions is needed to recover the named groups from each individual
variant in the multi-variant.
"""
