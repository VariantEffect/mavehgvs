from fqfa.constants.iupac.dna import DNA_BASES
from mavehgvs.patterns.util import combine_patterns, remove_named_groups
from mavehgvs.patterns.shared import pos, pos_cnr

dna_nt: str = rf"[{''.join(DNA_BASES)}]"
"""str: Pattern matching any uppercase DNA base.

This does not include IUPAC ambiguity characters.
"""

dna_sub_cn: str = rf"(?P<dna_sub_cn>(?:(?P<position>{pos_cnr})(?P<ref>{dna_nt})>(?P<new>{dna_nt}))|(?P<equal>=))"
"""str: Pattern matching a DNA substitution with numeric or relative-to-transcript positions.
"""

dna_del_cn: str = rf"(?P<dna_del_cn>(?:(?:(?P<start>{pos_cnr})_(?P<end>{pos_cnr}))|(?P<pos>{pos_cnr}))del)"
"""str: Pattern matching a DNA deletion with numeric or relative-to-transcript positions.
"""

dna_dup_cn: str = rf"(?P<dna_dup_cn>(?:(?:(?P<start>{pos_cnr})_(?P<end>{pos_cnr}))|(?P<pos>{pos_cnr}))dup)"
"""str: Pattern matching a DNA duplication with numeric or relative-to-transcript positions.
"""

dna_ins_cn: str = rf"(?P<dna_ins_cn>(?P<start>{pos_cnr})_(?P<end>{pos_cnr})ins(?P<seq>{dna_nt}+))"
"""str: Pattern matching a DNA insertion with numeric or relative-to-transcript positions.
"""

dna_delins_cn: str = rf"(?P<dna_delins_cn>(?:(?:(?P<start>{pos_cnr})_(?P<end>{pos_cnr}))|(?P<pos>{pos_cnr}))delins(?P<seq>{dna_nt}+))"
"""str: Pattern matching a DNA deletion-insertion with numeric or relative-to-transcript positions.
"""

dna_sub_gmo: str = dna_sub_cn.replace(pos_cnr, pos).replace(
    "(?P<dna_sub_cn>", "(?P<dna_sub_gmo>"
)
"""str: Pattern matching a DNA substitution with only numeric positions for genomic-style variants.

This pattern does not match substitutions that are relative to a transcript (e.g. UTR and intronic substitutions).
"""

dna_del_gmo: str = dna_del_cn.replace(pos_cnr, pos).replace(
    "(?P<dna_del_cn>", "(?P<dna_del_gmo>"
)
"""str: Pattern matching a DNA deletion with only numeric positions for genomic-style variants.

This pattern does not match deletions that are relative to a transcript (e.g. UTR and intronic deletions).
"""

dna_dup_gmo: str = dna_dup_cn.replace(pos_cnr, pos).replace(
    "(?P<dna_dup_cn>", "(?P<dna_dup_gmo>"
)
"""str: Pattern matching a DNA duplication with only numeric positions for genomic-style variants.

This pattern does not match duplications that are relative to a transcript (e.g. UTR and intronic duplications).
"""

dna_ins_gmo: str = dna_ins_cn.replace(pos_cnr, pos).replace(
    "(?P<dna_ins_cn>", "(?P<dna_ins_gmo>"
)
"""str: Pattern matching a DNA insertion with only numeric positions for genomic-style variants.

This pattern does not match deletions that are relative to a transcript (e.g. UTR and intronic substitutions).
"""

dna_delins_gmo: str = dna_delins_cn.replace(pos_cnr, pos).replace(
    "(?P<dna_delins_cn>", "(?P<dna_delins_gmo>"
)
"""str: Pattern matching a DNA deletion-insertion with only numeric positions for genomic-style variants.

This pattern does not match deletion-insertions that are relative to a transcript (e.g. UTR and intronic deletion-insertions).
"""

dna_variant_gmo: str = combine_patterns(
    [dna_sub_gmo, dna_del_gmo, dna_dup_gmo, dna_ins_gmo, dna_delins_gmo], None
)
"""str: Pattern matching any of the genomic-style DNA variants.

This pattern does not match duplications that are relative to a transcript (e.g. UTR and intronic duplications).
"""

dna_variant_cn: str = combine_patterns(
    [dna_sub_cn, dna_del_cn, dna_dup_cn, dna_ins_cn, dna_delins_cn], None
)
"""str: Pattern matching any of the transcript-style DNA variants.
"""

dna_single_variant: str = rf"(?P<dna_cn>[cn]\.{dna_variant_cn})|(?P<dna_gmo>[gmo]\.{dna_variant_gmo})"
"""str: Pattern matching any complete single DNA variant, including the prefix character.
"""

dna_multi_variant: str = rf"(?P<dna_cn_multi>[cn]\.\[{remove_named_groups(dna_variant_cn)}(?:;{remove_named_groups(dna_variant_cn)}){{1,}}\])|(?P<dna_gmo_multi>[gmo]\.\[{remove_named_groups(dna_variant_gmo)}(?:;{remove_named_groups(dna_variant_gmo)}){{1,}})\]"
"""str: Pattern matching any complete DNA multi-variant, including the prefix character.

Named capture groups have been removed from the variant patterns because of non-uniqueness.
Another applications of single-variant regular expressions is needed to recover the named groups from each individual
variant in the multi-variant.
"""
