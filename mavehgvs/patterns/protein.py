from fqfa.constants import AA_CODES
from mavehgvs.patterns.util import combine_patterns, remove_named_groups
from mavehgvs.patterns.position import pos

amino_acid: str = rf"(?:{'|'.join(AA_CODES.values())})"
"""str: Pattern matching any amino acid or Ter.

This does not include ambiguous amino acids such as Glx and Xaa.
"""

aa_pos: str = rf"(?:{amino_acid}{pos})"
"""str: Pattern matching an amino acid code followed by a position.
"""

pro_sub: str = rf"(?P<pro_sub>(?:(?P<position>{aa_pos})(?P<new>{amino_acid}|=))|(?P<equal>=)|(?P<equal_sy>\(=\)))"
"""str: Pattern matching a protein substitution.
"""

pro_del: str = rf"(?P<pro_del>(?:(?P<start>{aa_pos})_(?P<end>{aa_pos})del)|(?:(?P<pos>{aa_pos})del))"
"""str: Pattern matching a protein deletion.
"""

pro_dup: str = rf"(?P<pro_dup>(?:(?P<start>{aa_pos})_(?P<end>{aa_pos})dup)|(?:(?P<pos>{aa_pos})dup))"
"""str: Pattern matching a protein duplication.
"""

pro_ins: str = rf"(?P<pro_ins>(?P<start>{aa_pos})_(?P<end>{aa_pos})ins(?P<seq>{amino_acid}+))"
"""str: Pattern matching a protein insertion.
"""

pro_delins: str = rf"(?P<pro_delins>(?:(?:(?P<start>{aa_pos})_(?P<end>{aa_pos}))|(?P<pos>{aa_pos}))delins(?P<seq>{amino_acid}+))"
"""str: Pattern matching a protein deletion-insertion.
"""

pro_variant: str = combine_patterns(
    [pro_sub, pro_del, pro_dup, pro_ins, pro_delins], None
)
"""str: Pattern matching any single protein variant event.
"""

pro_single_variant: str = rf"(?P<pro>p\.{pro_variant})"
"""str: Pattern matching any complete protein variant, including the prefix character.
"""

pro_multi_variant: str = rf"(?P<pro_multi>p\.\[{remove_named_groups(pro_variant)}(?:;{remove_named_groups(pro_variant)}){{1,}}\])"
"""str: Pattern matching any complete protein multi-variant, including the prefix character.

Named capture groups have been removed from the variant patterns because of non-uniqueness.
Another applications of single-variant regular expressions is needed to recover the named groups from each individual
variant in the multi-variant.
"""
