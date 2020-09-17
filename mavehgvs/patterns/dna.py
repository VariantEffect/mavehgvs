import re
from fqfa.constants.iupac.dna import DNA_BASES
from mavehgvs.patterns.util import combine_patterns, remove_named_groups


dna_nt: str = rf"[{''.join(DNA_BASES)}]"
"""str: Pattern matching any uppercase base.

This does not include IUPAC ambiguity characters.
"""

pos: str = r"[1-9][0-9]*"
"""str: Pattern matching a positive integer not starting with 0.

This pattern is used for sequence positions, as position 0 does not exist.
"""

pos_cn: str = rf"[*-]?{pos}(?:[+-]{pos})?"
"""str: Pattern matching a position relative to a transcript.

This pattern is used for sequence positions in a spliced transcript or coding sequence, and can express positions in
the 5' and 3' UTR as well as intronic positions.
"""

substitution_cn: str = rf"(?P<substitution_cn>(?:(?P<position>{pos_cn})(?P<ref>{dna_nt})>(?P<alt>{dna_nt}))|(?P<equal>=))"
"""str: Pattern matching a DNA substitution with numeric or relative-to-transcript positions.
"""

deletion_cn: str = rf"(?P<deletion_cn>(?:(?P<start>{pos_cn})_(?P<end>{pos_cn})del)|(?:(?P<pos>{pos_cn})del))"
"""str: Pattern matching a DNA deletion with numeric or relative-to-transcript positions.
"""

duplication_cn: str = rf"(?P<duplication_cn>(?:(?P<start>{pos_cn})_(?P<end>{pos_cn})dup)|(?:(?P<pos>{pos_cn})dup))"
"""str: Pattern matching a DNA duplication with numeric or relative-to-transcript positions.
"""

insertion_cn: str = rf"(?P<insertion_cn>(?P<start>{pos_cn})_(?P<end>{pos_cn})ins(?P<bases>{dna_nt}+))"
"""str: Pattern matching a DNA insertion with numeric or relative-to-transcript positions.
"""

delins_cn: str = rf"(?P<delins>(?:(?P<start>{pos_cn})_(?P<end>{pos_cn})delins)|(?:(?P<pos>{pos_cn})delins)(?P<bases>{dna_nt}+))"
"""str: Pattern matching a DNA deletion-insertion with numeric or relative-to-transcript positions.
"""

substitution_gmo: str = substitution_cn.replace(pos_cn, pos).replace("(?P<substitution_cn>", "(?P<substitution_gmo>")
"""str: Pattern matching a DNA substitution with only numeric positions for genomic-style variants.

This pattern does not match substitutions that are relative to a transcript (e.g. UTR and intronic substitutions).
"""

deletion_gmo: str = deletion_cn.replace(pos_cn, pos).replace("(?P<deletion_cn>", "(?P<deletion_gmo>")
"""str: Pattern matching a DNA deletion with only numeric positions for genomic-style variants.

This pattern does not match deletions that are relative to a transcript (e.g. UTR and intronic deletions).
"""

duplication_gmo: str = duplication_cn.replace(pos_cn, pos).replace("(?P<duplication_cn>", "(?P<duplication_gmo>")
"""str: Pattern matching a DNA duplication with only numeric positions for genomic-style variants.

This pattern does not match duplications that are relative to a transcript (e.g. UTR and intronic duplications).
"""

insertion_gmo: str = insertion_cn.replace(pos_cn, pos).replace("(?P<insertion_cn>", "(?P<insertion_gmo>")
"""str: Pattern matching a DNA insertion with only numeric positions for genomic-style variants.

This pattern does not match deletions that are relative to a transcript (e.g. UTR and intronic substitutions).
"""

delins_gmo: str = delins_cn.replace(pos_cn, pos).replace("(?P<delins_cn>", "(?P<delins_gmo>")
"""str: Pattern matching a DNA deletion-insertion with only numeric positions for genomic-style variants.

This pattern does not match deletion-insertions that are relative to a transcript (e.g. UTR and intronic deletion-insertions).
"""

variant_gmo = combine_patterns([substitution_gmo, deletion_gmo, duplication_gmo, insertion_gmo, delins_gmo], None)

variant_cn = combine_patterns([substitution_cn, deletion_cn, duplication_cn, insertion_cn, delins_cn], None)

dna_single_variant = rf"(?P<dna_cn>[cn]\.{variant_cn})|(?P<dna_gmo>[gmo]\.{variant_gmo})"

dna_multi_variant = rf"(?P<dna_cn_multi>[cn]\.\[{remove_named_groups(variant_cn)}(?:;{remove_named_groups(variant_cn)}){{1,}}\])|(?P<dna_gmo_multi>[gmo]\.\[{remove_named_groups(variant_gmo)}(?:;{remove_named_groups(variant_gmo)}){{1,}})\]"
# Another pass of regexes will be needed to recover the various capture groups after initial validation

# ---- Compiled Regexes
# deletion_re = combine_patterns([rf"(?:[cn]\.{deletion_cn})", rf"(?:[gmo]\.{deletion})"])
# duplication_re = combine_patterns([rf"(?:[cn]\.{duplication_cn})", rf"(?:[gmo]\.{duplication})"])
# insertion_re = combine_patterns([rf"(?:[cn]\.{insertion_cn})", rf"(?:[gmo]\.{insertion})"])
# delins_re = combine_patterns([rf"(?:[cn]\.{delins_cn})", rf"(?:[gmo]\.{delins})"])
# substitution_re = combine_patterns([rf"(?:[cn]\.{substitution_cn})", rf"(?:[gmo]\.{substitution})"])

deletion_re = re.compile(deletion_cn)
duplication_re = re.compile(duplication_cn)
insertion_re = re.compile(insertion_cn)
delins_re = re.compile(delins_cn)
substitution_re = re.compile(substitution_cn)
