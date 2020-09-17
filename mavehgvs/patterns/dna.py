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

pos_tx: str = rf"[*-]?{pos}(?:[+-]{pos})?"
"""str: Pattern matching a position relative to a transcript.

This pattern is used for sequence positions in a spliced transcript or coding sequence, and can express positions in
the 5' and 3' UTR as well as intronic positions.
"""

substitution: str = rf"(?P<substitution>(?:(?P<position>{pos})(?P<ref>{dna_nt})>(?P<alt>{dna_nt}))|(?P<equal>=))"
"""str: Pattern matching a DNA substitution with only numeric positions.

This pattern does not match substitutions that are relative to a transcript (e.g. UTR and intronic substitutions).
"""

substitution_tx: str = rf"(?P<substitution_tx>(?:(?P<position>{pos_tx})(?P<ref>{dna_nt})>(?P<alt>{dna_nt}))|(?P<equal>=))"
"""str: Pattern matching a DNA substitution with numeric or relative-to-transcript positions.
"""

deletion: str = rf"(?P<deletion>(?:(?P<start>{pos})_(?P<end>{pos})del)|(?:(?P<pos>{pos})del))"
"""str: Pattern matching a DNA deletion with only numeric positions.

This pattern does not match deletions that are relative to a transcript (e.g. UTR and intronic deletions).
"""

deletion_tx: str = rf"(?P<deletion_tx>(?:(?P<start>{pos_tx})_(?P<end>{pos_tx})del)|(?:(?P<pos>{pos_tx})del))"
"""str: Pattern matching a DNA deletion with numeric or relative-to-transcript positions.
"""

duplication: str = rf"(?P<duplication>(?:(?P<start>{pos})_(?P<end>{pos})dup)|(?:(?P<pos>{pos})dup))"
"""str: Pattern matching a DNA duplication with only numeric positions.

This pattern does not match duplications that are relative to a transcript (e.g. UTR and intronic duplications).
"""

duplication_tx: str = rf"(?P<duplication_tx>(?:(?P<start>{pos_tx})_(?P<end>{pos_tx})dup)|(?:(?P<pos>{pos_tx})dup))"
"""str: Pattern matching a DNA duplication with numeric or relative-to-transcript positions.
"""

insertion: str = rf"(?P<insertion>(?P<start>{pos})_(?P<end>{pos})ins(?P<bases>{dna_nt}+))"
"""str: Pattern matching a DNA insertion with only numeric positions.

This pattern does not match deletions that are relative to a transcript (e.g. UTR and intronic substitutions).
"""

insertion_tx: str = rf"(?P<insertion_tx>(?P<start>{pos_tx})_(?P<end>{pos_tx})ins(?P<bases>{dna_nt}+))"
"""str: Pattern matching a DNA insertion with numeric or relative-to-transcript positions.
"""

delins: str = rf"(?P<delins>(?:(?P<start>{pos})_(?P<end>{pos})delins)|(?:(?P<pos>{pos})delins)(?P<bases>{dna_nt}+))"
"""str: Pattern matching a DNA deletion-insertion with only numeric positions.

This pattern does not match deletion-insertions that are relative to a transcript (e.g. UTR and intronic deletion-insertions).
"""

delins_tx: str = rf"(?P<delins>(?:(?P<start>{pos_tx})_(?P<end>{pos_tx})delins)|(?:(?P<pos>{pos_tx})delins)(?P<bases>{dna_nt}+))"
"""str: Pattern matching a DNA deletion-insertion with numeric or relative-to-transcript positions.
"""


# Remove capture groups used for use in joining regexes in
# multi-variants since capture groups cannot be defined more than once.
any_event_re = combine_patterns(
    [substitution, deletion, duplication, insertion], "any_event"
)

any_event_tx_re = combine_patterns(
    [substitution_tx, deletion_tx, duplication_tx, insertion_tx], "any_event_tx"
)

single_variant_re = re.compile(
    rf"(?P<dna_tx>[cn]\.{any_event_tx_re.pattern})|(?P<dna>[gmo]\.{any_event_re.pattern})"
)


multi_variant = rf"(?P<dna_tx_multi>[cn]\.\[{remove_named_groups(any_event_tx_re.pattern)}(?:;{remove_named_groups(any_event_tx_re.pattern)}){{1,}}\])|(?P<dna_multi>[gmo]\.\[{remove_named_groups(any_event_re.pattern)}(?:;{remove_named_groups(any_event_re.pattern)}){{1,}})\]"
multi_variant_re = re.compile(multi_variant)
# Another pass of regexes will be needed to recover the various capture groups after initial validation

# ---- Compiled Regexes
# deletion_re = combine_patterns([rf"(?:[cn]\.{deletion_tx})", rf"(?:[gmo]\.{deletion})"])
# duplication_re = combine_patterns([rf"(?:[cn]\.{duplication_tx})", rf"(?:[gmo]\.{duplication})"])
# insertion_re = combine_patterns([rf"(?:[cn]\.{insertion_tx})", rf"(?:[gmo]\.{insertion})"])
# delins_re = combine_patterns([rf"(?:[cn]\.{delins_tx})", rf"(?:[gmo]\.{delins})"])
# substitution_re = combine_patterns([rf"(?:[cn]\.{substitution_tx})", rf"(?:[gmo]\.{substitution})"])

deletion_re = re.compile(deletion_tx)
duplication_re = re.compile(duplication_tx)
insertion_re = re.compile(insertion_tx)
delins_re = re.compile(delins_tx)
substitution_re = re.compile(substitution_tx)
