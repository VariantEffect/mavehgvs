import re
from typing import Sequence, Optional, Pattern
from fqfa.constants.iupac.dna import DNA_BASES

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


def combine_patterns(
    patterns: Sequence[str], groupname: Optional[str] = None
) -> Pattern:
    """Combine multiple pattern strings and generate a compiled regular expression object.

    Because multiple identical group names are not allowed in a pattern, the resulting object renames all named match
    groups such they are prefixed with the first match group name in the pattern. For example,
    :code:`(?P<substitution>(?P<position>[1-9][0-9]*)...` becomes
    :code:`(?P<substitution>(?P<substitution_position>[1-9][0-9]*)...`.

    The function assumes that all input patterns are enclosed in parentheses.

    Parameters
    ----------
    patterns : Sequence[str]
        Sequence of pattern strings to combine.

    groupname : Optional[str]
        Name for the capture group surrounding the resulting pattern. If this is None, a non-capturing group will be
        used instead.

    Returns
    -------
    re.Pattern
        Compiled regular expression that matches any of the input patterns. Match groups are renamed as described above
        to attempt to ensure uniqueness across the combined pattern.

    """
    tag_re = re.compile(r"\(\?P<(\w+)>")
    stripped_patterns = list()
    for p in patterns:
        tags = list(tag_re.finditer(p))
        prefix = f"{tags[0].group(1)}_"
        new_p = p
        for t in tags[:0:-1]:
            start, end = t.span(1)
            new_p = "".join((new_p[:start], prefix, new_p[start:]))
        stripped_patterns.append(new_p)
    if groupname is None:
        combined = rf"(?:{r'|'.join(stripped_patterns)})"
    else:
        combined = rf"(?P<{groupname}>{r'|'.join(stripped_patterns)})"

    return re.compile(combined)


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


def remove_named_groups(pattern: str, noncapturing: bool = True) -> str:
    """Function that replaces named match groups in a regular expression pattern.

    Named groups are replaced with either regular parentheses or non-capturing parentheses.

    Parameters
    ----------
    pattern : str
        The pattern string to strip match groups from.

    noncapturing : bool
        If True, the named grouping parentheses are replaced by non-capturing parentheses.
        If False, regular parentheses are used.

    Returns
    -------
    str
        The pattern string without named match groups.

    """
    if noncapturing:
        new_parens = "(?:"
    else:
        new_parens = "("

    return re.sub(r'\(\?P<\w+>', new_parens, pattern)


multi_variant = rf"(?P<dna_tx_multi>[cn]\.\[{remove_named_groups(any_event_tx_re.pattern)}(?:;{remove_named_groups(any_event_tx_re.pattern)}){{1,}}\])|(?P<dna_multi>[gmo]\.\[{remove_named_groups(any_event_re.pattern)}(?:;{remove_named_groups(any_event_re.pattern)}){{1,}})\]"
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
