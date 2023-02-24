"""Utility functions for working with mavehgvs regex pattern strings.
"""

import re
from typing import Sequence, Optional


def combine_patterns(patterns: Sequence[str], groupname: Optional[str] = None) -> str:
    """Combine multiple pattern strings into a single pattern string.

    Because multiple identical group names are not allowed in a pattern, the resulting object renames all named match
    groups such they are prefixed with the first match group name in the pattern. For example,
    ``(?P<substitution>(?P<position>[1-9][0-9]*)...`` becomes
    ``(?P<substitution>(?P<substitution_position>[1-9][0-9]*)...``.

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
    str
        Pattern string that matches any of the input patterns. Match groups are renamed as described above to attempt
        to ensure uniqueness across the combined pattern.

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

    return combined


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

    return re.sub(r"\(\?P<\w+>", new_parens, pattern)
