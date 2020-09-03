import re

from fqfa.constants.iupac.dna import DNA_CHARACTERS

__all__ = [
    "deletion_re",
    "insertion_re",
    "delins_re",
    "substitution_re",
    "multi_variant_re",
    "single_variant_re",
    "any_event_re",
]

nucleotides = ''.join(DNA_CHARACTERS)
utr_descriptor = r"(?P<utr>[*-])"
position = r"(?:\d+|\?|(?:[*-]?\d+(?:[\+-]?(?:\d+|\?))?))"
interval = rf"(?:{position}_{position})"
fragment = rf"(?:\({interval}\))"
breakpoint_ = rf"(?:{fragment}_{fragment})"


# Expression with capture groups
deletion = (
    r"(?P<del>"
    rf"(?:(?P<interval>{interval})(?:(?:\=(?:/|//))|(?:del\=(?:/|//)))?del)"
    r"|"
    rf"(?:(?P<breakpoint>{breakpoint_})del)"
    r"|"
    rf"(?:(?P<position>{position})del(?P<base>[{nucleotides}])?)"
    r")"
)
insertion = (
    r"(?P<ins>"
    r"(?:"
    rf"(?:(?P<interval>{interval})ins)"
    r"|"
    rf"(?:(?P<fragment>{fragment})ins)"
    r")"
    rf"(?:(?P<bases>[{nucleotides}]+)|(?P<length>\(\d+\)))"
    r")"
)
delins = (
    r"(?P<delins>"
    r"(?:"
    rf"(?:(?P<interval>{interval})delins)"
    r"|"
    rf"(?:(?P<position>{position})delins)"
    r")"
    rf"(?:(?P<bases>[{nucleotides}]+)|(?P<length>\(\d+\)))"
    r")"
)
substitution = (
    r"(?P<sub>"
    rf"(?P<position>{position})"
    r"(?:"
    rf"(?:(?P<mosaic>(?:\=/)|(?:\=//))?(?P<ref>[{nucleotides}])>(?P<alt>[{nucleotides}]))"
    r"|"
    r"(?P<silent>\=)"
    r")"
    r")"
)


# Remove capture groups used for use in joining regexes in
# multi-variants since capture groups cannot be defined more than once.
any_event = rf"{utr_descriptor}?(?:{r'|'.join([insertion, deletion, delins, substitution])})"
any_event, _ = re.subn(r"P<\w+(_\w+)?>", ":", any_event)

single_variant = rf"[cngm]\.{any_event}"
multi_variant = rf"[cngm]\.\[(?:{any_event})(?:;{any_event}){{1,}}(?!;)\]"


# ---- Compiled Regexes
deletion_re = re.compile(rf"(?:[cngm]\.)?{utr_descriptor}?{deletion}")
insertion_re = re.compile(rf"(?:[cngm]\.)?{utr_descriptor}?{insertion}")
delins_re = re.compile(rf"(?:[cngm]\.)?{utr_descriptor}?{delins}")
substitution_re = re.compile(rf"(?:[cngm]\.)?{utr_descriptor}?{substitution}")

single_variant_re = re.compile(single_variant)
multi_variant_re = re.compile(multi_variant)
any_event_re = re.compile(any_event)
