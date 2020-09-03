import re

__all__ = [
    "deletion_re",
    "insertion_re",
    "delins_re",
    "substitution_re",
    "multi_variant_re",
    "single_variant_re",
    "any_event_re",
]

nucleotides = "acgubdhkmnrsvwyx"
utr_descriptor = r"(?P<utr>[*-])"
position = r"(?:\d+|\?)"
interval = rf"(?:{position}_{position})"
fragment = rf"(?:\({interval}\))"
intronic_position = r"(?:\d+|\?|\d+(?:[\+-]?(?:\d+|\?))?)"
intronic_interval = rf"(?:{intronic_position}_{intronic_position})"

# Expression with capture groups
edge_cases = r"(?:0|spl|\?)"
deletion = (
    r"(?P<del>"
    rf"(?:(?:\=(?:/|//))?(?P<interval>{interval})del)"
    r"|"
    rf"(?:(?P<fragment>{fragment})del)"
    r"|"
    rf"(?:(?P<position>{position})del(?P<base>[{nucleotides}])?)"
    r")"
)
insertion = (
    r"(?P<ins>"
    rf"(?:(?P<interval>{interval})|(?P<fragment>{fragment}))"
    r"ins"
    r"(?:"
    rf"(?P<intronic>\[(?:{intronic_interval})(?:;{intronic_interval}){{1,}}(?!;)\])"
    r"|"
    rf"(?:(?P<bases>[{nucleotides}]+)|(?P<length>\(\d+\)))"
    r")"
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
    r"(?:"
    r"(?:0|\?|spl)"
    r"|"
    r"(?:"
    rf"(?P<position>{position})"
    r"(?:"
    rf"(?:(?P<mosaic>(?:\=(?:/|//)))?(?P<ref>[{nucleotides}])>(?P<alt>[{nucleotides}]))"
    r"|"
    r"(?P<silent>\=)"
    r")"
    r")"
    r")"
    r")"
)


# Remove capture groups used for use in joining regexes in
# multi-variants since capture groups cannot be defined more than once.
any_event = rf"(?:{r'|'.join([insertion, deletion, delins, substitution])})"
any_event, _ = re.subn(r"P<\w+(_\w+)?>", ":", any_event)

single_variant = rf"r\.{any_event}"
comma_separated = rf"(?:(?:{any_event})(?:,{any_event}){{1,}}(?!,))"
semi_colon_separated = rf"(?:{any_event})(?:;{any_event}){{1,}}(?!;)"
multi_variant = rf"r\.(?:(?:\[{semi_colon_separated}\])|(?:\[{comma_separated}\]))"

# ---- Compiled Regexes
deletion_re = re.compile(rf"(?:r\.)?{utr_descriptor}?{deletion}")
insertion_re = re.compile(rf"(?:r\.)?{utr_descriptor}?{insertion}")
delins_re = re.compile(rf"(?:r\.)?{utr_descriptor}?{delins}")
substitution_re = re.compile(rf"(?:r\.)?{utr_descriptor}?{substitution}")

single_variant_re = re.compile(single_variant)
multi_variant_re = re.compile(multi_variant)
comma_separated_re = re.compile(comma_separated)
semi_colon_separated_re = re.compile(semi_colon_separated)
any_event_re = re.compile(any_event)
