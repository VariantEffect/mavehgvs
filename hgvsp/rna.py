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
interval = r"(?:{0}_{0})".format(position)
fragment = r"(?:\({0}\))".format(interval)
intronic_position = r"(?:\d+|\?|\d+(?:[\+-]?(?:\d+|\?))?)"
intronic_interval = r"(?:{0}_{0})".format(intronic_position)

# Expression with capture groups
edge_cases = r"(?:0|spl|\?)"
deletion = (
    r"(?P<del>"
    r"(?:(?:\=(?:/|//))?(?P<interval>{0})del)"
    r"|"
    r"(?:(?P<fragment>{1})del)"
    r"|"
    r"(?:(?P<position>{2})del(?P<base>[{3}])?)"
    r")".format(interval, fragment, position, nucleotides)
)
insertion = (
    r"(?P<ins>"
    r"(?:(?P<interval>{0})|(?P<fragment>{1}))"
    r"ins"
    r"(?:"
    r"(?P<intronic>\[(?:{2})(?:;{2}){{1,}}(?!;)\])"
    r"|"
    r"(?:(?P<bases>[{3}]+)|(?P<length>\(\d+\)))"
    r")"
    r")".format(interval, fragment, intronic_interval, nucleotides)
)
delins = (
    r"(?P<delins>"
    r"(?:"
    r"(?:(?P<interval>{0})delins)"
    r"|"
    r"(?:(?P<position>{1})delins)"
    r")"
    r"(?:(?P<bases>[{2}]+)|(?P<length>\(\d+\)))"
    r")".format(interval, position, nucleotides)
)
substitution = (
    r"(?P<sub>"
    r"(?:"
    r"(?:0|\?|spl)"
    r"|"
    r"(?:"
    r"(?P<position>{0})"
    r"(?:"
    r"(?:(?P<mosaic>(?:\=(?:/|//)))?(?P<ref>[{1}])>(?P<alt>[{1}]))"
    r"|"
    r"(?P<silent>\=)"
    r")"
    r")"
    r")"
    r")".format(position, nucleotides)
)


# Remove capture groups used for use in joining regexes in
# multi-variants since capture groups cannot be defined more than once.
any_event = r"(?:{})".format(r"|".join([insertion, deletion, delins, substitution]))
any_event, _ = re.subn(r"P<\w+(_\w+)?>", ":", any_event)

single_variant = r"r\.{0}".format(any_event)
comma_separated = r"(?:(?:{0})(?:,{0}){{1,}}(?!,))".format(any_event)
semi_colon_separated = r"(?:{0})(?:;{0}){{1,}}(?!;)".format(any_event)
multi_variant = r"r\.(?:(?:\[{0}\])|(?:\[{1}\]))".format(
    semi_colon_separated, comma_separated
)

# ---- Compiled Regexes
deletion_re = re.compile(r"(?:r\.)?{0}?{1}".format(utr_descriptor, deletion))
insertion_re = re.compile(r"(?:r\.)?{0}?{1}".format(utr_descriptor, insertion))
delins_re = re.compile(r"(?:r\.)?{0}?{1}".format(utr_descriptor, delins))
substitution_re = re.compile(r"(?:r\.)?{0}?{1}".format(utr_descriptor, substitution))

single_variant_re = re.compile(single_variant)
multi_variant_re = re.compile(multi_variant)
comma_separated_re = re.compile(comma_separated)
semi_colon_separated_re = re.compile(semi_colon_separated)
any_event_re = re.compile(any_event)
