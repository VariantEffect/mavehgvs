import re

__all__ = [
    'deletion_re', 'insertion_re', 'delins_re',
    'substitution_re', 'multi_variant_re', 'single_variant_re',
    'any_event_re',
]

nucleotides = 'ATCGXNH'
utr_descriptor = r"(?P<utr>[*-])"
position = r"(?:\d+|\?|(?:[*-]?\d+(?:[\+-]?(?:\d+|\?))?))"
interval = r"(?:{0}_{0})".format(position)
fragment = r"(?:\({0}\))".format(interval)
breakpoint_ = r"(?:{0}_{0})".format(fragment)


# Expression with capture groups
deletion = (
    r"(?P<del>"
        r"(?:(?P<interval>{0})(?:(?:\=(?:/|//))|(?:del\=(?:/|//)))?del)"
        r"|"
        r"(?:(?P<breakpoint>{1})del)"
        r"|"
        r"(?:(?P<position>{2})del(?P<base>[{3}])?)"
    r")".format(interval, breakpoint_, position, nucleotides)
)
insertion = (
    r"(?P<ins>"
        r"(?:"
            r"(?:(?P<interval>{0})ins)"
            r"|"
            r"(?:(?P<fragment>{1})ins)"
        r")"
        r"(?:(?P<bases>[{2}]+)|(?P<length>\(\d+\)))"
    r")".format(interval, fragment, nucleotides)
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
        r"(?P<position>{0})"
        r"(?:"
            r"(?:(?P<mosaic>(?:\=/)|(?:\=//))?(?P<ref>[{1}])>(?P<alt>[{1}]))"
            r"|"
            r"(?P<silent>\=)"
        r")"
    r")".format(position, nucleotides)
)


# Remove capture groups used for use in joining regexes in
# multi-variants since capture groups cannot be defined more than once.
any_event = r"{0}?(?:{1})".format(
    utr_descriptor,
    r"|".join([insertion, deletion, delins, substitution]))
any_event, _ = re.subn(r"P<\w+(_\w+)?>", ':', any_event)

single_variant = r"[cngm]\.{0}".format(any_event)
multi_variant = r"[cngm]\.\[(?:{0})(?:;{0}){{1,}}(?!;)\]".format(any_event)


# ---- Compiled Regexes
deletion_re = re.compile(
    r"(?:[cngm]\.)?{0}?{1}".format(utr_descriptor, deletion))
insertion_re = re.compile(
    r"(?:[cngm]\.)?{0}?{1}".format(utr_descriptor, insertion))
delins_re = re.compile(
    r"(?:[cngm]\.)?{0}?{1}".format(utr_descriptor, delins))
substitution_re = re.compile(
    r"(?:[cngm]\.)?{0}?{1}".format(utr_descriptor, substitution))

single_variant_re = re.compile(single_variant)
multi_variant_re = re.compile(multi_variant)
any_event_re = re.compile(any_event)
