import re

__all__ = [
    'deletion_re', 'insertion_re', 'delins_re',
    'substitution_re', 'multi_variant_re', 'single_variant_re',
    'any_event_re', 'frame_shift_re', 'predicted_variant_re'
]

#: Conversions between single- and three-letter amino acid codes
AA_CODES = {
        'Ala': 'A', 'A': 'Ala',
        'Arg': 'R', 'R': 'Arg',
        'Asn': 'N', 'N': 'Asn',
        'Asp': 'D', 'D': 'Asp',
        'Cys': 'C', 'C': 'Cys',
        'Glu': 'E', 'E': 'Glu',
        'Gln': 'Q', 'Q': 'Gln',
        'Gly': 'G', 'G': 'Gly',
        'His': 'H', 'H': 'His',
        'Ile': 'I', 'I': 'Ile',
        'Leu': 'L', 'L': 'Leu',
        'Lys': 'K', 'K': 'Lys',
        'Met': 'M', 'M': 'Met',
        'Phe': 'F', 'F': 'Phe',
        'Pro': 'P', 'P': 'Pro',
        'Ser': 'S', 'S': 'Ser',
        'Thr': 'T', 'T': 'Thr',
        'Trp': 'W', 'W': 'Trp',
        'Tyr': 'Y', 'Y': 'Tyr',
        'Val': 'V', 'V': 'Val',
        'Ter': '*', '*': 'Ter',
        '???': '?', '?': '???',
}

amino_acids = '(?:{})'.format(
    '|'.join(set(AA_CODES.keys())).replace('?', '\?').replace('*', '\*')
)

position = r"(?:(?:{0}\d+)|\?)".format(amino_acids)
interval = r"(?:{0}_{0})".format(position)
amino_acid_choice = r"(?:(?:{0}){{1}}(?:\^(?:{0}))+(?!\^))".format(amino_acids)


# Expression with capture groups
deletion = (
    r"(?P<del>"
        r"(?:"
            r"(?P<interval>{0})"
            r"|"
            r"(?:(?P<position>{1})(?P<mosaic>\=/)?)"
        r")"
        r"del"
    r")".format(interval, position)
)
insertion = (
    r"(?P<ins>"
        r"(?P<interval>{0})"
        r"ins"
        r"(?:"
            r"(?P<inserted>{1}+|{2})"
            r"|"
            r"(?P<length>\d+)"
            r"|"
            r"(?P<unknown>(?:\(\d+\))|X+)"
        r")"
    r")".format(interval, amino_acids, amino_acid_choice)
)
delins = (
    r"(?P<delins>"
        r"(?:"
            r"(?P<interval>{0})"
            r"|"
            r"(?P<position>{1})"
        r")"
        r"delins"
        r"(?:"
            r"(?P<inserted>{2}+|{3})"
            r"|"
            r"(?P<length>\d+)"
            r"|"
            r"(?P<unknown>(?:\(\d+\))|X+)"
        r")"
    r")".format(interval, position, amino_acids, amino_acid_choice)
)
substitution = (
    r"(?P<sub>"
        r"(?:(?P<no_protein>0)|(?P<not_predicted>\?)|(?P<equal>=))"
        r"|"
        r"(?:"
            r"(?P<pre>{0})(?P<position>\d+)"
            r"(?:"
                r"(?P<post>(?:(?P<mosaic>\=/)?(?:{0}))|(?P<choice>{1})|(?:\*))"
                r"|"
                r"(?P<silent>\=)"
                r"|"
                r"(?P<unknown>\?)"
            r")"
        r")"
    r")".format(amino_acids, amino_acid_choice)
)
frame_shift = (
    r"(?P<fs>"
        r"(?P<left_aa>{0})(?P<position>\d+)(?P<right_aa>{0})?fs"
        r"(?P<shift>"
            r"(?:"
                r"(?:{0}\d+)"
                r"|"
                r"(?:\*\?)"
                r"|"
                r"(?:\*\d+)"
            r")"
        r")?"
    r")"
).format(amino_acids)


# Remove capture groups used for use in joining regexes in
# multi-variants since capture groups cannot be defined more than once.
any_event = r"(?:{0})".format(
    r"|".join([insertion, deletion, delins, substitution, frame_shift])
)
any_event, _ = re.subn(r"P<\w+(_\w+)?>", ':', any_event)
predicted_event = r"\({0}\)".format(any_event)
predicted_variant = r"p.\({0}\)".format(any_event)
single_variant = r"(?:p\.{0})|(?:{1})".format(any_event, predicted_variant)

multi_any = r"({}|{})".format(predicted_event, any_event)
multi_variant = r"p\.\[(?:{0})(?:;{0}){{1,}}(?!;)\]".format(multi_any)


# ---- Compiled Regexes
deletion_re = re.compile(
    r"(?:p\.)?(?P<predicted>\()?{0}(?(predicted)\)|)".format(deletion))
insertion_re = re.compile(
    r"(?:p\.)?(?P<predicted>\()?{0}(?(predicted)\)|)".format(insertion))
delins_re = re.compile(
    r"(?:p\.)?(?P<predicted>\()?{0}(?(predicted)\)|)".format(delins))
substitution_re = re.compile(
    r"(?:p\.)?(?P<predicted>\()?{0}(?(predicted)\)|)".format(substitution))
frame_shift_re = re.compile(
    r"(?:p\.)?(?P<predicted>\()?{0}(?(predicted)\)|)".format(frame_shift))

single_variant_re = re.compile(single_variant)
multi_variant_re = re.compile(multi_variant)
predicted_variant_re = re.compile(predicted_variant)
any_event_re = re.compile(any_event)


def split_amino_acids(aa_str):
    return re.findall('[A-Z\?][^A-Z\?]*', aa_str)
