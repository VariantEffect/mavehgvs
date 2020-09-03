import re

from fqfa.constants.iupac.protein import AA_CODES_ALL

__all__ = [
    "deletion_re",
    "insertion_re",
    "delins_re",
    "substitution_re",
    "multi_variant_re",
    "single_variant_re",
    "any_event_re",
    "frame_shift_re",
    "predicted_variant_re",
]

amino_acids = rf"(?:{'|'.join(AA_CODES_ALL.values())})"

position = rf"(?:(?:{amino_acids}\d+)|\?)"
interval = rf"(?:{position}_{position})"
amino_acid_choice = rf"(?:(?:{amino_acids}){{1}}(?:\^(?:{amino_acids}))+(?!\^))"


# Expression with capture groups
deletion = (
    r"(?P<del>"
    r"(?:"
    rf"(?P<interval>{interval})"
    r"|"
    rf"(?:(?P<position>{position})(?P<mosaic>\=/)?)"
    r")"
    r"del"
    r")"
)
insertion = (
    r"(?P<ins>"
    rf"(?P<interval>{interval})"
    r"ins"
    r"(?:"
    rf"(?P<inserted>{amino_acids}+|{amino_acid_choice})"
    r"|"
    r"(?P<length>\d+)"
    r"|"
    r"(?P<unknown>(?:\(\d+\))|X+)"
    r")"
    r")"
)
delins = (
    r"(?P<delins>"
    r"(?:"
    rf"(?P<interval>{interval})"
    r"|"
    rf"(?P<position>{position})"
    r")"
    r"delins"
    r"(?:"
    rf"(?P<inserted>{amino_acids}+|{amino_acid_choice})"
    r"|"
    r"(?P<length>\d+)"
    r"|"
    r"(?P<unknown>(?:\(\d+\))|X+)"
    r")"
    r")"
)
substitution = (
    r"(?P<sub>"
    r"(?:(?P<no_protein>0)|(?P<not_predicted>\?)|(?P<equal>=))"
    r"|"
    r"(?:"
    rf"(?P<pre>{amino_acids})(?P<position>\d+)"
    r"(?:"
    rf"(?P<post>(?:(?P<mosaic>\=/)?(?:{amino_acids}))|(?P<choice>{amino_acid_choice})|(?:\*))"
    r"|"
    r"(?P<silent>\=)"
    r")"
    r")"
    r")"
)
frame_shift = (
    r"(?P<fs>"
    rf"(?P<left_aa>{amino_acids})(?P<position>\d+)(?P<right_aa>{amino_acids})?fs"
    r"(?P<shift>"
    r"(?:"
    rf"(?:{amino_acids}\d+)"
    r"|"
    r"(?:\*\?)"
    r"|"
    r"(?:\*\d+)"
    r")"
    r")?"
    r")"
)


# Remove capture groups used for use in joining regexes in
# multi-variants since capture groups cannot be defined more than once.
any_event = rf"(?:{'|'.join([insertion, deletion, delins, substitution, frame_shift])})"
any_event, _ = re.subn(r"P<\w+(_\w+)?>", ":", any_event)
predicted_event = rf"\({any_event}\)"
predicted_variant = rf"p.\({any_event}\)"
single_variant = rf"(?:p\.{any_event})|(?:{predicted_variant})"

multi_any = rf"({predicted_event}|{any_event})"
multi_variant = rf"p\.\[(?:{multi_any})(?:;{multi_any}){{1,}}(?!;)\]"


# ---- Compiled Regexes
deletion_re = re.compile(
    rf"(?:p\.)?(?P<predicted>\()?{deletion}(?(predicted)\)|)"
)
insertion_re = re.compile(
    rf"(?:p\.)?(?P<predicted>\()?{insertion}(?(predicted)\)|)"
)
delins_re = re.compile(rf"(?:p\.)?(?P<predicted>\()?{delins}(?(predicted)\)|)")
substitution_re = re.compile(
    rf"(?:p\.)?(?P<predicted>\()?{substitution}(?(predicted)\)|)"
)
frame_shift_re = re.compile(
    rf"(?:p\.)?(?P<predicted>\()?{frame_shift}(?(predicted)\)|)"
)

single_variant_re = re.compile(single_variant)
multi_variant_re = re.compile(multi_variant)
predicted_variant_re = re.compile(predicted_variant)
any_event_re = re.compile(any_event)


def split_amino_acids(aa_str):
    return re.findall("[A-Z\?][^A-Z\?]*", aa_str)
