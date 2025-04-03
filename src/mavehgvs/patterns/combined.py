from mavehgvs.patterns.dna import dna_single_variant as dsv, dna_multi_variant as dmv
from mavehgvs.patterns.rna import rna_single_variant as rsv, rna_multi_variant as rmv
from mavehgvs.patterns.protein import (
    pro_single_variant as psv,
    pro_multi_variant as pmv,
)

any_variant = (
    r"(?:(?P<target_id>[a-zA-Z0-9_.-]+):)?"
    + r"(?P<variant>"
    + rf"(?P<single_variant>{r'|'.join([dsv, rsv, psv])})|"
    + rf"(?P<multi_variant>{r'|'.join([dmv, rmv, pmv])})"
    + r")"
)
