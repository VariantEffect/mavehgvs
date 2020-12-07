from mavehgvs.patterns.dna import dna_single_variant, dna_multi_variant
from mavehgvs.patterns.rna import rna_single_variant, rna_multi_variant
from mavehgvs.patterns.protein import pro_single_variant, pro_multi_variant

any_variant = rf"(?:(?P<target_id>[a-zA-Z0-9_.-]+):)?(?P<single_variant>{r'|'.join([dna_single_variant, rna_single_variant, pro_single_variant])})|(?P<multi_variant>{r'|'.join([dna_multi_variant, rna_multi_variant, pro_multi_variant])})"
