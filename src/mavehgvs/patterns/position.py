pos: str = r"[1-9][0-9]*"
"""str: Pattern matching a positive integer not starting with 0.

This pattern is used for sequence positions, as position 0 does not exist.
"""

pos_intron: str = rf"{pos}(?:[+-]{pos})?"
"""str: Pattern matching a position with optional intron component.

This pattern is used for sequence positions in an RNA or noncoding sequence.
"""

pos_intron_utr: str = rf"[*-]?{pos}(?:[+-]{pos})?"
"""str: Pattern matching a position with optional intron and UTR components.

This pattern is used for sequence positions in a coding sequence.
"""
