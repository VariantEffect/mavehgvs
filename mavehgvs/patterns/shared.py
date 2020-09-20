pos: str = r"[1-9][0-9]*"
"""str: Pattern matching a positive integer not starting with 0.

This pattern is used for sequence positions, as position 0 does not exist.
"""

pos_cnr: str = rf"[*-]?{pos}(?:[+-]{pos})?"
"""str: Pattern matching a position relative to a transcript.

This pattern is used for sequence positions in a spliced transcript or coding sequence, and can express positions in
the 5' and 3' UTR as well as intronic positions.
"""
