import re
from functools import total_ordering
from mavehgvs.patterns.shared import pos

pos_extended: str = rf"(?P<position>[*-]?{pos})(?P<position_intron>[+-]{pos})?"
"""str: Pattern matching a position relative to a transcript with match groups.

This pattern is used for sequence positions in a spliced transcript or coding sequence that are found in a UTR or
intron. It does not match integer-only positions.
"""


@total_ordering
class VariantPosition:
    """Class for storing a variant position.

     The class includes special fields for variants using the extended position syntax.

     Attributes
     ----------
     position : Optional[int]
         The position as an integer.
         Negative positions are only expected for 5' UTR positions.
     intronic_position : Optional[int]
         The number of bases into the intron for intronic positions.
         None for non-intronic positions.

         Nucleotides in the 5' half of the intron have positive ``intronic_position`` and their position is that of
         the last base of the 5' exon.
         Nucleotides in the 3' half of the intron have negative ``intronic_position`` and their position is that of
         the first base of the 3' exon.
     utr : Optional[bool]
         True if the position is in the UTR. None for all other positions.

     """

    _fullmatch = re.compile(pos_extended, flags=re.ASCII).fullmatch
    """Callable[[str, int, int], Optional[Match[str]]]: fullmatch callable for parsing positions
    
    Returns a match object if the full string matches one of the position groups in :py:data:`pos_extended`.
    """

    def __init__(self, pos_str: str) -> None:
        """Parse a position string into a VariantPosition object.

        Parameters
        ----------
        pos_str : str
            The string to convert to a VariantPosition object.

        """
        try:
            gdict = VariantPosition._fullmatch(pos_str).groupdict()
        except AttributeError:
            raise ValueError(f"invalid variant position string '{pos_str}'")

        self.position = None
        self.intronic_position = None
        self.utr = None

        if gdict["position"] is not None:
            if gdict["position"].startswith("*"):
                self.utr = True
                self.position = int(gdict["position"][1:])
            else:
                if gdict["position"].startswith("-"):
                    self.utr = True
                self.position = int(gdict["position"])
        else:
            raise ValueError(f"invalid variant position string '{pos_str}'")

        if gdict["position_intron"] is not None:
            self.intronic_position = int(gdict["position_intron"])

    def is_utr(self) -> bool:
        """Return whether this is a UTR position.

        Returns
        -------
        bool
            True if the object describes a position in the UTR; else False.

        """
        return self.utr is not None

    def is_intronic(self) -> bool:
        """Return whether this is an intronic position.

        Returns
        -------
        bool
            True if the object describes a position in an intron; else False.

        """
        return self.intronic_position is not None

    def is_extended(self) -> bool:
        """Return whether this position was described using the extended syntax.

        Returns
        -------
        bool
            True if the position was described using the extended syntax; else False.

        """
        return self.is_intronic() or self.is_utr()

    # the string annotation used in the type hint below is required for Python 3.6 compatibility
    def is_adjacent(self, other: "VariantPosition") -> bool:
        """Return whether this variant and another are immediately adjacent in sequence space.

        The following special cases are not handled correctly:

        * The special case involving the last variant in a transcript sequence and the first base in the 3'
          UTR will be evaluated as not adjacent, as the object does not have sequence length information.
        * The special case involving the two middle bases in an intron where the numbering switches from
          positive with respect to the 5' end of the intron to negative with respect to the 3' end of the intron
          will be evaluated as not adjacent, as the object does not have intron length information.
        * This ignores the special case where there is an intron between the last base of the 5' UTR and the first
          base of the coding sequence because it is not biologically relevant to the best of my knowledge.

        Parameters
        ----------
        other : VariantPosition
            The object to calculate adjacency to.

        Returns
        -------
        bool
            True if the positions describe adjacent bases in sequence space; else False.

        """
        if self.utr == other.utr:
            if self.intronic_position is None and other.intronic_position is None:
                return abs(self.position - other.position) == 1
            elif (
                self.position == other.position
            ):  # intronic positions can only be adjacent if they are relative to the same base
                if (
                    self.intronic_position is not None
                    and other.intronic_position is not None
                ):
                    return abs(self.intronic_position - other.intronic_position) == 1
                else:  # handle special case for first/last base of intron and corresponding first/last base of exon
                    return (
                        self.intronic_position == -1
                        or self.intronic_position == 1
                        or other.intronic_position == -1
                        or other.intronic_position == 1
                    )
            else:
                return False
        else:  # handle special case for last base of 5' utr and first base of non-UTR sequence
            return (self.position == -1 and other.position == 1) or (
                other.position == -1 and self.position == 1
            )

    def __lt__(self, other: "VariantPosition") -> bool:
        """Less than comparison operator.

        Other comparison operators will be filled in using :py:func:`functools.total_ordering`.

        Parameters
        ----------
        other : VariantPosition
            The other VariantPosition to compare to.

        Returns
        -------
        bool
            True if this position evaluates as strictly less than the other position; else False.

        """
        if self.utr == other.utr:
            if self.position == other.position:
                if self.intronic_position == other.intronic_position:
                    return False
                elif self.intronic_position is None:
                    return other.intronic_position < 0
                elif other.intronic_position is None:
                    return self.intronic_position > 0
            else:
                return self.position < other.position
        else:  # 5' < non-UTR < 3'
            if self.utr:
                if self.position < 0:  # self is in 5' UTR
                    return True
                else:  # self is in 3' UTR
                    return False
            else:
                if other.position < 0:  # other is in 5' UTR
                    return False
                else:  # other is in 3' UTR
                    return True

    def __eq__(self, other: "VariantPosition") -> bool:
        """Equality comparison operator.

        Other comparison operators will be filled in using :py:func:`functools.total_ordering`.

        Parameters
        ----------
        other : VariantPosition
            The other VariantPosition to compare to.

        Returns
        -------
        bool
            True if this position is the same as the other position; else False.

        """
        return (self.position, self.intronic_position, self.utr) == (
            other.position,
            other.intronic_position,
            other.utr,
        )

    def __ne__(self, other: "VariantPosition") -> bool:
        """Not equal comparison operator.

        Other comparison operators will be filled in using :py:func:`functools.total_ordering`.

        Parameters
        ----------
        other : VariantPosition
            The other VariantPosition to compare to.

        Returns
        -------
        bool
            True if this position is not the same as the other position; else False.

        """
        return (self.position, self.intronic_position, self.utr) != (
            other.position,
            other.intronic_position,
            other.utr,
        )
