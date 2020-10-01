from typing import Optional, Union, List, Tuple, Sequence, TypeVar
from mavehgvs.position import VariantPosition

# TODO: revisit this naming convention
VariantSequence = TypeVar("VariantSequence", str, Tuple[str, str])
"""Type variable for a variant sequence.

The sequence can be a string or a pair of strings for reference and new bases in a substitution variant.
"""


class Variant:
    def __init__(self, s: str):
        """Convert a MAVE-HGVS variant string into a corresponding object with named fields.

        Parameters
        ----------
        s : str
            MAVE-HGVS variant string to convert into an object.

        """
        self.variant_string = s
        self.matchdict = None
        self.validation_failure_message = None
        self.variant_count = None
        self._position = None
        self._prefix = None
        self._reference_id = None
        self._sequence = None
        self._variant_type = None

    @staticmethod
    def _position_distance(pos1: Union[int, str], pos2: Union[int, str]) -> int:
        """Calculate the minimum distance between two sequence positions.

        Some distances between positions using the extended notation cannot be calculated exactly without knowing the
        full target sequence.
        In such cases this function returns the minimum possible distance.

        Parameters
        ----------
        pos1 : Union[int, str]
            The first position, either as an integer or a string using the extended position notation.
        pos2 : Union[int, str]
            The second position, either as an integer or a string using the extended position notation.

        Returns
        -------
        int
            The minimum distance between positions.

        Raises
        ------
        ValueError
            If an integer position is not a positive number.
        ValueError
            If a string position is not valid.

        """
        if isinstance(pos1, int) and isinstance(pos2, int):
            return pos2 - pos1 + 1
        else:
            return None  # TODO: implement this

    @staticmethod
    def _positions_are_sorted(positions: Sequence[Union[int, str]]) -> bool:
        """Determine whether positions are in sorted order according to the 3' rule (ascending order).

        Parameters
        ----------
        positions : Sequence[Union[int, str]]
            The positions to check the order of.

        Returns
        -------
        bool
            True if the positions are sorted; else False.

        """
        pass

    @staticmethod
    def _substitution_matches_reference(pos: int, ref: str, target: str) -> bool:
        """Determine whether the reference portion of a substitution matches the target sequence.

        # TODO: this needs to be aware of protein vs nucleotide references

        Parameters
        ----------
        pos :
        ref
        target

        Returns
        -------

        """
        pass

    def is_valid(self) -> bool:
        """Return whether the variant is considered valid MAVE-HGVS.

        Returns
        -------
        bool
            True if the variant string is valid MAVE-HGVS; else False.

        """
        if self.validation_failure_message is None:
            return False
        else:
            return True

    def is_multi_variant(self) -> Optional[bool]:
        """Return whether the variant is a multi-variant.

        A multi-variant is a single variant describing multiple events enclosed in '[]'.
        Multi-variants are referred to as alleles in the HGVS standard.

        Returns
        -------
        Optional[bool]
            True if the variant is a multi-variant; else False. Returns None if the variant is invalid.

        """
        if not self.is_valid():
            return None
        elif self.variant_count > 1:
            return True
        else:
            return False

    @property
    def prefix(self) -> Optional[str]:
        """The single-letter prefix for this variant.

        Returns
        -------
        Optional[str]
            Single-letter prefix corresponding to the sequence type or None of the variant is invalid.

            See the following table for sequence type prefixes and their meanings:

            .. csv-table::
               :file: ../docs/prefix.csv
               :header: "Prefix", "Description"
               :widths: 5, 20

        """
        if not self.is_valid():
            return None
        else:
            return self._prefix

    @property
    def variant_types(self) -> Optional[Union[str, List[str]]]:
        """The type for this variant.

        Valid variant types are:

        * ``'substitution'``
        * ``'deletion'``
        * ``'duplication'``
        * ``'insertion'``
        * ``'deletion-insertion'``

        Returns
        -------
        Optional[Union[str, List[str]]]
            Single-word string containing the variant type or None of the variant is invalid.
            Returns a list of single-word strings for a multi-variant.

        """
        if not self.is_valid():
            return None
        else:
            return self._variant_types

    def uses_extended_positions(self) -> Optional[bool]:
        """Return whether the variant uses the extended position notation to describe intronic or UTR positions.

        Examples of variants using the extended position notation include:

        * c.122-6T>A
        * r.*33a>c
        * c.43-6_595+12delinsCTT

        This should always be false for variants with a genomic or protein prefix, as variants with these prefixes
        cannot use positions relative to a transcript under the MAVE-HGVS specification.

        Returns
        -------
        Optional[bool]
            True if the variant (or any of the individual variants for a multi-variant) uses the extended position
            notation or None if the variant is invalid.

        """
        if not self.is_valid():
            return None
        elif self.is_multi_variant():
            all_positions = list()
            for p in self.positions:
                if isinstance(p, tuple):
                    all_positions.extend(p)
                else:
                    all_positions.append(p)
            return any(p.is_extended() for p in all_positions)
        else:
            if isinstance(self.positions, tuple):
                return any(p.is_extended() for p in self.positions)
            else:
                return self.positions.is_extended()

    @property
    def positions(
        self
    ) -> Optional[
        Union[VariantPosition, Tuple[VariantPosition, VariantPosition], List[Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]]]
    ]:
        """The variant position as a single position or tuple containing start and end positions.

        Each position is an instance of :py:class:`mavehgvs.position.VariantPosition`.

        Returns
        -------
        Union[VariantPosition, Tuple[VariantPosition, VariantPosition], List[Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]]]
            Variant position or tuple of start/end positions, or None of the variant is invalid.
            Returns a list of positions or start/end tuples for a multi-variant.

        """
        if not self.is_valid():
            return None
        else:
            return self._positions

    @property
    def sequence(
        self
    ) -> Optional[Union[str, Tuple[str, str], List[Optional[Union[str, Tuple[str, str]]]]]]:
        """The sequence portion of the variant.

        This can be a tuple of reference and new bases for a substitution, a single sequence for insertions or
        deletion-insertions, or the "=" character for variants that are identical to the target sequence.

        Returns
        -------
        Union[str, Tuple[str, str], List[Optional[Union[str, Tuple[str, str]]]]]]
            Tuple of ref/new bases for substitutions, string containing inserted sequence, or the "=" character.
            Returns None if the variant is invalid or does not have a sequence component (deletion or duplication).
            Returns a list for a multi-variant, which may contain None values for deletions or duplications.

        """
        if not self.is_valid():
            return None
        else:
            return self._sequences

    @property
    def reference_id(self) -> Optional[str]:
        """The reference identifier for the variant (if applicable).

        The reference identifier precedes the prefix and is followed by a ``:``.
        For example in ``NM_001130145.3:c.832C>T`` the reference identifier is "NM_001130145.3".

        Returns
        -------
        Optional[str]
            The reference identifier, or None if it is not set or the variant is invalid.

        """
        if not self.is_valid():
            return None
        else:
            return self._reference_id
