from typing import Optional, Union, List, Tuple, TypeVar

VariantPosition = TypeVar(
    "VariantPosition", int, str, Tuple[Union[int, str], Union[int, str]]
)
"""Type variable for a variant position.

The position can be an int or string (for extended positions relative to transcripts) and contain either a single value
or a tuple of start and end coordinates for variants affecting a position range (e.g. insertions and deletions).
"""

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
        pass

    def is_valid(self) -> bool:
        """Return whether the variant is considered valid MAVE-HGVS.

        Returns
        -------
        bool
            True if the variant string is valid MAVE-HGVS; else False.

        """
        pass

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
        else:
            pass

    def prefix(self) -> Optional[str]:
        """Return the single-letter prefix for this variant.

        # TODO: consider whether to use properties for this and similar methods

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
            pass

    def variant_type(self) -> Optional[Union[str, List[str]]]:
        """Return the type for this variant.

        Valid variant types are:

        * :code:`'substitution'`
        * :code:`'deletion'`
        * :code:`'duplication'`
        * :code:`'insertion'`
        * :code:`'deletion-insertion'`

        # TODO: consider whether to use properties for this and similar methods

        Returns
        -------
        Optional[Union[str, List[str]]]
            Single-word string containing the variant type or None of the variant is invalid.
            Returns a list of single-word strings for a multi-variant.

        """
        if not self.is_valid():
            return None
        else:
            pass

    def uses_transcript_positions(self) -> Optional[bool]:
        """Return whether the variant uses the extended position notation to describe intronic or UTR positions.

        Examples of variants using the extended position notation include:

        * c.122-6T>A
        * r.*33a>c
        * c.43-6_595+12delinsCTT

        This will always be false for variants with a genomic or protein prefix, as variants with these prefixes cannot
        use positions relative to a transcript.

        Returns
        -------
        Optional[bool]
            True if the variant (or any of the individual variants for a multi-variant) uses the extended position
            notation or None if the variant is invalid.

        """
        if not self.is_valid():
            return None
        else:
            pass

    def position(self) -> Optional[Union[VariantPosition, List[VariantPosition]]]:
        """Returns the variant position as a single position or tuple containing start and end positions.

        Each position is represented as an integer or a string (for variants using extended position notation).

        Returns
        -------
        Optional[Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]]
            Variant position or tuple of start/end positions, or None of the variant is invalid.
            Returns a list of positions or start/end tuples for a multi-variant.

        """
        if not self.is_valid():
            return None
        else:
            pass

    def sequence(
        self
    ) -> Optional[Union[VariantSequence, List[Optional[VariantSequence]]]]:
        """Returns the sequence portion of the variant.

        This can be a tuple of reference and new bases for a substitution, a single sequence for insertions or
        deletion-insertions, or the "=" character for variants that are identical to the target sequence.

        Returns
        -------
        Optional[Union[VariantSequence, List[Optional[VariantSequence]]]]
            Tuple of ref/new bases for substitutions, string containing inserted sequence, or the "=" character.
            Returns None if the variant is invalid or does not have a sequence component (deletion or duplication).
            Returns a list for a multi-variant, which may contain None values for deletions or duplications.

        """
        if not self.is_valid():
            return None
        else:
            pass
