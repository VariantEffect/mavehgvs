import re
from typing import Optional, Union, List, Tuple
from mavehgvs.position import VariantPosition
from mavehgvs.patterns.combined import any_variant


class Variant:
    __variant_fullmatch = re.compile(any_variant, flags=re.ASCII).fullmatch
    """Callable[[str, int, int], Optional[Match[str]]]: fullmatch callable for parsing a single MAVE-HGVS variant
    
    Returns an :py:obj:`re.Match` object if the full string defines a valid MAVE-HGVS variant.
    Match groups in the result can be used to extract components of the variant.
    """

    __vtypes = (
        "sub",     # substitution
        "del",     # deletion
        "dup",     # duplication
        "ins",     # insertion
        "delins",  # deletion-insertion"
    )
    """Tuple[str]: variant type tags used in MAVE-HGVS patterns and variant type names.
    """

    def __init__(self, s: str, targetseq: Optional[str] = None, relaxed_ordering: bool = False):
        """Convert a MAVE-HGVS variant string into a corresponding object with named fields.

        Parameters
        ----------
        s : str
            MAVE-HGVS variant string to convert into an object.

        targetseq : Optional[str]
            If provided, the variant will be validated for agreement with this sequence.
            Target sequence validation is not supported for variants using the extended position syntax.

            The type of the target sequence (DNA, RNA, or amino acid) will be inferred.
            DNA and amino acid sequences should be in uppercase, RNA in lowercase.

        relaxed_ordering : bool
            If True, variant strings that do not observe the 3-prime rule for variant position ordering are allowed.
            The object representation will observe the 3-prime rule, so it may differ from the input string in this
            case.

        """
        self.variant_string = s
        variant_match = self.__variant_fullmatch(s)
        if variant_match is None:
            self.validation_failure_message = "failed regular expression validation"
        else:
            self.validation_failure_message = None
            self._groupdict = variant_match.groupdict()

            # set reference id if present
            if self._groupdict["reference_id"] is not None:
                self._reference_id = self._groupdict["reference_id"]
            else:
                self._reference_id = None

            # set prefix and determine if this is a multi-variant
            if self._groupdict["single_variant"] is not None:
                self.variant_count = 1
                self._prefix = self._groupdict["single_variant"][0]
            elif self._groupdict["multi_variant"] is not None:
                self.variant_count = len(s.split(";"))
                self._prefix = self._groupdict["multi_variant"][0]
            else:  # pragma: no cover
                raise ValueError("invalid match type")

            if self.variant_count == 1:
                # determine which named groups to check
                if self._prefix == "p":
                    gdict_prefixes = [(f"pro_{t}", t) for t in self.__vtypes]
                elif self._prefix == "r":
                    gdict_prefixes = [(f"rna_{t}", t) for t in self.__vtypes]
                elif self._prefix in "cn":
                    gdict_prefixes = [
                        (f"dna_{t}_{self._prefix}", t) for t in self.__vtypes
                    ]
                elif self._prefix in "gmo":
                    gdict_prefixes = [(f"dna_{t}_gmo", t) for t in self.__vtypes]
                else:  # pragma: no cover
                    raise ValueError("unexpected prefix")

                # set the variant type
                self._variant_types = None
                vtype_set = False
                groupdict_prefix = None
                for groupname, vtype in gdict_prefixes:
                    if self._groupdict[groupname] is not None:
                        if vtype_set:  # pragma: no cover
                            raise ValueError(f"ambiguous match: '{groupname}' and '{groupdict_prefix}'")
                        self._variant_types = vtype
                        groupdict_prefix = groupname

                # set the position and sequence
                self._positions = None
                self._sequences = None
                if self._variant_types == "sub":
                    if self._groupdict[f"{groupdict_prefix}_equal"] is not None:  # special case for target identity
                        self._sequences = self._groupdict[f"{groupdict_prefix}_equal"]
                    else:
                        self._positions = VariantPosition(self._groupdict[f"{groupdict_prefix}_position"])
                        if self._prefix == "p":
                            self._sequences = (self._positions.amino_acid, self._groupdict[f"{groupdict_prefix}_new"])
                        elif self._prefix in "gmo" "cn" "r":
                            self._sequences = (self._groupdict[f"{groupdict_prefix}_ref"], self._groupdict[f"{groupdict_prefix}_new"])
                        else:  # pragma: no cover
                            raise ValueError("unexpected prefix")
                elif self._variant_types in ("del", "dup", "ins", "delins"):
                    # set position
                    if self._groupdict[f"{groupdict_prefix}_pos"] is not None:
                        self._positions = VariantPosition(self._groupdict[f"{groupdict_prefix}_pos"])
                    else:
                        self._positions = (VariantPosition(self._groupdict[f"{groupdict_prefix}_start"]), VariantPosition(self._groupdict[f"{groupdict_prefix}_end"]))
                        # extra validation on positions
                        if self._positions[0] >= self._positions[1]:
                            raise ValueError("start position must be before end position")
                        if self._variant_types == "ins":
                            if not self._positions[0].is_adjacent(self._positions[1]):
                                raise ValueError("insertion positions must be adjacent")

                    # set sequence if needed
                    if self._variant_types in ("ins", "delins"):
                        self._sequences = self._groupdict[f"{groupdict_prefix}_seq"]
                else:  # pragma: no cover
                    raise ValueError("unexpected variant type")

    def __repr__(self) -> str:
        """The object representation is equivalent to the input string.

        Returns
        -------
        str
            The object representation.

        """
        def format_variant(vtype: str, pos: Union[VariantPosition, Tuple[VariantPosition, VariantPosition]], seq: Optional[Union[str, Tuple[str, str]]]) -> str:
            """Helper function for building variant strings.

            Parameters
            ----------
            vtype : str
                The variant type, as described by :py:obj:`Variant.__vtypes`
            pos : Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]
                The position or pair of positions describing the variant.
            seq : Optional[Union[str, Tuple[str, str]]]
                The sequence or pair of sequences describing the variant.
                Only used for substitions, insertions, and deletion-insertions.

            Returns
            -------
            str
                A string representing this variant element.

            """
            if vtype == "sub":
                if self._prefix == "p":  # protein variant
                    return f"{pos}{seq[1]}"
                else:                    # nucleotide variant
                    return f"{pos}{seq[0]}>{seq[1]}"
            elif vtype in ("del", "dup"):
                if isinstance(pos, tuple):
                    return f"{pos[0]}_{pos[1]}{vtype}"
                else:
                    return f"{pos}{vtype}"
            elif vtype == ("ins", "delins"):
                if isinstance(pos, tuple):
                    return f"{pos[0]}_{pos[1]}{vtype}{seq}"
                else:
                    return f"{pos}{vtype}{seq}"
            else:   # pragma: no cover
                raise ValueError("invalid variant type")

        if self._reference_id is not None:
            prefix = f"{self._reference_id}:{self._prefix}"
        else:
            prefix = f"{self._prefix}"

        if not self.is_valid():
            return repr(None)
        elif self.is_target_identical():
            return f"{prefix}.="
        elif self.variant_count > 1:
            elements = list()
            for vtype, pos, seq in zip(self._variant_types, self._positions, self._sequences):
                elements.append(format_variant(vtype, pos, seq))
            return f"{prefix}.[{';'.join(elements)}]"
        else:
            return f"{prefix}.{format_variant(self._variant_types, self._positions, self._sequences)}"

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
        return self.validation_failure_message is None

    def is_target_identical(self) -> Optional[bool]:
        """Return whether the variant describes the "wild-type" sequence.

        This is the variant described with only the equals sign (e.g. ``c.=``).

        Returns
        -------
        Optional[bool]
            True if this variant describes the wild-type or target sequence; else False.
            Returns None if the variant is invalid.

        """
        if not self.is_valid():
            return None
        else:
            return self._positions is None

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

        * ``'sub'`` for substitutions
        * ``'del'`` for deletions
        * ``'dup'`` for duplications
        * ``'ins'`` for insertions
        * ``'delins'`` for deletion-insertions

        Returns
        -------
        Optional[Union[str, List[str]]]
            String containing the variant type or None of the variant is invalid.
            Returns a list of strings for a multi-variant.

        """
        if not self.is_valid():
            return None
        else:
            if self.is_multi_variant():
                return [self.__vtypes[t] for t in self._variant_types]
            else:
                return self.__vtypes[self._variant_types]

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
        Union[
            VariantPosition,
            Tuple[VariantPosition, VariantPosition],
            List[Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]],
        ]
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
    ) -> Optional[
        Union[str, Tuple[str, str], List[Optional[Union[str, Tuple[str, str]]]]]
    ]:
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
