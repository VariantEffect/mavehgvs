import re
import itertools
from typing import Optional, Union, List, Tuple, Mapping, Any, Sequence, Dict, Generator

from fqfa.constants import AA_CODES

from mavehgvs.position import VariantPosition
from mavehgvs.patterns.combined import any_variant
from mavehgvs.exceptions import MaveHgvsParseError

__all__ = ["Variant"]

AA_3_TO_1 = {value: key for key, value in AA_CODES.items()}
"""Dict[str, str]: for converting three-letter amino acid codes to single-letter codes.
"""


class Variant:
    fullmatch = re.compile(any_variant, flags=re.ASCII).fullmatch
    """Callable[[str, int, int], Optional[Match[str]]]: fullmatch callable for parsing a single MAVE-HGVS variant

    Returns an :py:obj:`re.Match` object if the full string defines a valid MAVE-HGVS variant.
    Match groups in the result can be used to extract components of the variant.
    """

    VTYPES = (
        "equal",  # equality
        "sub",  # substitution
        "fs",  # frame shift
        "del",  # deletion
        "dup",  # duplication
        "ins",  # insertion
        "delins",  # deletion-insertion
    )
    """Tuple[str]: variant type tags used in MAVE-HGVS patterns and variant type names.
    """

    def __init__(
        self,
        s: Union[str, Mapping[str, Any], Sequence[Mapping[str, Any]]],
        targetseq: Optional[str] = None,
        relaxed_ordering: bool = False,
    ):
        """Convert a MAVE-HGVS variant string into a corresponding object with named fields.

        Parameters
        ----------
        s : Union[str, Mapping[str, Any], Sequence[Mapping[str, Any]]]
            MAVE-HGVS variant string to convert into an object, dictionary type object containing key-value pairs
            corresponding to a MAVE-HGVS object, or list/tuple of dictionary type objects for a variant with
            multiple events.

        targetseq : Optional[str]
            If provided, the variant will be validated for agreement with this sequence.
            Target sequence validation is not supported for variants using the extended position syntax.

            This must be an amino acid sequence for protein variants or a nucleotide sequence for
            coding/noncoding/genomic variants.
            DNA and amino acid sequences should be in uppercase, RNA in lowercase.

        relaxed_ordering : bool
            If True, variants that do not observe the 3-prime rule for variant position ordering are allowed.
            The object representation will observe the 3-prime rule, so it may differ from the input string in this
            case.

        """
        if isinstance(s, str):  # variant string to parse
            variant_string = s
        elif isinstance(s, Mapping):  # dictionary-style single variant
            variant_string = self._variant_dictionary_to_string(s, include_prefix=True)
        elif isinstance(s, Sequence):  # dictionary-style multi-variant
            if not all(isinstance(v, Mapping) for v in s):
                raise ValueError("multi-variant iterable must contain Mapping objects")
            try:
                all_prefixes = [v["prefix"] for v in s]
            except KeyError:
                raise MaveHgvsParseError("variant dictionary missing required keys")
            if len(set(all_prefixes)) != 1:
                raise MaveHgvsParseError(
                    "cannot combine variants with different prefixes"
                )
            variant_string = f"{s[0]['prefix']}.[{';'.join(self._variant_dictionary_to_string(v, include_prefix=False) for v in s)}]"
        else:
            raise ValueError("can only create Variants from string or Mapping objects")

        variant_match = self.fullmatch(variant_string)
        if variant_match is None:
            raise MaveHgvsParseError("failed regular expression validation")
        else:
            match_dict = variant_match.groupdict()

            # set target id if present
            if match_dict["target_id"] is not None:
                self._target_id = match_dict["target_id"]
            else:
                self._target_id = None

            # set prefix and determine if this is a multi-variant
            if match_dict["single_variant"] is not None:
                self.variant_count = 1
                self._prefix = match_dict["single_variant"][0]
            elif match_dict["multi_variant"] is not None:
                self.variant_count = len(variant_string.split(";"))
                self._prefix = match_dict["multi_variant"][0]
            else:  # pragma: no cover
                raise ValueError("invalid match type")

            if self.variant_count == 1:
                (
                    self._variant_types,
                    self._positions,
                    self._sequences,
                ) = self._process_string_variant(
                    match_dict, relaxed_ordering=relaxed_ordering
                )
            elif self.variant_count > 1:
                self._variant_types = list()
                self._positions = list()
                self._sequences = list()

                # format each individual variant event as a single variant and parse it
                for variant_substring in match_dict["multi_variant"][3:-1].split(";"):
                    groupdict = self.fullmatch(
                        f"{self._prefix}.{variant_substring}"
                    ).groupdict()
                    vt, p, s = self._process_string_variant(
                        groupdict, relaxed_ordering=relaxed_ordering
                    )
                    if vt == "equal":
                        raise MaveHgvsParseError(
                            "multi-variants cannot contain target-identical variants"
                        )

                    self._variant_types.append(vt)
                    self._positions.append(p)
                    self._sequences.append(s)

                # ensure that multiple variants aren't defined for the same positions
                for vp1, vp2 in itertools.combinations(self._positions, 2):
                    if isinstance(vp1, VariantPosition) and isinstance(
                        vp2, VariantPosition
                    ):  # both single position
                        if vp1 == vp2:
                            raise MaveHgvsParseError(
                                "multi-variant has multiple changes for the same position"
                            )
                    elif isinstance(vp1, VariantPosition) and isinstance(vp2, Tuple):
                        if vp2[0] <= vp1 <= vp2[1]:
                            raise MaveHgvsParseError(
                                "multi-variant has overlapping changes"
                            )
                    elif isinstance(vp1, Tuple) and isinstance(vp2, VariantPosition):
                        if vp1[0] <= vp2 <= vp1[1]:
                            raise MaveHgvsParseError(
                                "multi-variant has overlapping changes"
                            )
                    elif isinstance(vp1, Tuple) and isinstance(vp2, Tuple):
                        if (
                            vp1[0] <= vp2[0] <= vp1[1]
                            or vp1[0] <= vp2[1] <= vp1[1]
                            or vp2[0] <= vp1[0] <= vp2[1]
                            or vp2[0] <= vp1[1] <= vp2[1]
                        ):
                            raise MaveHgvsParseError(
                                "multi-variant has overlapping changes"
                            )
                    else:  # pragma: no cover
                        raise ValueError("invalid position type")

                # re-order variants and validate
                def sort_key(x):
                    if isinstance(x[1], VariantPosition):
                        return x[1]
                    elif isinstance(x[1], Tuple):
                        return x[1][0]
                    else:  # pragma: no cover
                        raise ValueError("invalid position type")

                variant_list = list(self.variant_tuples())
                ordered_list = sorted(variant_list, key=sort_key)
                if variant_list != ordered_list:
                    if relaxed_ordering:
                        self._variant_types = [x[0] for x in ordered_list]
                        self._positions = [x[1] for x in ordered_list]
                        self._sequences = [x[2] for x in ordered_list]
                    else:
                        raise MaveHgvsParseError("multi-variants not in sorted order")

                # make sure there is at most one frame shift
                if sum(x == "fs" for x in self._variant_types) > 1:
                    raise MaveHgvsParseError("maximum of one frame shift is permitted")

                # make sure the frame shift is last if present
                if any(x == "fs" for x in self._variant_types):
                    if self._variant_types[-1] != "fs":
                        raise MaveHgvsParseError(
                            "no variants are permitted to follow a frame shift"
                        )

            else:  # pragma: no cover
                raise ValueError("invalid variant count")

        if targetseq is not None:
            for vtype, pos, seq in self.variant_tuples():
                if self._prefix != "p" and vtype == "sub":
                    self._target_validate(pos, seq[0], targetseq)
                elif (
                    pos is None and vtype == "equal"
                ):  # special case for full-length target identical variants
                    pass
                else:
                    self._target_validate(pos, None, targetseq)

    def variant_tuples(
        self,
    ) -> Generator[
        Tuple[
            str,
            Optional[Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]],
            Optional[Union[str, Tuple[str, str]]],
        ],
        None,
        None,
    ]:
        """Generator that yields tuples containing the variant components.

        Yields
        ------
        Tuple
            Tuple of the variant type, position(s), and sequence(s) for each element in the variant.

        """
        if self.is_multi_variant():
            for vtype, pos, seq in zip(
                self._variant_types, self._positions, self._sequences
            ):
                yield vtype, pos, seq
        else:
            yield self._variant_types, self._positions, self._sequences

    def _process_string_variant(
        self, match_dict: Dict[str, str], relaxed_ordering: bool
    ) -> Tuple[
        str,
        Optional[Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]],
        Optional[Union[str, Tuple[str, str]]],
    ]:
        """Process the match dictionary from a single variant into its components.

        Parameters
        ----------
        match_dict : Dict[str, str]
            Match dictionary from the MAVE-HGVS regular expression.
        relaxed_ordering : bool
            If True, variants that do not observe the 3-prime rule for variant position ordering are allowed.

        Returns
        -------
        Tuple[str, Optional[Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]], Optional[Union[str, Tuple[str, str]]]]
            Returns a 3-tuple containing the variant type, optional position (or start/end positions),
            and optional before/after substitution sequences or inserted sequence.

        """
        variant_type = None
        positions = None
        sequences = None

        # determine which named groups to check
        if self._prefix == "p":
            pattern_group_tuples = [(f"pro_{t}", t) for t in self.VTYPES]
        elif self._prefix == "r":
            pattern_group_tuples = [(f"rna_{t}", t) for t in self.VTYPES if t != "fs"]
        elif self._prefix in tuple("cn"):
            pattern_group_tuples = [
                (f"dna_{t}_{self._prefix}", t) for t in self.VTYPES if t != "fs"
            ]
        elif self._prefix in tuple("gmo"):
            pattern_group_tuples = [
                (f"dna_{t}_gmo", t) for t in self.VTYPES if t != "fs"
            ]
        else:  # pragma: no cover
            raise ValueError("unexpected prefix")

        # set the variant type
        vtype_set = False
        pattern_group = None
        for pg, vtype in pattern_group_tuples:
            if match_dict[pg] is not None:
                if vtype_set:  # pragma: no cover
                    raise ValueError(f"ambiguous match: '{pg}' and '{pattern_group}'")
                variant_type = vtype
                pattern_group = pg
                vtype_set = True

        # set the position and sequence
        if variant_type == "sub":
            positions = VariantPosition(match_dict[f"{pattern_group}_position"])
            if self._prefix == "p":
                sequences = (positions.amino_acid, match_dict[f"{pattern_group}_new"])
            elif self._prefix in tuple("gmo" "cn" "r"):
                sequences = (
                    match_dict[f"{pattern_group}_ref"],
                    match_dict[f"{pattern_group}_new"],
                )
            else:  # pragma: no cover
                raise ValueError("unexpected prefix")
        elif variant_type in ("equal", "fs", "del", "dup", "ins", "delins"):
            # set position
            if (
                match_dict.get(f"{pattern_group}_position") is not None
            ):  # use get() since ins pattern doesn't have pos
                positions = VariantPosition(match_dict[f"{pattern_group}_position"])
            elif (
                match_dict.get(f"{pattern_group}_start") is not None
                and match_dict.get(f"{pattern_group}_end") is not None
            ):
                positions = (
                    VariantPosition(match_dict[f"{pattern_group}_start"]),
                    VariantPosition(match_dict[f"{pattern_group}_end"]),
                )
                # extra validation on positions
                if positions[0] >= positions[1]:
                    if relaxed_ordering:
                        positions = (positions[1], positions[0])
                    else:
                        raise MaveHgvsParseError(
                            "start position must be before end position"
                        )
                if variant_type == "ins":
                    if not positions[0].is_adjacent(positions[1]):
                        raise MaveHgvsParseError("insertion positions must be adjacent")
            else:  # pragma: no cover
                if variant_type != "equal":
                    raise MaveHgvsParseError("variant position not found")

            # set sequence if needed
            if variant_type in ("ins", "delins"):
                sequences = match_dict[f"{pattern_group}_seq"]
            elif variant_type == "equal":
                if (
                    match_dict[f"{pattern_group}_equal"] is not None
                ):  # special case for target identity
                    sequences = match_dict[f"{pattern_group}_equal"]
                elif match_dict[f"pro_equal_equal_sy"] is not None:
                    sequences = match_dict[f"pro_equal_equal_sy"]

        return variant_type, positions, sequences

    # TODO: API documentation for the dictionary objects
    @staticmethod
    def _variant_dictionary_to_string(
        vdict: Mapping[str, Any], include_prefix: bool
    ) -> str:
        """Convert a match dictionary from a single variant into a string for further validation.

        This method performs minimal validation of the values provided in the input, and instead converts it into a
        variant string that is validated using the regular expression based validators.

        Parameters
        ----------
        vdict : Mapping[str, Any]
            Key-value pairs describing a single variant.
        include_prefix: bool
            If True, the variant prefix and '.' will be included in the string; else it is omitted (for use with
            multi-variants).

        Returns
        -------
        str
            A string representing this variant.

        Raises
        ------
        MaveHgvsParseError
            If the dictionary does not have a valid set of keys.

        """
        try:
            variant_type = vdict["variant_type"]
            prefix = vdict["prefix"]
        except KeyError:
            raise MaveHgvsParseError("variant dictionary missing required keys")

        if variant_type == "equal":
            expected_keys = ["variant_type", "prefix"]
            if prefix == "p":
                expected_keys.extend(["position", "target"])
            else:
                expected_keys.extend(["start_position", "end_position"])
            if sorted(vdict.keys()) != sorted(expected_keys):
                raise MaveHgvsParseError("variant dictionary contains invalid keys")
            if prefix == "p":
                variant_string = f"{vdict['target']}{vdict['position']}="
            elif vdict["start_position"] == vdict["end_position"]:
                variant_string = f"{vdict['start_position']}="
            else:
                variant_string = f"{vdict['start_position']}_{vdict['end_position']}="
        elif variant_type == "sub":
            if sorted(vdict.keys()) != sorted(
                ["variant_type", "prefix", "position", "target", "variant"]
            ):
                raise MaveHgvsParseError("variant dictionary contains invalid keys")
            if prefix == "p":
                variant_string = (
                    f"{vdict['target']}{vdict['position']}{vdict['variant']}"
                )
            else:
                variant_string = (
                    f"{vdict['position']}{vdict['target']}>{vdict['variant']}"
                )
        elif variant_type == "fs":
            if sorted(vdict.keys()) != sorted(
                ["variant_type", "prefix", "position", "target"]
            ):
                raise MaveHgvsParseError("variant dictionary contains invalid keys")
            if prefix == "p":
                variant_string = f"{vdict['target']}{vdict['position']}fs"
            else:
                raise MaveHgvsParseError(
                    "frame shifts are only supported for protein variants"
                )
        elif variant_type in ("del", "dup"):
            expected_keys = ["variant_type", "prefix", "start_position", "end_position"]
            if prefix == "p":
                expected_keys.extend(["start_target", "end_target"])
            if sorted(vdict.keys()) != sorted(expected_keys):
                raise MaveHgvsParseError("variant dictionary contains invalid keys")
            if prefix == "p":
                start = f"{vdict['start_target']}{vdict['start_position']}"
                end = f"{vdict['end_target']}{vdict['end_position']}"
            else:
                start = vdict["start_position"]
                end = vdict["end_position"]
            if start == end:
                variant_string = f"{start}{variant_type}"
            else:
                variant_string = f"{start}_{end}{variant_type}"
        elif variant_type in ("ins", "delins"):
            expected_keys = [
                "variant_type",
                "prefix",
                "start_position",
                "end_position",
                "variant",
            ]
            if prefix == "p":
                expected_keys.extend(["start_target", "end_target"])
            if sorted(vdict.keys()) != sorted(expected_keys):
                raise MaveHgvsParseError("variant dictionary contains invalid keys")
            if prefix == "p":
                start = f"{vdict['start_target']}{vdict['start_position']}"
                end = f"{vdict['end_target']}{vdict['end_position']}"
            else:
                start = vdict["start_position"]
                end = vdict["end_position"]
            variant_string = f"{start}_{end}{variant_type}{vdict['variant']}"
        else:
            raise MaveHgvsParseError("invalid variant type")

        if include_prefix:
            return f"{vdict['prefix']}.{variant_string}"
        else:
            return variant_string

    def __eq__(self, other: "Variant") -> bool:
        """Equality comparison operator.

        Parameters
        ----------
        other : Variant
            The other Variant to compare to.

        Returns
        -------
        bool
            True if this variant is the same as the other position; else False.

        """
        return (
            self._target_id,
            self.variant_count,
            self._prefix,
            self._variant_types,
            self._positions,
            self._sequences,
        ) == (
            other._target_id,
            other.variant_count,
            other._prefix,
            other._variant_types,
            other._positions,
            other._sequences,
        )

    def __repr__(self) -> str:
        """The object representation is equivalent to the input string.

        Returns
        -------
        str
            The object representation.

        """

        def format_variant(
            vtype: str,
            pos: Union[VariantPosition, Tuple[VariantPosition, VariantPosition]],
            seq: Optional[Union[str, Tuple[str, str]]],
        ) -> str:
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
                else:  # nucleotide variant
                    return f"{pos}{seq[0]}>{seq[1]}"
            elif vtype == "fs":
                return f"{pos}fs"
            elif vtype in ("del", "dup"):
                if isinstance(pos, tuple):
                    return f"{pos[0]}_{pos[1]}{vtype}"
                else:
                    return f"{pos}{vtype}"
            elif vtype in ("ins", "delins"):
                if isinstance(pos, tuple):
                    return f"{pos[0]}_{pos[1]}{vtype}{seq}"
                else:
                    return f"{pos}{vtype}{seq}"
            elif vtype == "equal":
                if pos is None:
                    return f"{seq}"
                elif isinstance(pos, tuple):
                    return f"{pos[0]}_{pos[1]}{seq}"
                else:
                    return f"{pos}{seq}"
            else:  # pragma: no cover
                raise ValueError("invalid variant type")

        if self._target_id is not None:
            prefix = f"{self._target_id}:{self._prefix}"
        else:
            prefix = f"{self._prefix}"

        elements = [format_variant(*t) for t in self.variant_tuples()]
        if self.is_multi_variant():
            return f"{prefix}.[{';'.join(elements)}]"
        else:
            return f"{prefix}.{elements[0]}"

    @staticmethod
    def _target_validate(
        pos: Union[VariantPosition, Tuple[VariantPosition, VariantPosition]],
        ref: Optional[str],
        target: str,
    ) -> None:
        """Determine whether the target portion of a variant matches the target sequence.

        Note that variants using extended syntax cannot be validated with this method.
        If an extended syntax variant is encountered, it will be interpreted as valid/matching.

        Parameters
        ----------
        pos : Union[VariantPosition, Tuple[VariantPosition, VariantPosition]]
            Single variant position or start/end tuple for an indel.
        ref : Optional[str]
            Reference base to validate for nucleotide substitutions.
            This should be None for amino acid substitutions, since the reference is included in the VariantPosition.
        target : str
            Target sequence. This must be an amino acid sequence for protein variants or a nucleotide sequence
            for coding/noncoding/genomic variants.
            RNA sequences should be in lowercase, DNA sequences should be in uppercase.

        Returns
        -------
        None

        Raises
        ------
        MaveHgvsParseError
            If the reference base or amino acid does not match the target at the given position
        MaveHgvsParseError
            If the position is outside the bounds of the target.

        """
        if not isinstance(pos, tuple):
            pos = (pos,)

        if any(p.is_extended() for p in pos):
            return
        elif any(p.position > len(target) for p in pos):
            raise MaveHgvsParseError("variant coordinate out of bounds")
        else:
            if ref is not None and len(pos) == 1:  # nucleotide substitution
                if target[pos[0].position - 1] != ref:
                    raise MaveHgvsParseError("variant reference does not match target")
            elif pos[0].amino_acid is not None:  # protein variant
                for p in pos:
                    if target[p.position - 1] != AA_3_TO_1[p.amino_acid]:
                        raise MaveHgvsParseError(
                            "variant reference does not match target"
                        )
            else:
                return

    def is_target_identical(self) -> bool:
        """Return whether the variant describes the "wild-type" sequence or is the special synonymous variant.

        This is the variant described with only the equals sign (e.g. ``c.=``)
        or the uncertain equals protein variant (e.g. ``p.(=)``).

        Coding or genomic variants that specify an identical region (e.g. ``c.1_3=`` are also considered target
        identical.

        Synonymous protein variants (e.g. ``p.Leu12=``) are not considered target identical.

        Returns
        -------
        bool
            True if this variant describes the wild-type or target sequence; else False.

        """
        if self._variant_types == "equal":
            if self._prefix == "p":
                return self._positions is None
            else:
                return True
        else:
            return False

    def is_synonymous(self) -> bool:
        """Return whether the variant describes a synonymous protein variant or is the special synonymous variant.

        Returns
        -------
        bool
            True if this variant describes a synonymous protein variant; else False.

        """
        return self._variant_types == "equal" and self._prefix == "p"

    def is_multi_variant(self) -> bool:
        """Return whether the variant is a multi-variant.

        A multi-variant is a single variant describing multiple events enclosed in '[]'.
        Multi-variants are referred to as alleles in the HGVS standard.

        Returns
        -------
        bool
            True if the variant is a multi-variant; else False.

        """
        return self.variant_count > 1

    @property
    def prefix(self) -> str:
        """The single-letter prefix for this variant.

        Returns
        -------
        str
            Single-letter prefix corresponding to the sequence type.

            See the following table for sequence type prefixes and their meanings:

            .. csv-table::
               :file: ../docs/prefix.csv
               :header: "Prefix", "Description"
               :widths: 5, 20

        """
        return self._prefix

    @property
    def variant_type(self) -> Union[str, List[str]]:
        """The type for this variant.

        Valid variant types are:

        * ``'sub'`` for substitutions
        * ``'del'`` for deletions
        * ``'dup'`` for duplications
        * ``'ins'`` for insertions
        * ``'delins'`` for deletion-insertions

        Returns
        -------
        Union[str, List[str]]
            String containing the variant type. Returns a list of strings for a multi-variant.

        """
        return self._variant_types

    def uses_extended_positions(self) -> bool:
        """Return whether the variant uses the extended position notation to describe intronic or UTR positions.

        Examples of variants using the extended position notation include:

        * c.122-6T>A
        * r.*33a>c
        * c.43-6_595+12delinsCTT

        This should always be false for variants with a genomic or protein prefix, as variants with these prefixes
        cannot use positions relative to a transcript under the MAVE-HGVS specification.

        Returns
        -------
        bool
            True if the variant (or any of the individual variants for a multi-variant) uses the extended position
            notation.

        """
        if self.is_multi_variant():
            all_positions = list()
            for p in self.positions:
                if isinstance(p, tuple):
                    all_positions.extend(p)
                else:
                    all_positions.append(p)
            return any(p.is_extended() for p in all_positions)
        else:
            if self._positions is None:  # special case for target identity
                return False
            elif isinstance(self.positions, tuple):
                return any(p.is_extended() for p in self.positions)
            else:
                return self.positions.is_extended()

    @property
    def positions(
        self,
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
            Variant position or tuple of start/end positions.
            Returns a list of positions or start/end tuples for a multi-variant.

        """
        return self._positions

    @property
    def sequence(
        self,
    ) -> Optional[
        Union[str, Tuple[str, str], List[Optional[Union[str, Tuple[str, str]]]]]
    ]:
        """The sequence portion of the variant.

        This can be a tuple of target and new bases for a substitution, a single sequence for insertions or
        deletion-insertions, or the "=" character for variants that are identical to the target sequence.

        Returns
        -------
        Union[str, Tuple[str, str], List[Optional[Union[str, Tuple[str, str]]]]]]
            Tuple of ref/new bases for substitutions, string containing inserted sequence, or the "=" character.
            Returns None if the variant does not have a sequence component (deletion or duplication).
            Returns a list for a multi-variant, which may contain None values for deletions or duplications.

        """
        return self._sequences

    @property
    def target_id(self) -> Optional[str]:
        """The target identifier for the variant (if applicable).

        The target identifier precedes the prefix and is followed by a ``:``.
        For example in ``NM_001130145.3:c.832C>T`` the target identifier is "NM_001130145.3".

        Returns
        -------
        Optional[str]
            The target identifier, or None if it is not set.

        """
        return self._target_id
