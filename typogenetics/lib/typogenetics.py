import logging
from dataclasses import dataclass
from enum import StrEnum, auto
from typing import Iterator, List, Optional, Tuple

logger = logging.getLogger(__name__)


class BaseType(StrEnum):
    PURINE = auto()
    PYRIMIDINE = auto()


class Base(StrEnum):
    C = auto()
    G = auto()
    T = auto()
    A = auto()

    @classmethod
    def from_str(cls, base_str: str) -> "Base":
        return {
            "C": cls.C,
            "G": cls.G,
            "T": cls.T,
            "A": cls.A,
        }[base_str]

    def __repr__(self) -> str:
        return self.value.upper()

    def __str__(self) -> str:
        return self.__repr__()

    def is_type(self, base_type: BaseType) -> bool:
        if base_type == BaseType.PURINE:
            return self.is_purine()
        return self.is_pyrimidine()

    def is_purine(self) -> bool:
        return self in [Base.A, Base.G]

    def is_pyrimidine(self) -> bool:
        return self in [Base.C, Base.T]

    def get_complement(self) -> "Base":
        return {
            Base.C: Base.G,
            Base.G: Base.C,
            Base.T: Base.A,
            Base.A: Base.T,
        }[self]


class AminoAcid(StrEnum):
    CUT = auto()
    DEL = auto()
    SWI = auto()
    MVR = auto()
    MVL = auto()
    COP = auto()
    OFF = auto()
    INA = auto()
    INC = auto()
    ING = auto()
    INT = auto()
    RPY = auto()
    RPU = auto()
    LPY = auto()
    LPU = auto()


Duplet = Tuple[Base, Base]


@dataclass
class Strand:
    bases: List[Base]

    @classmethod
    def from_str(cls, strand_str: str) -> "Strand":
        bases = []
        for base_str in strand_str:
            if base_str == " ":
                continue
            base = Base.from_str(base_str)
            bases.append(base)
        return cls(bases)

    def iter_bases(self) -> Iterator[Base]:
        yield from self.bases

    def iter_duplets(self) -> Iterator[Duplet]:
        unit = 0
        while True:
            if unit + 1 >= len(self):
                break

            yield (self[unit], self[unit + 1])

            unit += 2

    def __repr__(self) -> str:
        return "".join([str(b) for b in self.bases])

    def __str__(self) -> str:
        return self.__repr__()

    def __getitem__(self, unit: int) -> Base:
        return self.bases[unit]

    def __len__(self) -> int:
        return len(self.bases)


@dataclass
class Enzyme:
    amino_acids: List[AminoAcid]

    def iter_amino_acids(self) -> Iterator[AminoAcid]:
        yield from self.amino_acids

    @classmethod
    def from_str(cls, enzyme_str: str) -> "Enzyme":
        amino_acids = []
        for amino_acid_str in enzyme_str.split("-"):
            amino_acid = AminoAcid(amino_acid_str)
            amino_acids.append(amino_acid)
        return cls(amino_acids)

    def __repr__(self) -> str:
        return "-".join([str(b) for b in self.amino_acids])

    def __str__(self) -> str:
        return self.__repr__()

    def __getitem__(self, index: int) -> AminoAcid:
        return self.amino_acids[index]

    def __len__(self) -> int:
        return len(self.amino_acids)


# pylint: disable=too-few-public-methods
class Translator:
    """
    |     | A   | C   | G   | T   |
    | --- | --- | --- | --- | --- |
    | A   |     | cut | del | swi |
    | C   | mvr | mvl | cop | off |
    | G   | ina | inc | ing | int |
    | T   | rpy | rpu | lpy | lpu |
    """

    @classmethod
    def translate(cls, strand: Strand) -> List[Enzyme]:
        enzymes = []
        amino_acids: List[AminoAcid] = []
        for duplet in strand.iter_duplets():
            amino_acid = cls._translate_duplet(duplet)
            if amino_acid is None and len(amino_acids) > 0:
                enzyme = Enzyme(amino_acids)
                enzymes.append(enzyme)
                amino_acids = []
            elif amino_acid is not None:
                amino_acids.append(amino_acid)

        if len(amino_acids) > 0:
            enzyme = Enzyme(amino_acids)
            enzymes.append(enzyme)

        return enzymes

    @classmethod
    def _translate_duplet(cls, duplet: Duplet) -> Optional[AminoAcid]:
        return {
            (Base.A, Base.A): None,
            (Base.A, Base.C): AminoAcid.CUT,
            (Base.A, Base.G): AminoAcid.DEL,
            (Base.A, Base.T): AminoAcid.SWI,
            (Base.C, Base.A): AminoAcid.MVR,
            (Base.C, Base.C): AminoAcid.MVL,
            (Base.C, Base.G): AminoAcid.COP,
            (Base.C, Base.T): AminoAcid.OFF,
            (Base.G, Base.A): AminoAcid.INA,
            (Base.G, Base.C): AminoAcid.INC,
            (Base.G, Base.G): AminoAcid.ING,
            (Base.G, Base.T): AminoAcid.INT,
            (Base.T, Base.A): AminoAcid.RPY,
            (Base.T, Base.C): AminoAcid.RPU,
            (Base.T, Base.G): AminoAcid.LPY,
            (Base.T, Base.T): AminoAcid.LPU,
        }[duplet]


class Orientation(StrEnum):
    U = auto()
    D = auto()
    L = auto()
    R = auto()

    @classmethod
    def from_turning_number(cls, turning_number: int) -> "Orientation":
        return {
            0: cls.R,
            1: cls.D,
            2: cls.L,
            3: cls.U,
        }[turning_number % 4]


class Turn(StrEnum):
    L = auto()
    S = auto()
    R = auto()

    def to_int(self) -> int:
        return {
            Turn.L: -1,
            Turn.S: 0,
            Turn.R: 1,
        }[self]


class Folder:
    """
    | ins | dir |
    | --- | --- |
    | cut | s   |
    | del | s   |
    | swi | r   |
    | mvr | s   |
    | mvl | s   |
    | cop | r   |
    | off | l   |
    | ina | s   |
    | inc | r   |
    | ing | r   |
    | int | l   |
    | rpy | r   |
    | rpu | l   |
    | lpy | l   |
    | lpu | l   |

    | first | last | base |
    | ----- | ---- | ---- |
    | R     | R    | A    |
    | R     | U    | C    |
    | R     | D    | G    |
    | R     | L    | T    |
    """

    @classmethod
    def fold(cls, enzyme: Enzyme) -> Orientation:
        turning_number = 0
        for amino_acid in enzyme.iter_amino_acids():
            turn = cls.get_turn(amino_acid)
            turning_number += turn.to_int()
        return Orientation.from_turning_number(turning_number)

    @classmethod
    def get_turn(cls, amino_acid: AminoAcid) -> Turn:
        return {
            AminoAcid.CUT: Turn.S,
            AminoAcid.DEL: Turn.S,
            AminoAcid.SWI: Turn.R,
            AminoAcid.MVR: Turn.S,
            AminoAcid.MVL: Turn.S,
            AminoAcid.COP: Turn.R,
            AminoAcid.OFF: Turn.L,
            AminoAcid.INA: Turn.S,
            AminoAcid.INC: Turn.R,
            AminoAcid.ING: Turn.R,
            AminoAcid.INT: Turn.L,
            AminoAcid.RPY: Turn.R,
            AminoAcid.RPU: Turn.L,
            AminoAcid.LPY: Turn.L,
            AminoAcid.LPU: Turn.L,
        }[amino_acid]

    @classmethod
    def get_binding_site(cls, enzyme: Enzyme, strand: Strand) -> Optional[int]:
        orientation = cls.fold(enzyme)
        binding_affinity = cls.get_binding_affinity(orientation)
        for unit, base in enumerate(strand.iter_bases()):
            if base == binding_affinity:
                return unit
        return None

    @classmethod
    def get_binding_affinity(cls, orientation: Orientation) -> Base:
        return {
            Orientation.R: Base.A,
            Orientation.U: Base.C,
            Orientation.D: Base.G,
            Orientation.L: Base.T,
        }[orientation]


@dataclass
class BasePair:
    bind: Optional[Base]
    comp: Optional[Base]

    def swap(self) -> None:
        self.bind, self.comp = self.comp, self.bind

    def add_comp(self) -> None:
        assert self.bind is not None
        self.comp = self.bind.get_complement()


class Rewriter:
    """
    | ins | action                                         |
    | --- | ---------------------------------------------- |
    | cut | cut strand(s)                                  |
    | del | delete a base from strand                      |
    | swi | switch enzyme to other strand                  |
    | mvr | move one unit to the right                     |
    | mvl | move one unit to the left                      |
    | cop | turn on Copy mode                              |
    | off | turn off Copy mode                             |
    | ina | insert A to the right of this unit             |
    | inc | insert C to the right of this unit             |
    | ing | insert G to the right of this unit             |
    | int | insert T to the right of this unit             |
    | rpy | search for the nearest pyrimidine to the right |
    | rpu | search for the nearest purine to the right     |
    | lpy | search for the nearest pyrimidine to the left  |
    | lpu | search for the nearest purine to the left      |
    """

    # pylint: disable=too-many-branches
    @classmethod
    def rewrite(cls, enzyme: Enzyme, strand: Strand) -> List[Strand]:
        copy_mode = False

        unit = Folder.get_binding_site(enzyme, strand)
        logger.info(f"Rewriting strand {strand} with enzyme {enzyme}, unit={unit}")
        if unit is None:
            return [strand]

        pairs = [BasePair(base, None) for base in strand.iter_bases()]

        logger.debug(f"Init @ {unit}, copy={copy_mode}")
        logger.debug(cls.pairs_to_string(pairs))

        strands = []
        for amino_acid in enzyme.iter_amino_acids():
            logger.debug(f"Applying {amino_acid} @ {unit}, copy={copy_mode}")

            if amino_acid == AminoAcid.CUT:
                # TODO: Handle case where we cut at the very end of the strand
                cut_pairs = pairs[unit + 1 :]
                strands += cls.strands_from_pairs(cut_pairs)
                pairs = pairs[: unit + 1]
            elif amino_acid == AminoAcid.DEL:
                pairs[unit].bind = None
                unit -= 1
                # NOTE: It's not clear from the specification which direction we should move
                # after a deletion, we here we choose left arbitrarily.
                if unit < 0:
                    logger.debug("Reached end of strand")
                    break
            elif amino_acid == AminoAcid.SWI:
                if pairs[unit].comp is None:
                    logger.debug("Tried to switch to empty base pair complement")
                    break
                for pair in pairs:
                    pair.swap()
                pairs = pairs[::-1]
                unit = len(pairs) - unit - 1
            elif amino_acid in [AminoAcid.MVR, AminoAcid.MVL]:
                unit += cls.amino_acid_to_direction(amino_acid)
                if unit < 0 or unit >= len(pairs):
                    logger.debug("Reached end of strand")
                    break
                if pairs[unit].bind is None:
                    logger.debug("Reached end of strand")
                    break
                if copy_mode:
                    pairs[unit].add_comp()
            elif amino_acid == AminoAcid.COP:
                copy_mode = True
                pair = pairs[unit]
                assert pair.bind is not None
                pair.comp = pair.bind.get_complement()
            elif amino_acid == AminoAcid.OFF:
                copy_mode = False
            elif amino_acid in [AminoAcid.INA, AminoAcid.INC, AminoAcid.ING, AminoAcid.INT]:
                bind = cls.amino_acid_to_base(amino_acid)
                comp = bind.get_complement() if copy_mode else None
                pairs.insert(unit + 1, BasePair(bind, comp))
            elif amino_acid in [AminoAcid.RPY, AminoAcid.RPU, AminoAcid.LPY, AminoAcid.LPU]:
                end_of_strand = False
                while True:
                    unit += cls.amino_acid_to_direction(amino_acid)
                    if unit < 0 or unit >= len(pairs):
                        end_of_strand = True
                        break
                    pair = pairs[unit]
                    bind_base = pair.bind
                    if bind_base is None:
                        end_of_strand = True
                        break
                    if copy_mode:
                        pair.add_comp()
                    if bind_base.is_type(cls.amino_acid_to_base_type(amino_acid)):
                        break
                if end_of_strand:
                    logger.debug("Reached end of strand")
                    break

            logger.debug(cls.pairs_to_string(pairs))

        strands += cls.strands_from_pairs(pairs)
        return strands

    @classmethod
    def strands_from_pairs(cls, pairs: List[BasePair]) -> List[Strand]:
        strands = []
        bind_bases = []
        comp_bases = []

        for pair in pairs:
            if pair.bind is not None:
                bind_bases.append(pair.bind)
            elif len(bind_bases) > 0:
                strands.append(Strand(bind_bases))
                bind_bases = []

            if pair.comp is not None:
                comp_bases.append(pair.comp)
            elif len(comp_bases) > 0:
                strands.append(Strand(comp_bases[::-1]))
                comp_bases = []

        if len(bind_bases) > 0:
            strands.append(Strand(bind_bases))
        if len(comp_bases) > 0:
            strands.append(Strand(comp_bases[::-1]))

        return strands

    @classmethod
    def amino_acid_to_base(cls, amino_acid: AminoAcid) -> Base:
        return {
            AminoAcid.INA: Base.A,
            AminoAcid.INC: Base.C,
            AminoAcid.ING: Base.G,
            AminoAcid.INT: Base.T,
        }[amino_acid]

    @classmethod
    def amino_acid_to_base_type(cls, amino_acid: AminoAcid) -> BaseType:
        return {
            AminoAcid.RPY: BaseType.PYRIMIDINE,
            AminoAcid.RPU: BaseType.PURINE,
            AminoAcid.LPY: BaseType.PYRIMIDINE,
            AminoAcid.LPU: BaseType.PURINE,
        }[amino_acid]

    @classmethod
    def amino_acid_to_direction(cls, amino_acid: AminoAcid) -> int:
        return {
            AminoAcid.RPY: 1,
            AminoAcid.RPU: 1,
            AminoAcid.LPY: -1,
            AminoAcid.LPU: -1,
            AminoAcid.MVR: 1,
            AminoAcid.MVL: -1,
        }[amino_acid]

    @classmethod
    def pairs_to_string(cls, pairs: List[BasePair]) -> str:
        res = "[ "
        comp_map = {Base.A: "∀", Base.C: "Ↄ", Base.G: "⅁", Base.T: "⊥"}
        for pair in pairs:
            if pair.comp is None:
                res += "  "
            else:
                res += str(comp_map[pair.comp]) + " "
        res += "]\n[ "
        for pair in pairs:
            if pair.bind is None:
                res += "  "
            else:
                res += str(pair.bind) + " "
        res += "]"
        return res
