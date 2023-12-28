import logging
from dataclasses import dataclass
from enum import StrEnum, auto
from typing import Iterator, List, Optional, Tuple

logger = logging.getLogger(__name__)


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

    def is_purine(self) -> bool:
        return self in [Base.A, Base.G]

    def is_pyrimidine(self) -> bool:
        return self in [Base.C, Base.T]


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
        idx = 0
        while True:
            if idx + 1 >= len(self):
                break

            yield (self[idx], self[idx + 1])

            idx += 2

    def __repr__(self) -> str:
        return "".join([str(b) for b in self.bases])

    def __str__(self) -> str:
        return self.__repr__()

    def __getitem__(self, index: int) -> Base:
        return self.bases[index]

    def __len__(self) -> int:
        return len(self.bases)


@dataclass
class Enzyme:
    amino_acids: List[AminoAcid]

    def __repr__(self) -> str:
        return "-".join([str(b) for b in self.amino_acids])

    def __str__(self) -> str:
        return self.__repr__()

    def __getitem__(self, index: int) -> AminoAcid:
        return self.amino_acids[index]

    def __len__(self) -> int:
        return len(self.amino_acids)


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
            if amino_acid is None:
                enzyme = Enzyme(amino_acids)
                enzymes.append(enzyme)
                amino_acids = []
            else:
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
    """

    @classmethod
    def fold(cls, amino_acids: List[AminoAcid]) -> Orientation:
        turning_number = 0
        for amino_acid in amino_acids:
            turn = cls._get_turn(amino_acid)
            turning_number += turn.to_int()
        return Orientation.from_turning_number(turning_number)

    @classmethod
    def _get_turn(cls, amino_acid: AminoAcid) -> Turn:
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
