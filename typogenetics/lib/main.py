import logging
from dataclasses import dataclass
from enum import StrEnum, auto
from typing import Iterator, List, Optional, Tuple, Type

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
    def fold(cls, enzyme: Enzyme) -> Orientation:
        turning_number = 0
        for amino_acid in enzyme.iter_amino_acids():
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


@dataclass
class RewriterState:
    """Result of instruction application."""

    strand: Strand
    # TODO: Update
    complement: Strand
    unit: int

    def get_strands(self) -> List[Strand]:
        # TODO: Update
        return [self.strand, self.complement]


@dataclass
class Inst:
    """Base class for instructions."""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Cut(Inst):
    """Cut strand(s)"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Del(Inst):
    """Delete a base from strand"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Swi(Inst):
    """Switch enzyme to other strand"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Mvr(Inst):
    """Move one unit to the right"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Mvl(Inst):
    """Move one unit to the left"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Cop(Inst):
    """Turn on Copy mode"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Off(Inst):
    """Turn off Copy mode"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Ina(Inst):
    """Insert A to the right of this unit"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Inc(Inst):
    """Insert C to the right of this unit"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Ing(Inst):
    """Insert G to the right of this unit"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Int(Inst):
    """Insert T to the right of this unit"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Rpy(Inst):
    """Search for the nearest pyrimidine to the right"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Rpu(Inst):
    """Search for the nearest purine to the right"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Lpy(Inst):
    """Search for the nearest pyrimidine to the left"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Lpu(Inst):
    """Search for the nearest purine to the left"""

    def __call__(self, state: RewriterState) -> RewriterState:
        raise NotImplementedError()


class Binder:
    @classmethod
    def get_binding_site(cls, binding_affinity: Base, strand: Strand) -> Optional[int]:
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


class Rewriter:
    @classmethod
    def rewrite(cls, enzyme: Enzyme, strand: Strand) -> List[Strand]:
        orientation = Folder.fold(enzyme)
        binding_affinity = Binder.get_binding_affinity(orientation)
        unit = Binder.get_binding_site(binding_affinity, strand)
        if unit is None:
            return []

        # TODO: Fix this to not have two of the same strand
        state = RewriterState(strand, strand, unit)
        for amino_acid in enzyme.iter_amino_acids():
            inst_class = cls._get_instruction_class(amino_acid)
            inst = inst_class()
            state = inst(state)
        return state.get_strands()

    @classmethod
    def _get_instruction_class(cls, amino_acid: AminoAcid) -> Type[Inst]:
        return {
            AminoAcid.CUT: Cut,
            AminoAcid.DEL: Del,
            AminoAcid.SWI: Swi,
            AminoAcid.MVR: Mvr,
            AminoAcid.MVL: Mvl,
            AminoAcid.COP: Cop,
            AminoAcid.OFF: Off,
            AminoAcid.INA: Ina,
            AminoAcid.INC: Inc,
            AminoAcid.ING: Ing,
            AminoAcid.INT: Int,
            AminoAcid.RPY: Rpy,
            AminoAcid.RPU: Rpu,
            AminoAcid.LPY: Lpy,
            AminoAcid.LPU: Lpu,
        }[amino_acid]
