import numpy as np

from typogenetics.search import Editor, EditType
from typogenetics.typogenetics import Strand


class TestSearch:
    def test_select_edit_type(self) -> None:
        rng = np.random.default_rng(42)
        assert Editor.select_edit_type(rng) == EditType.MUTATE

    def test_mutate(self) -> None:
        rng = np.random.default_rng(42)
        strand = Strand.from_str("ACGT")
        new_strand = Editor.mutate(strand, rng)
        assert new_strand == Strand.from_str("TCGT")

    def test_insert(self) -> None:
        rng = np.random.default_rng(42)
        strand = Strand.from_str("ACGT")
        new_strand = Editor.insert(strand, rng)
        assert new_strand == Strand.from_str("TACGT")

    def test_delete(self) -> None:
        rng = np.random.default_rng(42)
        strand = Strand.from_str("ACGT")
        new_strand = Editor.delete(strand, rng)
        assert new_strand == Strand.from_str("CGT")
