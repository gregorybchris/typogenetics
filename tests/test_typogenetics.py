import pytest

from typogenetics.lib.main import AminoAcid, Base, Binder, Enzyme, Folder, Orientation, Rewriter, Strand, Translator


class TestTypogenetics:
    def test_strand_from_str(self) -> None:
        assert Strand.from_str("CG GA TA CT AA AC CG A") == Strand(
            [
                Base.C,
                Base.G,
                Base.G,
                Base.A,
                Base.T,
                Base.A,
                Base.C,
                Base.T,
                Base.A,
                Base.A,
                Base.A,
                Base.C,
                Base.C,
                Base.G,
                Base.A,
            ]
        )

    def test_translator(self) -> None:
        strand = Strand.from_str("CG GA TA CT AA AC CG A")
        assert Translator.translate(strand) == [
            Enzyme([AminoAcid.COP, AminoAcid.INA, AminoAcid.RPY, AminoAcid.OFF]),
            Enzyme([AminoAcid.CUT, AminoAcid.COP]),
        ]

    def test_folder(self) -> None:
        enzyme = Enzyme([AminoAcid.COP, AminoAcid.INA, AminoAcid.RPY, AminoAcid.OFF])
        assert Folder.fold(enzyme) == Orientation.D

    def test_binder(self) -> None:
        enzyme = Enzyme([AminoAcid.COP, AminoAcid.INA, AminoAcid.RPY, AminoAcid.OFF])
        strand = Strand.from_str("CG GA TA CT AA AC CG A")
        orientation = Folder.fold(enzyme)
        assert orientation == Orientation.D
        binding_affinity = Binder.get_binding_affinity(orientation)
        assert Binder.get_binding_site(binding_affinity, strand) == 1

    @pytest.mark.xfail(reason="Instructions not implemented")
    def test_rewriter(self) -> None:
        enzyme = Enzyme([AminoAcid.COP, AminoAcid.INA, AminoAcid.RPY, AminoAcid.OFF])
        strand = Strand.from_str("CG GA TA CT AA AC CG A")
        new_strands = Rewriter.rewrite(enzyme, strand)
        assert new_strands == [
            Strand.from_str("CG"),
        ]
