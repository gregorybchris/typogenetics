from typogenetics.lib.main import AminoAcid, Base, Enzyme, Folder, Orientation, Strand, Translator


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
        assert Folder.fold(enzyme.amino_acids) == Orientation.D
