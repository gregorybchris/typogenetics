from typogenetics.lib.typogenetics import Base, Enzyme, Folder, Orientation, Rewriter, Strand, Translator


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
            Enzyme.from_str("cop-ina-rpy-off"),
            Enzyme.from_str("cut-cop"),
        ]

    def test_folder(self) -> None:
        enzyme = Enzyme.from_str("cop-ina-rpy-off")
        strand = Strand.from_str("CG GA TA CT AA AC CG A")
        assert Folder.fold(enzyme) == Orientation.D
        assert Folder.get_binding_site(enzyme, strand) == 1

    def test_rewriter(self) -> None:
        enzyme = Enzyme.from_str("cop-ina-rpy-off")
        strand = Strand.from_str("CGGATACTAAACCGA")
        new_strands = Rewriter.rewrite(enzyme, strand)
        assert new_strands == [
            Strand.from_str("ATCTC"),
            Strand.from_str("CGAGATACTAAACCGA"),
        ]
