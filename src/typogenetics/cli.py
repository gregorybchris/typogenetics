import logging
from typing import Annotated, Optional

import numpy as np
from rich.console import Console
from rich.logging import RichHandler
from typer import Argument, Option, Typer

from typogenetics.search import Search
from typogenetics.typogenetics import AminoAcid, Base, Enzyme, Folder, Rewriter, Strand, Translator

logger = logging.getLogger(__name__)

app = Typer(pretty_exceptions_enable=False)
console = Console()


def set_logger_config(info: bool, debug: bool) -> None:
    handlers = [RichHandler(rich_tracebacks=True)]
    log_format = "%(message)s"

    if info:
        logging.basicConfig(level=logging.INFO, handlers=handlers, format=log_format)
    if debug:
        logging.basicConfig(level=logging.DEBUG, handlers=handlers, format=log_format)


def amino_acid_to_console(amino_acid: AminoAcid) -> str:
    color = "white"
    if amino_acid in [AminoAcid.LPU, AminoAcid.LPY, AminoAcid.RPU, AminoAcid.RPY]:
        color = "blue"
    elif amino_acid in [AminoAcid.CUT, AminoAcid.DEL, AminoAcid.SWI]:
        color = "red"
    elif amino_acid in [AminoAcid.MVR, AminoAcid.MVL]:
        color = "yellow"
    elif amino_acid in [AminoAcid.INA, AminoAcid.INC, AminoAcid.ING, AminoAcid.INT]:
        color = "green"
    elif amino_acid in [AminoAcid.COP]:
        color = "magenta"
    elif amino_acid in [AminoAcid.OFF]:
        color = "orange"
    return f"[{color}]{amino_acid.value}[/]"


def enzyme_to_console(enzyme: Enzyme) -> str:
    return "[white]-[/]".join(amino_acid_to_console(amino_acid) for amino_acid in enzyme.iter_amino_acids())


def base_to_console(base: Base) -> str:
    match base:
        case Base.C:
            return "[green]C[/]"
        case Base.G:
            return "[blue]G[/]"
        case Base.T:
            return "[yellow]T[/]"
        case Base.A:
            return "[red]A[/]"


def strand_to_console(strand: Strand) -> str:
    return "".join(base_to_console(base) for base in strand.iter_bases())


@app.command()
def translate(
    strand_str: Annotated[str, Argument(...)],
    info: Annotated[bool, Option("--info/--no-info")] = True,
    debug: Annotated[bool, Option("--debug/--no-debug")] = False,
) -> None:
    set_logger_config(info, debug)

    strand = Strand.from_str(strand_str)
    enzymes = Translator.translate(strand)
    for enzyme in enzymes:
        console.print(enzyme_to_console(enzyme))


@app.command()
def rewrite(
    enzyme_str: Annotated[str, Argument(...)],
    strand_str: Annotated[str, Argument(...)],
    info: Annotated[bool, Option("--info/--no-info")] = True,
    debug: Annotated[bool, Option("--debug/--no-debug")] = False,
) -> None:
    set_logger_config(info, debug)

    enzyme = Enzyme.from_str(enzyme_str)
    strand = Strand.from_str(strand_str)
    new_strands = Rewriter.rewrite(enzyme, strand)
    console.print("New strands:")
    for new_strand in new_strands:
        console.print(f"- {strand_to_console(new_strand)}")


@app.command()
def simulate(  # noqa: PLR0913
    init_strand_str: Annotated[str, Argument(...)],
    n_iterations: Annotated[int, Option("--iter")] = 100_000,
    seed: Annotated[Optional[int], Option("--seed")] = None,
    print_strands: Annotated[bool, Option("--print-strands/--no-print-strands")] = False,
    info: Annotated[bool, Option("--info/--no-info")] = True,
    debug: Annotated[bool, Option("--debug/--no-debug")] = False,
) -> None:
    set_logger_config(info, debug)

    rng = np.random.default_rng(seed)
    init_strand = Strand.from_str(init_strand_str)

    Search.random(init_strand, n_iterations, rng, print_strands=print_strands)


@app.command()
def search(  # noqa: PLR0913
    init_strand_str: Annotated[str, Argument(...)],
    apply_strand_str: Annotated[str, Argument(...)],
    target_depth: Annotated[int, Option("--depth")] = 10,
    seed: Annotated[Optional[int], Option("--seed")] = None,
    n_edits: Annotated[int, Option("--edits")] = 10,
    print_strands: Annotated[bool, Option("--print-strands/--no-print-strands")] = False,
    info: Annotated[bool, Option("--info/--no-info")] = True,
    debug: Annotated[bool, Option("--debug/--no-debug")] = False,
) -> None:
    set_logger_config(info, debug)

    rng = np.random.default_rng(seed)
    init_strand = Strand.from_str(init_strand_str)
    apply_strand = Strand.from_str(apply_strand_str)

    Search.bfs(init_strand, apply_strand, target_depth, n_edits, rng, print_strands=print_strands)


@app.command()
def go(
    info: Annotated[bool, Option("--info/--no-info")] = True,
    debug: Annotated[bool, Option("--debug/--no-debug")] = False,
) -> None:
    set_logger_config(info, debug)

    strand = Strand.from_str("CAAGGGTATACCCCATATCCT")
    console.print(f"Strand: {strand_to_console(strand)}")

    enzymes = Translator.translate(strand)
    console.print(f"Enzymes: {', '.join(enzyme_to_console(enzyme) for enzyme in enzymes)}")

    enzyme = enzymes[0]
    unit = Folder.get_binding_site(enzyme, strand)
    console.print(f"Binding site: {unit}")

    new_strands = Rewriter.rewrite(enzyme, strand)
    console.print(f"New strands: {', '.join(strand_to_console(new_strand) for new_strand in new_strands)}")
