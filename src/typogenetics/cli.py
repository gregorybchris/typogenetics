import logging
from typing import Annotated, Optional

import numpy as np
from rich.console import Console
from rich.logging import RichHandler
from typer import Argument, Option, Typer

from typogenetics.search import Search
from typogenetics.typogenetics import Enzyme, Folder, Rewriter, Strand, Translator

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
        console.print(enzyme)


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
        console.print(f"- {new_strand}")


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
    console.print(strand)

    enzymes = Translator.translate(strand)
    console.print(enzymes)

    enzyme = enzymes[0]
    unit = Folder.get_binding_site(enzyme, strand)
    console.print(unit)

    new_strands = Rewriter.rewrite(enzyme, strand)
    console.print(new_strands)
