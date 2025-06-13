import logging
from typing import Optional

import click
import numpy as np

from typogenetics.search import Search
from typogenetics.typogenetics import Enzyme, Folder, Rewriter, Strand, Translator

logger = logging.getLogger(__name__)


def set_logging_config(debug: bool) -> None:
    log_level = logging.INFO
    if debug:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level, format="%(message)s")


@click.group()
def main() -> None:
    """Main CLI entrypoint."""


@main.command(name="translate")
@click.argument("strand-str", type=str)
@click.option("--debug", type=bool, is_flag=True)
def translate_command(
    strand_str: str,
    debug: bool = False,
) -> None:
    set_logging_config(debug)

    strand = Strand.from_str(strand_str)
    enzymes = Translator.translate(strand)
    for enzyme in enzymes:
        print(enzyme)


@main.command(name="rewrite")
@click.argument("enzyme-str", type=str)
@click.argument("strand-str", type=str)
@click.option("--debug", type=bool, is_flag=True)
def rewrite_command(
    enzyme_str: str,
    strand_str: str,
    debug: bool = False,
) -> None:
    set_logging_config(debug)

    enzyme = Enzyme.from_str(enzyme_str)
    strand = Strand.from_str(strand_str)
    new_strands = Rewriter.rewrite(enzyme, strand)
    print("New strands:")
    for new_strand in new_strands:
        print(f"- {new_strand}")


@main.command(name="simulate")
@click.argument("init-strand-str", type=str)
@click.option("--iter", "n_iterations", type=int, default=100_000)
@click.option("--seed", "random_state", type=int, default=None)
@click.option("--debug", type=bool, is_flag=True)
@click.option("--print-strands", type=bool, is_flag=True)
def simulate_command(
    init_strand_str: str,
    n_iterations: int,
    random_state: Optional[int] = None,
    debug: bool = False,
    print_strands: bool = False,
) -> None:
    set_logging_config(debug)

    rng = np.random.default_rng(random_state)
    init_strand = Strand.from_str(init_strand_str)

    Search.random(init_strand, n_iterations, rng, print_strands=print_strands)


@main.command(name="search")
@click.argument("init-strand-str", type=str)
@click.argument("apply-strand-str", type=str)
@click.option("--depth", "target_depth", type=int, default=10)
@click.option("--edits", "n_edits", type=int, default=10)
@click.option("--seed", "random_state", type=int, default=None)
@click.option("--debug", type=bool, is_flag=True)
@click.option("--print-strands", type=bool, is_flag=True)
def search_command(
    init_strand_str: str,
    apply_strand_str: str,
    target_depth: int,
    n_edits: int,
    random_state: Optional[int] = None,
    debug: bool = False,
    print_strands: bool = False,
) -> None:
    set_logging_config(debug)

    rng = np.random.default_rng(random_state)
    init_strand = Strand.from_str(init_strand_str)
    apply_strand = Strand.from_str(apply_strand_str)

    Search.bfs(init_strand, apply_strand, target_depth, n_edits, rng, print_strands=print_strands)


@main.command(name="go")
@click.option("--debug", type=bool, is_flag=True)
def go_command(
    debug: bool = False,
) -> None:
    set_logging_config(debug)

    strand = Strand.from_str("CAAGGGTATACCCCATATCCT")
    print(strand)

    enzymes = Translator.translate(strand)
    print(enzymes)

    enzyme = enzymes[0]
    unit = Folder.get_binding_site(enzyme, strand)
    print(unit)

    new_strands = Rewriter.rewrite(enzyme, strand)
    print(new_strands)
