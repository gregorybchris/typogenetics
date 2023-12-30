import logging
import random
from typing import Optional

import click

from typogenetics.typogenetics import Enzyme, Rewriter, Strand, Translator

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
@click.option("--seed", "random_seed", type=int, default=None)
@click.option("--debug", type=bool, is_flag=True)
@click.option("--print-strands", type=bool, is_flag=True)
def simulate_command(
    init_strand_str: str,
    n_iterations: int,
    random_seed: Optional[int] = None,
    debug: bool = False,
    print_strands: bool = False,
) -> None:
    set_logging_config(debug)

    random.seed(random_seed)
    init_strand = Strand.from_str(init_strand_str)
    strands = [init_strand]
    known_set = set([str(init_strand)])
    for _ in range(n_iterations):
        enzyme_strand = strands[random.randint(0, len(strands) - 1)]
        enzymes = Translator.translate(enzyme_strand)
        if len(enzymes) == 0:
            continue
        enzyme = enzymes[random.randint(0, len(enzymes) - 1)]
        rewrite_strand = strands[random.randint(0, len(strands) - 1)]
        new_strands = Rewriter.rewrite(enzyme, rewrite_strand)
        for strand in new_strands:
            if str(strand) not in known_set:
                strands.append(strand)
            known_set.add(str(strand))

    if print_strands:
        print("Unique strands:")
        sorted_strands = sorted(known_set)
        for strand_str in sorted_strands:
            print(f"- {strand_str}")

    print(f"Discovered {len(known_set)} unique strands while simulating for {n_iterations} iterations")
