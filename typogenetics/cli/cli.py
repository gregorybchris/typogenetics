import logging
import random
from typing import Optional

import click

from typogenetics.lib.typogenetics import Enzyme, Rewriter, Strand, Translator

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

    enzymes = Translator.translate(Strand.from_str(strand_str))
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
@click.argument("init_strand", type=str)
@click.option("--iter", "n_iterations", type=int, default=100_000)
@click.option("--seed", "random_seed", type=int, default=None)
@click.option("--debug", type=bool, is_flag=True)
@click.option("--print-strands", type=bool, is_flag=True)
def simulate_command(
    init_strand: str,
    n_iterations: int,
    random_seed: Optional[int] = None,
    debug: bool = False,
    print_strands: bool = False,
) -> None:
    set_logging_config(debug)

    random.seed(random_seed)
    strands = [Strand.from_str(init_strand)]
    for _ in range(n_iterations):
        enzyme_strand = strands[random.randint(0, len(strands) - 1)]
        enzymes = Translator.translate(enzyme_strand)
        if len(enzymes) == 0:
            continue
        enzyme = enzymes[random.randint(0, len(enzymes) - 1)]
        rewrite_strand = strands[random.randint(0, len(strands) - 1)]
        new_strands = Rewriter.rewrite(enzyme, rewrite_strand)
        strands.extend(new_strands)

    unique_strands = set()
    for strand in strands:
        unique_strands.add(str(strand))

    if print_strands:
        print("Unique strands:")
        sorted_strands = sorted(unique_strands)
        for strand_str in sorted_strands:
            print(f"- {strand_str}")

    print(f"Discovered {len(unique_strands)} unique strands while simulating for {n_iterations} iterations")
