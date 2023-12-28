import logging

import click

from typogenetics.lib.typogenetics import Rewriter, Strand, Translator

logger = logging.getLogger(__name__)


@click.group()
def main() -> None:
    """Main CLI entrypoint."""


@main.command(name="run")
@click.option("--n-translations", type=int, default=10)
@click.option("--debug", type=bool, is_flag=True)
def run_command(
    n_translations: int,
    debug: bool,
) -> None:
    log_level = logging.INFO
    if debug:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level, format="%(message)s")

    print("n_translations", n_translations)

    strand = Strand.from_str("CG GA TA CT AA AC CG A")
    enzymes = Translator.translate(strand)
    # [cop, ina, rpy, off] and [cut, cop]
    enzyme = enzymes[0]
    logger.debug(f"Rewriting strand {strand} with enzyme {enzyme}")
    Rewriter.rewrite(enzyme, strand)
