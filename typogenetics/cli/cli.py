import logging

import click

from typogenetics.lib.main import Strand, Translator

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
    print("n_translations", n_translations)

    strand = Strand.from_str("CG GA TA CT AA AC CG A")
    # [cop, ina, rpy, off] and [cut, cop]
    print(strand)
    enzymes = Translator.translate(strand)
    for enzyme in enzymes:
        print(enzyme)
