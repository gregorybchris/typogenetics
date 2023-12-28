import logging

import click

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
