from datetime import datetime
from typing import Literal

import click
from gcn_kafka import Consumer
from rich.console import Console
from rich.prompt import Prompt

from doppelganger.notices import fermigbm as gbm
from doppelganger.notices import parse_notice
from doppelganger.types import CelestialCoords

console = Console()


TOPICS_PRESETS = {
    "fermi-gbm": [
        'gcn.classic.voevent.FERMI_GBM_ALERT',
        'gcn.classic.voevent.FERMI_GBM_FIN_POS',
        'gcn.classic.voevent.FERMI_GBM_FLT_POS',
        'gcn.classic.voevent.FERMI_GBM_GND_POS',
    ],
}


def info_from_notice(content: str):
    d = parse_notice(content)
    return (
        gbm.get_issue_time(d),
        gbm.get_obs_time(d),
        gbm.get_localization_coords(d),
        gbm.get_localization_error(d),
    )


def formatted(issuet: datetime, obst: datetime, loc: CelestialCoords, err: float) -> str:
    return f"""
    \tIssued: {issuet.strftime("%m/%d/%Y, %H:%M:%S")}
    \tTime:   {obst.strftime("%m/%d/%Y, %H:%M:%S")}
    \tCoords: {str(loc)}
    \tError:  {err:.2f} deg,
    """


@click.command()
@click.argument(
    "instruments",
    type=click.Choice(["fermi-gbm"]),
    nargs=-1,
)
@click.option(
    "--userid",
    "-u",
    type=str,
    help="A NASA GCN Kafka userid.",
)
@click.option(
    "--secret",
    "-s",
    type=str,
    help="A NASA GCN Kafka secret.",
)
@click.option(
    "--offset",
    type=click.Choice(["earliest", "latest"]),
    default="earliest",
    help="Sets the points in time from when the stream will start. "
    "Setting to `earliest` will stream all past events reachable, before start streaming new ones."
)
def main(
    instruments: list[Literal["fermi-gbm"]],
    userid: str,
    secret: str,
    offset: Literal["earliest", "latest"]
):
    if not userid:
        userid = Prompt.ask("Enter NASA-GCN userid")
    if not secret:
        secret = Prompt.ask("Enter NASA-GCN secret")

    console.print("Creating consumer.")
    consumer = Consumer(
        client_id=userid,
        client_secret=secret,
        config={
            "auto.offset.reset": offset,
        },
    )
    for instrument in instruments:
        console.print(f"Subscribing to {instrument} topics.")
        consumer.subscribe(TOPICS_PRESETS[instrument])

    console.print("Starting streaming events.")
    while True:
        for message in consumer.consume(timeout=1):
            if message.error():
                console.print(message.error())
                continue
            topic = message.topic().split('.')[-1].lower()
            offset_str = str(message.offset())
            content = message.value().decode()

            console.print("------------")
            console.print(f"Received VOEvent from topic={topic}, at offset={offset_str}.")
            console.print(formatted(*info_from_notice(content)), end="\n\n")


if __name__ == "__main__":
    main()