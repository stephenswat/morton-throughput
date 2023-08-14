import csv
import json
import pathlib
import re

import yaml

with open(snakemake.output[0], "w") as o:
    writer = csv.DictWriter(
        o, fieldnames=["file", "mca", "uarch", "compiler", "layout", "dims", "rthr"]
    )
    writer.writeheader()

    for f in snakemake.input:
        path = pathlib.PurePath(f)

        if path.suffix == ".yaml":
            if m := re.fullmatch(
                r"(?P<file>[a-zA-Z0-9_]+)\.osaca-(?P<uarch>[a-zA-Z0-9_]+)-(?P<compiler>[a-zA-Z0-9_]+)-(?P<layout>[a-zA-Z0-9_]+)-(?P<dims>\d+)",
                path.stem,
            ):
                with open(path, "r") as f:
                    data = yaml.safe_load(f)

                    writer.writerow(
                        {
                            "file": m.group("file"),
                            "mca": "osaca",
                            "uarch": m.group("uarch"),
                            "compiler": m.group("compiler"),
                            "layout": m.group("layout"),
                            "dims": int(m.group("dims")),
                            "rthr": max(data["Summary"]["PortPressure"].values()),
                        }
                    )
            else:
                raise Exception("YAML file but not an OSACA file!")
        elif path.suffix == ".json":
            if m := re.fullmatch(
                r"(?P<file>[a-zA-Z0-9_]+)\.llvm-mca-(?P<uarch>[a-zA-Z0-9_]+)-(?P<compiler>[a-zA-Z0-9_]+)",
                path.stem,
            ):
                with open(path, "r") as f:
                    data = json.load(f)

                    for r in data["CodeRegions"]:
                        if n := re.fullmatch(
                            r"(?P<layout>[a-zA-Z0-9_]+)\$(?P<dims>\d+)", r["Name"]
                        ):
                            writer.writerow(
                                {
                                    "file": m.group("file"),
                                    "mca": "llvm-mca",
                                    "uarch": m.group("uarch"),
                                    "compiler": m.group("compiler"),
                                    "layout": n.group("layout"),
                                    "dims": int(n.group("dims")),
                                    "rthr": r["SummaryView"]["TotalCycles"]
                                    / r["SummaryView"]["Iterations"],
                                }
                            )
                        else:
                            raise Exception("Invalid region name!")
            else:
                raise Exception("JSON file but not an LLVM-MCA file!")
        else:
            raise Exception("Invalid file extension!")
