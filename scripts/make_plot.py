import pathlib
import re

import matplotlib.pyplot
import numpy
import pandas

matplotlib.rcParams.update(
    {
        "font.size": 8,
        "font.family": "serif",
        "text.usetex": True,
        "text.latex.preamble": """
        \\usepackage{libertine}
        \\usepackage[libertine]{newtxmath}
        """,
    }
)

df = pandas.read_csv(snakemake.input[0]).sort_values(by=["dims"])

compilers = ["clang", "gcc"]
layouts = ["Canonical", "Morton"]

for o in snakemake.output:
    path = pathlib.PurePath(o)

    if m := re.fullmatch(
        r"(?P<file>[a-zA-Z0-9_]+).(?P<mca>[a-zA-Z0-9_\-]+)", path.stem
    ):
        fig, ax = matplotlib.pyplot.subplots(
            figsize=(3.333, 1.4), constrained_layout=True
        )

        ax.set_xlabel("Number of dimensions")
        ax.set_ylabel("Cycles per iteration")

        ax.grid(True, axis="y")

        for layout in layouts:
            for compiler in compilers:
                dff = df[
                    (df["compiler"] == compiler)
                    & (df["layout"] == layout)
                    & (df["uarch"] == "haswell")
                    & (df["mca"] == m.group("mca"))
                    & (df["file"] == m.group("file"))
                ]
                layout_name = "Canon." if layout == "Canonical" else layout
                ax.plot(
                    dff["dims"],
                    dff["rthr"],
                    "+-" if compiler == "clang" else "x-",
                    label=f"{layout_name} (\\textsc{{{compiler}}})",
                    linewidth=1,
                )

        ax.xaxis.set_ticks(numpy.arange(1, 11, 1))
        ax.legend(ncol=2, handletextpad=0.2, columnspacing=0.5, handlelength=1)

        fig.savefig(path, bbox_inches="tight", pad_inches=0.02)
    else:
        raise Exception("Invalid file name!")
