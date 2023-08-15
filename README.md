# Morton Layout Throughput Analysis

As an alternative to the canonical row-major and column-major array layouts, we
can use the Morton curve layout to improve cache locality in some applications.
Morton layours rely on bit-scattering operations and bit-wise ORs rather than
addition and multiplication; how fast are these operations? Can we calculate
Morton indices with similar performance to canonical layouts?

This repository contains a microcode analysis based on LLVM-MCA and OSACA which
calculates the throughput of index calculations for canonical and Morton
layouts for _gcc_ and _clang_ across a range of dimensionalities.

## Requirements

This project requires Poetry to be installed. It can be installed through your
favourite Python package manager. It also requires llvm-mca, clang++, and g++
to be available.

## Usage

To use this repository, perform a one-time installation of its dependencies:

```bash
poetry install
```

Then run the analysis through Snakemake:

```bash
poetry run snakemake -c
```

You will then find the plots in `out/plots/`.
