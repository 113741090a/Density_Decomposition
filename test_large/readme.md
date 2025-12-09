# Density Decomposition

This project provides an efficient algorithm for computing precise density decomposition, along with implementations of various iterative algorithms for hypergraph density decomposition. It includes tools for data collection to evaluate approximation quality, measurement of time and maximum memory usage, and visualization of results. Currently, the implementation is tested on a single small graph.

**For convenience, we provide master scripts to automatically run the experiments and generate results. See the last section for details.**

Implemented iterative algorithms include: Elist++, FISTA*, Greedy++, PR, PR_EXP, and PR_LIN.

The code has been tested on Linux-based systems.

The workflow consists of approximately five components. To generate results, navigate to the `test_small/` directory as your working folder and execute `master_script.sh`. This will produce an `images` folder. Copy all contents from this folder into `test_small/source_latex/sigmod_final/images`. Then, navigate to `test_small/source_latex/sigmod_final` and compile the LaTeX document to generate a PDF. The figures in the PDF should match those in `test_small/source_latex/sample.pdf`.

## Table of Contents

- [Installation](#installation)
- [Download and Preprocessing](#download-and-preprocessing)
- [Iterative Algorithms](#iterative-algorithms)
- [Parameters](#parameters)
- [Time and Memory](#time-and-memory)
- [License](#license)

## Installation

Run `install.sh` to install necessary dependencies.

Required packages:
- Python package `matplotlib` for figure generation.
- The Boost library (available at [https://www.boost.org/](https://www.boost.org/)) for compiling and running the exact decomposition algorithm. For convenience, the Boost package is included in the repository and will be unzipped during installation.

## Download and Preprocessing

Execute `download_process.sh` to download the graph dataset from the internet, preprocess it into the required format, and place the processed files into the appropriate directory.

## Iterative Algorithms

Use `run_iterative.sh` to execute iterative algorithms, collect data on errors and inversion counts, and generate corresponding figures.

## Parameters

Run `parameters.sh` to conduct an ablation study on the parameter $\gamma_t$.

## Time and Memory

Execute `time_mem_script.sh` to collect performance data on runtime and memory usage across different algorithms.

## Hardware Specifications

Experiments were conducted on a Linux server with an Intel Xeon(R) Silver 4114 CPU running at 2.20 GHz, with main memory limited to 57 GB.
