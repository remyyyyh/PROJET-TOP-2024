#!/usr/bin/python3

import argparse
import numpy as np
import pandas as pd
from typing import Tuple, List


class StencilResults:
    def __init__(self, dims: Tuple[int, int, int], results: pd.DataFrame, runtime: List[float]):
        self.dims = dims
        self.results = results
        self.runtime = runtime

    def __str__(self):
        return f"Stencil results:\n\tDimensions: {self.dims}\n\tResults:\n{self.results}\n\tRuntimes:\n{self.runtime}\n\n"


def retrieve_results_and_runtime(file_path: str) -> Tuple[Tuple[int, int, int], pd.DataFrame, List[float]]:
    raw_data = pd.read_csv(file_path, header=None, delim_whitespace=True)
    simdims = tuple(map(int, raw_data.iloc[0, -3:].values))
    results = raw_data.iloc[:, 0].values
    runtime = raw_data.iloc[:, 1].values
    return StencilResults(simdims, results, runtime)


def compare(ref: StencilResults, res: StencilResults) -> None:
    # Verify dimensions
    if ref.dims != res.dims:
        raise ValueError("Reference and result dimensions do not match.")

    # Initialize a flag to track whether any differences were found
    any_difference = False
    # Initialize a list to store the positions of incorrect cells
    incorrect_cells = []
    # Compare each cell in the results
    for i in range(len(ref.results)):
            diff = abs(ref.results[i] - res.results[i])
            if diff > 1e-12:
                print(f"\x1b[1;33mwarning:\x1b[0m difference found at iteration {i + 1}: {diff:e}")
                any_difference = True
                incorrect_cells.append(i)

    # If any differences found, raise an error
    if any_difference:
        error_message = "Results diverge too much from reference. Incorrect iterations: "
        error_message += ", ".join([f"{cell + 1}" for cell in incorrect_cells])
        raise ValueError(error_message)

    # Compare runtimes
    avg_runtime_ref = np.mean(ref.runtime)
    avg_runtime_res = np.mean(res.runtime)
    
    acc = ((avg_runtime_ref / avg_runtime_res) - 1.0) * 100.0
    if avg_runtime_ref < avg_runtime_res:
        print(f"\x1b[31mReference is {-acc:.2f}% faster than result\x1b[0m")
    elif avg_runtime_ref > avg_runtime_res:
        print(f"\x1b[32mResult is {acc:.2f}% faster than reference\x1b[0m")
    else:
        print("Reference and result have the same average runtime")


def main():
    parser = argparse.ArgumentParser(description='Compare results of a stencil run to a reference.')
    parser.add_argument('reference', type=str, help='Path to the reference file')
    parser.add_argument('results', type=str, help='Path to the results file')
    args = parser.parse_args()

    reference = retrieve_results_and_runtime(args.reference)
    results = retrieve_results_and_runtime(args.results)
    compare(reference, results)


if __name__ == "__main__":
    main()
