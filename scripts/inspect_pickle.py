#!/usr/bin/env python3
"""Small read-only helper for auditing serialized scientific resources."""

from __future__ import annotations

import argparse
import pickle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    args = parser.parse_args()

    with open(args.path, "rb") as handle:
        obj = pickle.load(handle)

    print(type(obj))
    if isinstance(obj, dict):
        print("keys", list(obj)[:30], "n", len(obj))
        for key in list(obj)[:5]:
            value = obj[key]
            print(key, type(value), getattr(value, "shape", None), getattr(value, "dtype", None))
            if hasattr(value, "nnz"):
                print("nnz", value.nnz, "min", value.min(), "max", value.max())
            elif hasattr(value, "size"):
                import numpy as np

                array = np.asarray(value)
                print(
                    "min",
                    float(np.nanmin(array)),
                    "max",
                    float(np.nanmax(array)),
                    "mean",
                    float(np.nanmean(array)),
                    "sum",
                    float(np.nansum(array)),
                    "nonzero",
                    int(np.count_nonzero(array)),
                )


if __name__ == "__main__":
    main()
