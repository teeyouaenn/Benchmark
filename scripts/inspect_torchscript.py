#!/usr/bin/env python3
"""Read-only TorchScript signature and smoke-test inspector."""

from __future__ import annotations

import argparse

import torch


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    parser.add_argument("--shape", nargs="+", type=int)
    args = parser.parse_args()
    model = torch.jit.load(args.path, map_location="cpu")
    print(model)
    print(model.graph)
    if args.shape:
        x = torch.zeros(tuple(args.shape), dtype=torch.float32)
        with torch.inference_mode():
            y = model(x)
        print("output", type(y), getattr(y, "shape", None))
        if isinstance(y, (tuple, list)):
            print([getattr(item, "shape", None) for item in y])


if __name__ == "__main__":
    main()
