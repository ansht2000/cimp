"""
Batch testing script to isolate color error when seam carving the tower-full-scale.ppm image
"""

import subprocess, sys
from pathlib import Path

CIMP_EXE = "./cimp"
INPUT_IMAGE = Path("data/tower-full-scale.ppm")
OUT_DIR = Path("out")

def run_cimp(cols_scale: float, row_scale: float, output: Path):
    cmd = [
        CIMP_EXE,
        str(INPUT_IMAGE),
        str(output),
        "seam",
        f"{cols_scale}",
        f"{row_scale}"
    ]
    print("->", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode:
        print(res.stderr, file=sys.stderr)
        raise RuntimeError(f"cimp failed with code {res.returncode}")

def main():
    for i in range(40, 71):
        for j in range(40, 71):
            out_name = OUT_DIR.joinpath(Path(f"tower-full-scale-seam-{i / 100}-{j / 100}.ppm"))
            run_cimp(i / 100, j / 100, out_name)

if __name__ == "__main__":
    main()