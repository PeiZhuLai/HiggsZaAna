#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path


DEFAULT_MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "output" / "latexTable_BrHZa.txt"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute Br(H->Za) for the Run 3 ALP mass grid and write a LaTeX table."
    )
    parser.add_argument(
        "--masses",
        nargs="+",
        type=float,
        default=DEFAULT_MASSES,
        help="ALP masses in GeV.",
    )
    parser.add_argument(
        "--mh",
        type=float,
        default=125.38,
        help="Higgs boson mass in GeV.",
    )
    parser.add_argument(
        "--mz",
        type=float,
        default=91.18,
        help="Z boson mass in GeV.",
    )
    parser.add_argument(
        "--c-zh-over-lambda",
        type=float,
        default=2.0e-5,
        help="Effective coupling c_Zh/Lambda in GeV^-1.",
    )
    parser.add_argument(
        "--gamma-sm",
        type=float,
        default=4.1e-3,
        help="SM Higgs width in GeV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output LaTeX table path.",
    )
    return parser.parse_args()


def kallen_lambda(x: float, y: float) -> float:
    return (1.0 - x - y) ** 2 - 4.0 * x * y


def gamma_h_to_za(ma: float, mh: float, mz: float, coupling: float) -> float:
    x = (mz / mh) ** 2
    y = (ma / mh) ** 2
    phase_space = kallen_lambda(x, y)
    if phase_space <= 0.0:
        raise ValueError(f"Mass point mA={ma:g} GeV is outside the kinematic phase space.")
    return (mh**3 / (16.0 * math.pi)) * (coupling**2) * (phase_space ** 1.5)


def branching_ratio(ma: float, mh: float, mz: float, coupling: float, gamma_sm: float) -> float:
    gamma_bsm = gamma_h_to_za(ma, mh, mz, coupling)
    return gamma_bsm / (gamma_bsm + gamma_sm)


def latex_scientific(value: float, decimals: int = 2) -> str:
    exponent = int(math.floor(math.log10(abs(value))))
    coefficient = value / (10.0**exponent)
    return rf"${coefficient:.{decimals}f} \times 10^{{{exponent}}}$"


def format_mass_label(mass: float) -> str:
    return str(int(mass)) if float(mass).is_integer() else f"{mass:g}"


def build_table(
    masses: list[float],
    mh: float,
    mz: float,
    coupling: float,
    gamma_sm: float,
) -> str:
    rows = []
    for mass in masses:
        br = branching_ratio(mass, mh, mz, coupling, gamma_sm)
        rows.append(f"{format_mass_label(mass)} & {latex_scientific(br)} \\\\")

    body = "\n".join(rows)
    return (
        "\\begin{table}[htbp]\n"
        "  \\centering\n"
        "  \\caption{Branching fraction of \\hToZa for the Run 3 signal mass grid, "
        "assuming $\\cZhLambda = 2 \\times 10^{-5}\\GeV^{-1}$.}\n"
        "  \\label{tab:BrHZa}\n"
        "  \\begin{tabular}{|c|c|}\n"
        "    \\hline\n"
        "    $m_{a}$ [GeV] & $\\mathrm{Br}(\\mathrm{H} \\rightarrow \\mathrm{Z}a)$ \\\\\n"
        "    \\hline\\hline\n"
        f"{body}\n"
        "    \\hline\n"
        "  \\end{tabular}\n"
        "\\end{table}\n"
    )


def main() -> None:
    args = parse_args()
    output_path = args.output.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    latex = build_table(
        masses=list(args.masses),
        mh=args.mh,
        mz=args.mz,
        coupling=args.c_zh_over_lambda,
        gamma_sm=args.gamma_sm,
    )
    output_path.write_text(latex)
    print(f"[write] {output_path}")


if __name__ == "__main__":
    main()
