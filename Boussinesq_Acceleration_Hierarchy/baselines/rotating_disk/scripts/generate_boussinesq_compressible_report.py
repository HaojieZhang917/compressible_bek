#!/usr/bin/env python3
"""Generate base-flow comparison data and a short report.

This script follows the basic-flow construction used by Bone.ipynb:

* incompressible + Boussinesq: Bone.py equations, solved with continuation;
* compressible path: CRD_STA.jl/baseflow_var + T_ca + physical interpolation,
  reimplemented here with SciPy to avoid Julia precompilation overhead.

Environment variables:
    TW_VALUES       Default: 1.0:0.1:2.0
    MR              Default: 0.3
    RO              Default: -1.0
    GAMMA           Default: 1.4
    PR              Default: 0.72
    SIGMA           Default: 0.72
    N_COMMON        Default: 2001
    N_RAW_INC       Default: 2001
    N_CHEB          Default: 199
    OUT_DIR         Default: baseflow_comparison_data
    BOUSS_STEP      Default: 0.02
    BVP_TOL         Default: 1e-6
"""

from __future__ import annotations

import csv
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.integrate import solve_bvp


ROOT = Path(__file__).resolve().parent.parent
FEATURE_Z_MIN = 0.2


def parse_float(name: str, default: float) -> float:
    return float(os.environ.get(name, str(default)))


def parse_int(name: str, default: int) -> int:
    return int(os.environ.get(name, str(default)))


def float_range(start: float, step: float, stop: float) -> list[float]:
    if step == 0.0:
        raise ValueError("step cannot be zero")
    out: list[float] = []
    x = start
    if step > 0:
        while x <= stop + 0.5 * abs(step):
            out.append(round(x, 12))
            x += step
    else:
        while x >= stop - 0.5 * abs(step):
            out.append(round(x, 12))
            x += step
    return out


def parse_tw_values(spec: str) -> list[float]:
    spec = spec.strip()
    if ":" in spec:
        parts = [float(p) for p in spec.split(":")]
        if len(parts) == 2:
            return float_range(parts[0], 0.1, parts[1])
        if len(parts) == 3:
            return float_range(parts[0], parts[1], parts[2])
        raise ValueError("TW_VALUES must be start:stop or start:step:stop")
    return [float(p) for p in spec.split(",") if p.strip()]


def write_rows(path: Path, header: Iterable[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(list(header))
        writer.writerows(rows)


def fmt_tw(tw: float) -> str:
    return f"{tw:.4f}".replace(".", "p")


@dataclass
class Settings:
    tw_values: list[float]
    mr: float
    ro: float
    gamma: float
    pr: float
    sigma: float
    n_common: int
    n_raw_inc: int
    n_cheb: int
    zmax_inc: float
    zmax_common: float
    bouss_step: float
    bvp_tol: float
    out_dir: Path


@dataclass
class Flow:
    z: np.ndarray
    f: np.ndarray
    g: np.ndarray
    h: np.ndarray
    t: np.ndarray


@dataclass
class BoussinesqFlow(Flow):
    df: np.ndarray
    dg: np.ndarray
    dh: np.ndarray
    dt: np.ndarray
    info: dict[str, float]


def load_settings() -> Settings:
    out_dir = Path(os.environ.get("OUT_DIR", str(ROOT / "baseflow_comparison_data")))
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    return Settings(
        tw_values=parse_tw_values(os.environ.get("TW_VALUES", "1.0:0.1:2.0")),
        mr=parse_float("MR", 0.3),
        ro=parse_float("RO", -1.0),
        gamma=parse_float("GAMMA", 1.4),
        pr=parse_float("PR", 0.72),
        sigma=parse_float("SIGMA", 0.72),
        n_common=parse_int("N_COMMON", 2001),
        n_raw_inc=parse_int("N_RAW_INC", 2001),
        n_cheb=parse_int("N_CHEB", 199),
        zmax_inc=parse_float("ZMAX_INC", 20.0),
        zmax_common=parse_float("Z_COMMON_MAX", 20.0),
        bouss_step=parse_float("BOUSS_STEP", 0.02),
        bvp_tol=parse_float("BVP_TOL", 1e-6),
        out_dir=out_dir,
    )


def boussinesq_guess(z: np.ndarray) -> np.ndarray:
    g = np.zeros((7, z.size))
    g[0] = -0.8845 * (1.0 - np.exp(-0.8 * z))
    g[1] = 0.5102 * (1.0 - 0.8 * z) * np.exp(-0.8 * z)
    g[2] = 0.5102 * z * np.exp(-0.8 * z)
    g[3] = 0.6159 * np.exp(-0.8 * z)
    g[4] = 1.0 - np.exp(-0.8 * z)
    g[5] = 0.0
    g[6] = 1.0
    return g


def boussinesq_ode_factory(pr: float, lambda_c: float = 1.0):
    def ode(_z: np.ndarray, y: np.ndarray) -> np.ndarray:
        fpp = y[2] ** 2 + y[0] * y[1] - (y[4] - 1.0) ** 2 - lambda_c * (y[6] - 1.0)
        gpp = 2.0 * y[2] * y[4] + y[0] * y[3] - 2.0 * y[2]
        return np.vstack(
            [
                -2.0 * y[2],
                fpp,
                y[1],
                gpp,
                y[3],
                pr * y[0] * y[5],
                y[5],
            ]
        )

    return ode


def boussinesq_bc_factory(tw: float):
    def bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
        return np.array(
            [
                ya[0],
                ya[2],
                ya[4],
                ya[6] - tw,
                yb[2],
                yb[4] - 1.0,
                yb[6] - 1.0,
            ]
        )

    return bc


def solve_boussinesq_sequence(settings: Settings) -> dict[float, BoussinesqFlow]:
    targets = sorted(set(settings.tw_values))
    if any(t < 1.0 for t in targets):
        raise ValueError("This continuation helper is configured for Tw >= 1.0.")

    z_mesh = np.linspace(0.0, settings.zmax_inc, 500)
    ode = boussinesq_ode_factory(settings.pr, lambda_c=1.0)
    sol = solve_bvp(
        ode,
        boussinesq_bc_factory(1.0),
        z_mesh,
        boussinesq_guess(z_mesh),
        tol=settings.bvp_tol,
        max_nodes=200000,
    )
    if not sol.success:
        raise RuntimeError(f"Boussinesq Tw=1 solve failed: {sol.message}")

    out: dict[float, BoussinesqFlow] = {}
    current_tw = 1.0
    current_sol = sol

    def store(tw: float, sol_obj) -> None:
        z = np.linspace(0.0, settings.zmax_inc, settings.n_raw_inc)
        y = sol_obj.sol(z)
        info = {
            "Tw": float(tw),
            "Pr": settings.pr,
            "lambda_c": 1.0,
            "Hinf": float(sol_obj.sol(settings.zmax_inc)[0]),
            "Fp0": float(sol_obj.sol(0.0)[1]),
            "Gp0": float(sol_obj.sol(0.0)[3]),
            "Tp0": float(sol_obj.sol(0.0)[5]),
        }
        out[round(tw, 12)] = BoussinesqFlow(
            z=z,
            h=y[0],
            df=y[1],
            f=y[2],
            dg=y[3],
            g=y[4],
            dt=y[5],
            t=y[6],
            dh=-2.0 * y[2],
            info=info,
        )

    if 1.0 in targets:
        store(1.0, current_sol)

    for target in targets:
        if abs(target - current_tw) < 1e-12:
            continue

        while current_tw < target - 1e-12:
            step = min(settings.bouss_step, target - current_tw)
            accepted = False
            while step >= 1e-4 and not accepted:
                next_tw = current_tw + step
                guess = current_sol.sol(z_mesh)
                trial = solve_bvp(
                    ode,
                    boussinesq_bc_factory(next_tw),
                    z_mesh,
                    guess,
                    tol=settings.bvp_tol,
                    max_nodes=200000,
                )
                if trial.success:
                    current_tw = next_tw
                    current_sol = trial
                    accepted = True
                else:
                    step *= 0.5
            if not accepted:
                raise RuntimeError(f"Boussinesq continuation failed near Tw={current_tw:.6f}")

        store(target, current_sol)
        print(f"  Boussinesq Tw={target:.4f} solved", flush=True)

    return out


def compressible_velocity_ode_factory(ro: float):
    co = 2.0 - ro - ro**2

    def ode(_z: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.vstack(
            [
                -2.0 * y[2],
                ro * (y[2] ** 2 + y[0] * y[1] - (y[4] ** 2 - 1.0)) - co * (y[4] - 1.0),
                y[1],
                ro * (2.0 * y[2] * y[4] + y[0] * y[3]) + co * y[2],
                y[3],
            ]
        )

    return ode


def compressible_velocity_bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
    return np.array([ya[0], ya[2], ya[4], yb[2], yb[4] - 1.0])


@dataclass
class CompressibleBasePieces:
    eta: np.ndarray
    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    du: np.ndarray
    dv: np.ndarray
    phi: np.ndarray
    f_heat: np.ndarray
    q_wall: np.ndarray


def solve_compressible_pieces(settings: Settings) -> CompressibleBasePieces:
    z = np.linspace(0.0, 40.0, 500)
    guess = np.zeros((5, z.size))
    guess[0] = 1.2
    guess[4] = 1.0
    sol = solve_bvp(
        compressible_velocity_ode_factory(settings.ro),
        compressible_velocity_bc,
        z,
        guess,
        tol=1e-8,
        max_nodes=500000,
    )
    if not sol.success:
        raise RuntimeError(f"Compressible velocity BVP failed: {sol.message}")

    eta = np.linspace(0.0, 40.0, 10001)
    y = sol.sol(eta)
    u = -y[2]
    v = -y[4]
    w = -y[0]
    du = -y[1]
    dv = -y[3]
    if settings.ro >= 0.0:
        u = -u
        v = -v
        w = -w
        du = -du
        dv = -dv

    delta = eta[1] - eta[0]
    phi = np.cumsum(u) * delta
    f_heat, q_wall = solve_heat_functions(settings, eta, u, du, dv, phi)
    return CompressibleBasePieces(eta=eta, u=u, v=v, w=w, du=du, dv=dv, phi=phi, f_heat=f_heat, q_wall=q_wall)


def solve_heat_functions(
    settings: Settings,
    eta: np.ndarray,
    u: np.ndarray,
    du: np.ndarray,
    dv: np.ndarray,
    phi: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mesh = np.linspace(0.0, 40.0, 600)

    def interp(arr: np.ndarray, x: np.ndarray) -> np.ndarray:
        return np.interp(x, eta, arr)

    def f_ode(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        ux = interp(u, x)
        dux = interp(du, x)
        dvx = interp(dv, x)
        phix = interp(phi, x)
        return np.vstack(
            [
                y[1],
                2.0 * settings.sigma * (dux**2 + dvx**2 + ux * y[0] - phix * y[1]),
            ]
        )

    def f_bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
        return np.array([ya[0], yb[0]])

    f_guess = np.zeros((2, mesh.size))
    f_sol = solve_bvp(f_ode, f_bc, mesh, f_guess, tol=1e-6, max_nodes=100000)
    if not f_sol.success:
        raise RuntimeError(f"Compressible heat f BVP failed: {f_sol.message}")

    def q_ode(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        phix = interp(phi, x)
        return np.vstack([y[1], -2.0 * settings.sigma * phix * y[1]])

    def q_bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
        return np.array([ya[0] - 1.0, yb[0]])

    q_guess = np.zeros((2, mesh.size))
    q_guess[0] = 1.0 - mesh / mesh[-1]
    q_guess[1] = -1.0 / mesh[-1]
    q_sol = solve_bvp(q_ode, q_bc, mesh, q_guess, tol=1e-6, max_nodes=100000)
    if not q_sol.success:
        raise RuntimeError(f"Compressible heat q BVP failed: {q_sol.message}")

    return f_sol.sol(eta)[0], q_sol.sol(eta)[0]


def physical_coordinate(t: np.ndarray, delta: float) -> np.ndarray:
    z = np.zeros_like(t)
    csum = np.cumsum(delta * t)
    z[1:] = csum[1:]
    return z


def compressible_flow_for_tw(settings: Settings, pieces: CompressibleBasePieces, tw: float) -> Flow:
    temp = 1.0 - ((settings.gamma - 1.0) / 2.0) * settings.mr**2 * pieces.f_heat + (tw - 1.0) * pieces.q_wall
    h = pieces.w * temp
    z_phys = physical_coordinate(temp, pieces.eta[1] - pieces.eta[0])
    return Flow(z=z_phys, f=pieces.u.copy(), g=pieces.v.copy(), h=h, t=temp)


def unique_xy(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(x)
    xs = np.asarray(x)[order]
    ys = np.asarray(y)[order]
    keep_x = []
    keep_y = []
    i = 0
    while i < xs.size:
        j = i + 1
        while j < xs.size and xs[j] == xs[i]:
            j += 1
        keep_x.append(xs[j - 1])
        keep_y.append(ys[j - 1])
        i = j
    return np.asarray(keep_x), np.asarray(keep_y)


def interp_to(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    xs, ys = unique_xy(x, y)
    return np.interp(x_new, xs, ys, left=ys[0], right=ys[-1])


def cheb_physical_points(n: int) -> np.ndarray:
    theta = np.linspace(0.0, math.pi, n + 1)
    xi = -np.cos(theta)
    a = 6.0
    b = 0.6
    c = 0.5
    numerator = a * (1.0 + b * xi + (1.0 - b) * (xi**3 + c * (1.0 - xi**2)))
    denominator = 1.0 - b * xi - (1.0 - b) * (xi**3 + c * (1.0 - xi**2))
    with np.errstate(divide="ignore", invalid="ignore"):
        x = numerator / denominator
    x[x > 20.0] = 20.0
    return x


def profile_derivatives(z: np.ndarray, y: np.ndarray, z_eval: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = interp_to(z, y, z_eval)
    d1 = np.gradient(values, z_eval, edge_order=2)
    d2 = np.gradient(d1, z_eval, edge_order=2)
    return values, d1, d2


def zero_crossings(z: np.ndarray, y: np.ndarray, max_count: int = 8) -> list[float]:
    vals: list[float] = []
    scale = max(float(np.nanmax(np.abs(y))), 1e-14)
    if scale < 1e-8:
        return vals
    threshold = scale * 5e-3
    last_z = None
    last_y = None
    last_sign = 0
    for zi, yi in zip(z, y):
        if zi <= z[0] + FEATURE_Z_MIN or zi >= z[-1] - 0.05:
            continue
        if not np.isfinite(yi) or abs(yi) < threshold:
            continue
        sign = 1 if yi > 0.0 else -1
        if last_sign != 0 and sign != last_sign and last_z is not None and last_y is not None:
            root = last_z - last_y * (zi - last_z) / (yi - last_y)
            vals.append(float(root))
            if len(vals) >= max_count:
                break
        last_z = zi
        last_y = yi
        last_sign = sign
    cleaned: list[float] = []
    for val in vals:
        if not cleaned or abs(val - cleaned[-1]) > 0.05:
            cleaned.append(val)
    return cleaned


def feature_rows_for_model(tw: float, model: str, z: np.ndarray, values: dict[str, np.ndarray]) -> list[list[object]]:
    rows: list[list[object]] = []
    interior = (z >= z[0] + FEATURE_Z_MIN) & (z <= z[-1] - 0.05)
    for name, arr in values.items():
        y, dy, d2y = profile_derivatives(z, arr, z)
        inflections = zero_crossings(z, d2y)
        idx_dy = int(np.argmax(np.abs(dy)))
        if np.any(interior):
            interior_indices = np.flatnonzero(interior)
            idx_d2y = int(interior_indices[np.argmax(np.abs(d2y[interior]))])
        else:
            idx_d2y = int(np.argmax(np.abs(d2y)))
        rows.append(
            [
                tw,
                model,
                name,
                y[0],
                y[-1],
                dy[0],
                dy[-1],
                np.max(np.abs(dy)),
                z[idx_dy],
                np.max(np.abs(d2y)),
                z[idx_d2y],
                len(inflections),
                ";".join(f"{x:.6f}" for x in inflections),
            ]
        )
    return rows


def build_report(
    settings: Settings,
    summary_rows: list[dict[str, object]],
    feature_rows: list[list[object]],
    report_path: Path,
) -> None:
    by_tw = {float(r["Tw"]): r for r in summary_rows}
    selected = [tw for tw in [1.0, 1.1, 1.2, 1.5, 2.0] if tw in by_tw]

    def metric_row(tw: float) -> str:
        r = by_tw[tw]
        return (
            f"| {tw:.2f} | {float(r['maxabs_dF_common']):.4e} | "
            f"{float(r['maxabs_dG_signmatched_common']):.4e} | "
            f"{float(r['maxabs_dH_common']):.4e} | "
            f"{float(r['maxabs_dT_common']):.4e} | "
            f"{float(r['rel_rho_wall_change']):.3f} | "
            f"{float(r['rho_linearization_rel_error_wall']):.3f} |"
        )

    feature_header = [
        "Tw",
        "model",
        "variable",
        "wall_value",
        "edge_value",
        "wall_d1",
        "edge_d1",
        "max_abs_d1",
        "z_max_abs_d1",
        "max_abs_d2",
        "z_max_abs_d2",
        "inflection_count",
        "inflection_z",
    ]
    idx = {name: i for i, name in enumerate(feature_header)}

    def find_feature(tw: float, model: str, variable: str) -> list[object] | None:
        for row in feature_rows:
            if abs(float(row[idx["Tw"]]) - tw) < 1e-12 and row[idx["model"]] == model and row[idx["variable"]] == variable:
                return row
        return None

    inflection_lines: list[str] = []
    for tw in selected:
        for variable in ["F", "G_compare", "H", "T"]:
            bi = find_feature(tw, "Boussinesq", variable)
            co = find_feature(tw, "Compressible", variable)
            if bi is None or co is None:
                continue
            inflection_lines.append(
                "| "
                + " | ".join(
                    [
                        f"{tw:.2f}",
                        variable,
                        f"{float(bi[idx['wall_d1']]):.4e}",
                        f"{float(co[idx['wall_d1']]):.4e}",
                        f"{float(bi[idx['max_abs_d1']]):.4e}",
                        f"{float(co[idx['max_abs_d1']]):.4e}",
                        str(bi[idx["inflection_z"]]) or "-",
                        str(co[idx["inflection_z"]]) or "-",
                    ]
                )
                + " |"
            )

    eps05 = [r for r in summary_rows if float(r["max_rel_velocity_diff"]) <= 0.05]
    eps10 = [r for r in summary_rows if float(r["max_rel_velocity_diff"]) <= 0.10]
    range05 = max((float(r["Tw"]) for r in eps05), default=None)
    range10 = max((float(r["Tw"]) for r in eps10), default=None)

    lines = [
        "# Base-flow comparison report",
        "",
        "## Scope",
        "",
        "This report was generated from the scripts in `Bone.ipynb`/`CRD_STA.jl` logic, but uses a SciPy reimplementation to avoid the long Julia precompile path.",
        "",
        f"* Wall-temperature sweep: `{min(settings.tw_values):.3f} <= Tw <= {max(settings.tw_values):.3f}`.",
        f"* Compressible settings: `Ro={settings.ro}`, `Mr={settings.mr}`, `gamma={settings.gamma}`, `sigma={settings.sigma}`.",
        f"* Boussinesq settings: `Pr={settings.pr}`, ideal-gas centrifugal buoyancy coefficient `lambda_c=1`.",
        f"* Common comparison coordinate: physical `z in [0, {settings.zmax_common}]`, `N={settings.n_common}`.",
        "",
        "Because the current compressible notebook uses `Ro=-1`, its azimuthal velocity approaches `G_c(infinity)=-1`. The comparison therefore uses `G_compare=-G_i` for the Boussinesq solution while still writing the raw `G_i` columns to the CSV files.",
        "",
        "## Main Differences",
        "",
        "The Boussinesq model keeps density constant everywhere except the centrifugal buoyancy term. In the present nondimensionalization, the radial-momentum feedback is proportional to `-(T-1)`, so changing wall temperature directly changes the velocity field.",
        "",
        "The compressible basic-flow path currently used by `Bone.ipynb` solves the isothermal BEK velocity field first, then constructs temperature, density, axial velocity and physical coordinate through `T_ca` and `z=int T deta`. Therefore the wall temperature mainly changes density, axial scaling and coordinate stretching; the raw tangential/radial similarity velocity is not recomputed from a temperature-coupled momentum system in this script.",
        "",
        "## Difference Metrics",
        "",
            "| Tw | max abs dF | max abs dG compare | max abs dH | max abs dT | wall density change | rho linearization error |",
            "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    lines.extend(metric_row(tw) for tw in selected)
    lines.extend(
        [
            "",
            "Here `wall density change = |1/Tw - 1|`, the ideal-gas relative density change at the wall. `rho linearization error` compares exact ideal-gas `rho_w/rho_inf=1/Tw` with the Boussinesq linearization `1-(Tw-1)`.",
            "",
            "## Derivatives And Inflection Points",
            "",
            "| Tw | variable | Bouss wall d1 | Comp wall d1 | Bouss max abs d1 | Comp max abs d1 | Bouss inflection z | Comp inflection z |",
            "|---:|:---|---:|---:|---:|---:|:---|:---|",
        ]
    )
    lines.extend(inflection_lines)
    lines.extend(
        [
            "",
            f"The inflection locations are detected from finite-difference second derivatives on the common uniform physical grid. Second-derivative extrema are written to `feature_summary.csv`, with their maximum taken from the interior grid (`z>={FEATURE_Z_MIN}`) to avoid endpoint artifacts. These values should be used as diagnostics of structural change, not as theorem-level locations.",
            "",
            "## Applicability Range",
            "",
            "For an ideal gas with `T*=T/T_inf`, the Boussinesq density linearization is `rho/rho_inf ~= 1-(T-1)`. Its formal small parameter is therefore `epsilon=|Tw-1|` near the wall.",
            "",
            "A conservative interpretation is:",
            "",
            "* `|Tw-1| <= 0.05`: usually a safe Boussinesq range for quantitative stability comparisons.",
            "* `0.05 < |Tw-1| <= 0.10`: useful for trend studies, but basic-flow and derivative differences should be checked case by case.",
            "* `|Tw-1| > 0.10-0.20`: no longer a clean Boussinesq limit for an ideal gas, especially in a rotating flow where the buoyancy term enters the radial momentum balance directly.",
            "* `Tw=2` gives `epsilon=1` and a wall density ratio `rho_w/rho_inf=0.5`; this is outside the traditional Boussinesq asymptotic range.",
            "",
        ]
    )
    if range05 is not None or range10 is not None:
        lines.append("Using the current data-driven maximum relative velocity difference metric:")
        if range05 is not None:
            lines.append(f"* the 5% criterion is satisfied up to approximately `Tw={range05:.2f}` in this sweep;")
        else:
            lines.append("* the 5% criterion is not satisfied even at the first non-isothermal point in this sweep;")
        if range10 is not None:
            lines.append(f"* the 10% criterion is satisfied up to approximately `Tw={range10:.2f}` in this sweep.")
        else:
            lines.append("* the 10% criterion is not satisfied even at the first non-isothermal point in this sweep.")
        lines.append("")
    lines.extend(
        [
            "For the paper logic, this suggests treating Boussinesq as the low-temperature-difference baseline and then using the compressible model once density variation, coordinate stretching, derivative changes or inflection-point migration become non-negligible.",
            "",
            "Important caveat: the data-driven velocity-difference range above is a comparison against the current `Bone.ipynb` compressible basic-flow construction. Since that construction does not recompute the radial/tangential velocity field with temperature-coupled density and viscosity in the momentum equations, it should not be quoted as a universal physical failure point of the Boussinesq approximation.",
            "",
            "## Output Files",
            "",
            "* `incompressible_boussinesq_physical_raw.csv`: raw Boussinesq profiles in physical/similarity `z`.",
            "* `compressible_physical_raw.csv`: raw compressible profiles in physical `z=int T deta`.",
            "* `physical_common_grid_interpolated.csv`: both models interpolated to the same physical grid, with differences.",
            "* `cheb_grid_interpolated.csv`: both models interpolated to the notebook Chebyshev physical grid.",
            "* `derivatives_common_grid.csv`: first and second derivatives on the common physical grid.",
            "* `feature_summary.csv`: wall/edge values, derivative extrema and inflection-point summaries.",
            "* `summary.csv`: one-line difference metrics for each `Tw`.",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    settings = load_settings()
    settings.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {settings.out_dir}", flush=True)
    print("Solving compressible base pieces...", flush=True)
    comp_pieces = solve_compressible_pieces(settings)
    print("Solving Boussinesq continuation...", flush=True)
    bouss = solve_boussinesq_sequence(settings)

    z_common = np.linspace(0.0, settings.zmax_common, settings.n_common)
    x_cheb = cheb_physical_points(settings.n_cheb)

    raw_inc_rows: list[list[object]] = []
    raw_comp_rows: list[list[object]] = []
    common_rows: list[list[object]] = []
    cheb_rows: list[list[object]] = []
    deriv_rows: list[list[object]] = []
    feature_rows: list[list[object]] = []
    summary_rows: list[dict[str, object]] = []

    for tw in settings.tw_values:
        key = round(tw, 12)
        bi = bouss[key]
        co = compressible_flow_for_tw(settings, comp_pieces, tw)
        rho_c = 1.0 / co.t

        for j in range(bi.z.size):
            raw_inc_rows.append(
                [
                    tw,
                    bi.z[j],
                    bi.f[j],
                    bi.g[j],
                    -bi.g[j],
                    bi.h[j],
                    bi.t[j],
                    bi.df[j],
                    bi.dg[j],
                    bi.dh[j],
                    bi.dt[j],
                    bi.info["Hinf"],
                    bi.info["Fp0"],
                    bi.info["Gp0"],
                    bi.info["Tp0"],
                ]
            )
        for j in range(co.z.size):
            raw_comp_rows.append([tw, co.z[j], co.f[j], co.g[j], co.h[j], co.t[j], rho_c[j]])

        f_i = interp_to(bi.z, bi.f, z_common)
        g_i = interp_to(bi.z, bi.g, z_common)
        h_i = interp_to(bi.z, bi.h, z_common)
        t_i = interp_to(bi.z, bi.t, z_common)
        f_c = interp_to(co.z, co.f, z_common)
        g_c = interp_to(co.z, co.g, z_common)
        h_c = interp_to(co.z, co.h, z_common)
        t_c = interp_to(co.z, co.t, z_common)
        rho_c_common = 1.0 / t_c

        for j in range(z_common.size):
            common_rows.append(
                [
                    tw,
                    z_common[j],
                    f_i[j],
                    g_i[j],
                    -g_i[j],
                    h_i[j],
                    t_i[j],
                    f_c[j],
                    g_c[j],
                    h_c[j],
                    t_c[j],
                    rho_c_common[j],
                    f_c[j] - f_i[j],
                    g_c[j] - g_i[j],
                    g_c[j] + g_i[j],
                    h_c[j] - h_i[j],
                    t_c[j] - t_i[j],
                ]
            )

        f_i_ch = interp_to(bi.z, bi.f, x_cheb)
        g_i_ch = interp_to(bi.z, bi.g, x_cheb)
        h_i_ch = interp_to(bi.z, bi.h, x_cheb)
        t_i_ch = interp_to(bi.z, bi.t, x_cheb)
        f_c_ch = interp_to(co.z, co.f, x_cheb)
        g_c_ch = interp_to(co.z, co.g, x_cheb)
        h_c_ch = interp_to(co.z, co.h, x_cheb)
        t_c_ch = interp_to(co.z, co.t, x_cheb)
        rho_c_ch = 1.0 / t_c_ch
        for j in range(x_cheb.size):
            cheb_rows.append(
                [
                    tw,
                    x_cheb[j],
                    f_i_ch[j],
                    g_i_ch[j],
                    -g_i_ch[j],
                    h_i_ch[j],
                    t_i_ch[j],
                    f_c_ch[j],
                    g_c_ch[j],
                    h_c_ch[j],
                    t_c_ch[j],
                    rho_c_ch[j],
                    f_c_ch[j] - f_i_ch[j],
                    g_c_ch[j] - g_i_ch[j],
                    g_c_ch[j] + g_i_ch[j],
                    h_c_ch[j] - h_i_ch[j],
                    t_c_ch[j] - t_i_ch[j],
                ]
            )

        b_values = {"F": f_i, "G_compare": -g_i, "H": h_i, "T": t_i}
        c_values = {"F": f_c, "G_compare": g_c, "H": h_c, "T": t_c}
        feature_rows.extend(feature_rows_for_model(tw, "Boussinesq", z_common, b_values))
        feature_rows.extend(feature_rows_for_model(tw, "Compressible", z_common, c_values))

        for model, values in [("Boussinesq", b_values), ("Compressible", c_values)]:
            for name, arr in values.items():
                y, dy, d2y = profile_derivatives(z_common, arr, z_common)
                for j in range(z_common.size):
                    deriv_rows.append([tw, model, name, z_common[j], y[j], dy[j], d2y[j]])

        vel_scale = max(
            float(np.max(np.abs(f_c))),
            float(np.max(np.abs(g_c))),
            float(np.max(np.abs(h_c))),
            1e-12,
        )
        max_rel_vel_diff = max(
            float(np.max(np.abs(f_c - f_i))),
            float(np.max(np.abs(g_c + g_i))),
            float(np.max(np.abs(h_c - h_i))),
        ) / vel_scale

        summary_rows.append(
            {
                "Tw": tw,
                "Ro": settings.ro,
                "Mr": settings.mr,
                "N_common": settings.n_common,
                "N_cheb": settings.n_cheb,
                "Hinf_i": bi.info["Hinf"],
                "Fp0_i": bi.info["Fp0"],
                "Gp0_i": bi.info["Gp0"],
                "Tp0_i": bi.info["Tp0"],
                "H_c_at_zmax": h_c[-1],
                "F_c_at_zmax": f_c[-1],
                "G_c_at_zmax": g_c[-1],
                "T_c_at_zmax": t_c[-1],
                "maxabs_dF_common": float(np.max(np.abs(f_c - f_i))),
                "maxabs_dG_raw_common": float(np.max(np.abs(g_c - g_i))),
                "maxabs_dG_signmatched_common": float(np.max(np.abs(g_c + g_i))),
                "maxabs_dH_common": float(np.max(np.abs(h_c - h_i))),
                "maxabs_dT_common": float(np.max(np.abs(t_c - t_i))),
                "maxabs_dF_cheb": float(np.max(np.abs(f_c_ch - f_i_ch))),
                "maxabs_dG_raw_cheb": float(np.max(np.abs(g_c_ch - g_i_ch))),
                "maxabs_dG_signmatched_cheb": float(np.max(np.abs(g_c_ch + g_i_ch))),
                "maxabs_dH_cheb": float(np.max(np.abs(h_c_ch - h_i_ch))),
                "maxabs_dT_cheb": float(np.max(np.abs(t_c_ch - t_i_ch))),
                "max_rel_velocity_diff": max_rel_vel_diff,
                "rel_rho_wall_change": abs(1.0 / tw - 1.0),
                "rho_linearization_rel_error_wall": abs((1.0 / tw) - (1.0 - (tw - 1.0))) / (1.0 / tw),
            }
        )
        print(f"  Data assembled for Tw={tw:.4f}", flush=True)

    write_rows(
        settings.out_dir / "incompressible_boussinesq_physical_raw.csv",
        [
            "Tw",
            "z",
            "F_i",
            "G_i",
            "minus_G_i",
            "H_i",
            "T_i",
            "dF_i",
            "dG_i",
            "dH_i",
            "dT_i",
            "Hinf_i",
            "Fp0_i",
            "Gp0_i",
            "Tp0_i",
        ],
        raw_inc_rows,
    )
    write_rows(
        settings.out_dir / "compressible_physical_raw.csv",
        ["Tw", "z_physical", "F_c", "G_c", "H_c", "T_c", "rho_c"],
        raw_comp_rows,
    )
    write_rows(
        settings.out_dir / "physical_common_grid_interpolated.csv",
        [
            "Tw",
            "z",
            "F_i",
            "G_i",
            "minus_G_i",
            "H_i",
            "T_i",
            "F_c",
            "G_c",
            "H_c",
            "T_c",
            "rho_c",
            "dF_c_minus_i",
            "dG_c_minus_i_raw",
            "dG_c_minus_minus_Gi",
            "dH_c_minus_i",
            "dT_c_minus_i",
        ],
        common_rows,
    )
    write_rows(
        settings.out_dir / "cheb_grid_interpolated.csv",
        [
            "Tw",
            "x_cheb",
            "F_i",
            "G_i",
            "minus_G_i",
            "H_i",
            "T_i",
            "F_c",
            "G_c",
            "H_c",
            "T_c",
            "rho_c",
            "dF_c_minus_i",
            "dG_c_minus_i_raw",
            "dG_c_minus_minus_Gi",
            "dH_c_minus_i",
            "dT_c_minus_i",
        ],
        cheb_rows,
    )
    write_rows(
        settings.out_dir / "derivatives_common_grid.csv",
        ["Tw", "model", "variable", "z", "value", "d1_dz", "d2_dz2"],
        deriv_rows,
    )
    write_rows(
        settings.out_dir / "feature_summary.csv",
        [
            "Tw",
            "model",
            "variable",
            "wall_value",
            "edge_value",
            "wall_d1",
            "edge_d1",
            "max_abs_d1",
            "z_max_abs_d1",
            "max_abs_d2",
            "z_max_abs_d2",
            "inflection_count",
            "inflection_z",
        ],
        feature_rows,
    )

    summary_header = list(summary_rows[0].keys())
    write_rows(
        settings.out_dir / "summary.csv",
        summary_header,
        ([row[name] for name in summary_header] for row in summary_rows),
    )

    metadata = [
        "Generated by generate_boussinesq_compressible_report.py",
        f"TW_VALUES={settings.tw_values}",
        f"MR={settings.mr}",
        f"RO={settings.ro}",
        f"GAMMA={settings.gamma}",
        f"PR={settings.pr}",
        f"SIGMA={settings.sigma}",
        f"N_COMMON={settings.n_common}",
        f"N_RAW_INC={settings.n_raw_inc}",
        f"N_CHEB={settings.n_cheb}",
        f"BOUSS_STEP={settings.bouss_step}",
        f"BVP_TOL={settings.bvp_tol}",
        "G comparison uses G_c - (-G_i), stored as dG_c_minus_minus_Gi.",
    ]
    (settings.out_dir / "metadata.txt").write_text("\n".join(metadata) + "\n", encoding="utf-8")
    build_report(settings, summary_rows, feature_rows, settings.out_dir / "baseflow_model_comparison_report.md")
    print(f"Done. Report: {settings.out_dir / 'baseflow_model_comparison_report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
