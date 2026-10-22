"""Generuje ilustracje do raportu o rzucie ukośnym.

Uruchom: python generuj_wykresy.py
Wymagane pakiety: numpy, matplotlib, scipy.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from scipy.integrate import solve_ivp


# Podstaw te wartości demonstracyjne swoimi cyframi indeksu, jeśli wymagają
# tego instrukcje zadania: V0 w m/s i ALPHA_DEG w stopniach.
V0 = 25.0
ALPHA_DEG = 40.0
G = 9.81
MASS = 0.145
RHO = 1.225
CD = 0.47
AREA = 0.0042
B_LINEAR = 0.020

OUT = Path(__file__).resolve().parent / "files" / "rzut_ukosny"
OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, GREEN, INK = "#2463A6", "#D47722", "#168578", "#233044"
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.2,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def integrate_drag(angle_deg, model, sample_count=1200):
    angle = np.deg2rad(angle_deg)
    vx0, vy0 = V0 * np.cos(angle), V0 * np.sin(angle)

    def rhs(_time, state):
        x, y, vx, vy = state
        speed = np.hypot(vx, vy)
        if model == "linear":
            ax = -(B_LINEAR / MASS) * vx
            ay = -G - (B_LINEAR / MASS) * vy
        else:
            coefficient = 0.5 * RHO * CD * AREA * speed / MASS
            ax = -coefficient * vx
            ay = -G - coefficient * vy
        return vx, vy, ax, ay

    def landing(_time, state):
        return state[1]

    landing.terminal = True
    landing.direction = -1
    flight_estimate = 2 * V0 * np.sin(angle) / G
    solution = solve_ivp(
        rhs,
        (0.0, max(2.0, 4 * flight_estimate)),
        (0.0, 0.0, vx0, vy0),
        events=landing,
        dense_output=True,
        rtol=1e-9,
        atol=1e-11,
        max_step=0.02,
    )
    if not solution.success or len(solution.t_events[0]) != 1:
        raise RuntimeError(f"Nie znaleziono zdarzenia lądowania ({model}, {angle_deg}°).")
    flight_time = solution.t_events[0][0]
    initial_state = solution.sol(0.0)
    if not np.allclose(initial_state, (0.0, 0.0, vx0, vy0), rtol=0, atol=1e-12):
        raise AssertionError("Solver nie zachował zadanych warunków początkowych.")
    if abs(solution.y_events[0][0][1]) > 1e-8:
        raise AssertionError("Zdarzenie lądowania nie spełniło warunku y = 0.")
    times = np.linspace(0.0, flight_time, sample_count)
    states = solution.sol(times)
    states[1, -1] = 0.0
    return times, states


def vacuum(angle_deg, times):
    angle = np.deg2rad(angle_deg)
    vx0, vy0 = V0 * np.cos(angle), V0 * np.sin(angle)
    x = vx0 * times
    y = vy0 * times - 0.5 * G * times**2
    vx = np.full_like(times, vx0)
    vy = vy0 - G * times
    return np.vstack((x, y, vx, vy))


def flight_summary(angle_deg, model):
    if model == "vacuum":
        angle = np.deg2rad(angle_deg)
        flight_time = 2 * V0 * np.sin(angle) / G
        return flight_time, V0**2 * np.sin(2 * angle) / G
    times, states = integrate_drag(angle_deg, model, sample_count=2)
    return times[-1], states[0, -1]


def make_schematic():
    angle = np.deg2rad(ALPHA_DEG)
    time = np.linspace(0, 2 * V0 * np.sin(angle) / G, 300)
    states = vacuum(ALPHA_DEG, time)
    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    ax.plot(states[0], states[1], color=BLUE, label="Tor bez oporu")
    ax.axhline(0, color=INK, linewidth=1)
    ax.add_patch(Circle((0, 0), 0.8, color=ORANGE, zorder=4))
    scale = max(states[0, -1] * 0.18, 3.0)
    ax.quiver(0, 0, scale * np.cos(angle), scale * np.sin(angle),
              angles="xy", scale_units="xy", scale=1, color=INK,
              width=0.006, label="$v_0$")
    ax.quiver(0, 0, scale * np.cos(angle), 0, angles="xy",
              scale_units="xy", scale=1, color=GREEN, width=0.005)
    ax.quiver(0, 0, 0, scale * np.sin(angle), angles="xy",
              scale_units="xy", scale=1, color=ORANGE, width=0.005)
    ax.text(scale * 0.48, scale * 0.15, rf"$\alpha={ALPHA_DEG:g}^\circ$", color=INK)
    ax.text(scale * 0.52, -scale * 0.12, "$v_{0x}$", color=GREEN)
    ax.text(-scale * 0.12, scale * 0.55, "$v_{0y}$", color=ORANGE)
    ax.set(xlabel="Położenie x [m]", ylabel="Wysokość y [m]",
           title="Rzut ukośny: prędkość początkowa i jej składowe")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(loc="upper right")
    save(fig, "schemat")


def make_trajectories():
    angle = np.deg2rad(ALPHA_DEG)
    tf = 2 * V0 * np.sin(angle) / G
    t_vacuum = np.linspace(0, tf, 1000)
    ideal = vacuum(ALPHA_DEG, t_vacuum)
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    ax.plot(ideal[0], ideal[1], color=BLUE, label="Bez oporu")
    for model, color, label in [
        ("linear", ORANGE, "Opór liniowy"),
        ("quadratic", GREEN, "Opór kwadratowy"),
    ]:
        _, states = integrate_drag(ALPHA_DEG, model)
        ax.plot(states[0], states[1], color=color, label=label)
    ax.axhline(0, color=INK, linewidth=0.8)
    ax.set(xlabel="Położenie poziome x [m]", ylabel="Wysokość y [m]",
           title=rf"Trajektorie dla $v_0={V0:g}\,\mathrm{{m/s}}$, "
                 rf"$\alpha={ALPHA_DEG:g}^\circ$")
    ax.legend()
    save(fig, "tory")


def make_velocities():
    angle = np.deg2rad(ALPHA_DEG)
    tf = 2 * V0 * np.sin(angle) / G
    t_vacuum = np.linspace(0, tf, 1000)
    ideal = vacuum(ALPHA_DEG, t_vacuum)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True,
                             layout="constrained")
    series = [(t_vacuum, ideal, "Bez oporu", BLUE)]
    for model, label, color in [
        ("linear", "Opór liniowy", ORANGE),
        ("quadratic", "Opór kwadratowy", GREEN),
    ]:
        times, states = integrate_drag(ALPHA_DEG, model)
        series.append((times, states, label, color))
    for ax, (times, states, label, color) in zip(axes, series):
        ax.plot(times, states[2], color=color, label="$v_x$")
        ax.plot(times, states[3], color=INK, linestyle="--", label="$v_y$")
        ax.axhline(0, color=INK, linewidth=0.7, alpha=0.5)
        ax.set(title=label, xlabel="Czas t [s]")
    axes[0].set_ylabel("Prędkość [m/s]")
    axes[-1].legend(fontsize=9)
    fig.suptitle("Składowe prędkości w czasie", fontweight="bold")
    save(fig, "predkosci")


def make_range_plot():
    angles = np.linspace(2, 88, 44)
    ranges = {"vacuum": [], "linear": [], "quadratic": []}
    for angle in angles:
        for model in ranges:
            _, distance = flight_summary(angle, model)
            ranges[model].append(distance)
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    for model, color, label in [
        ("vacuum", BLUE, "Bez oporu"),
        ("linear", ORANGE, "Opór liniowy"),
        ("quadratic", GREEN, "Opór kwadratowy"),
    ]:
        ax.plot(angles, ranges[model], color=color, label=label)
        peak = int(np.argmax(ranges[model]))
        ax.scatter(angles[peak], ranges[model][peak], color=color, zorder=3)
    ax.axvline(45, color=INK, linestyle=":", linewidth=1,
               label="45° w próżni")
    ax.set(xlabel="Kąt wyrzutu α [°]", ylabel="Zasięg R [m]",
           title="Zasięg w zależności od kąta wyrzutu")
    ax.legend()
    save(fig, "zasieg_kat")


def make_energy_plot():
    angle = np.deg2rad(ALPHA_DEG)
    tf = 2 * V0 * np.sin(angle) / G
    t_vacuum = np.linspace(0, tf, 1000)
    ideal = vacuum(ALPHA_DEG, t_vacuum)
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    for times, states, color, label in [
        (t_vacuum, ideal, BLUE, "Bez oporu"),
        (*integrate_drag(ALPHA_DEG, "linear"), ORANGE, "Opór liniowy"),
        (*integrate_drag(ALPHA_DEG, "quadratic"), GREEN, "Opór kwadratowy"),
    ]:
        energy = 0.5 * MASS * (states[2]**2 + states[3]**2) + MASS * G * states[1]
        ax.plot(times, energy, color=color, label=label)
        if label != "Bez oporu":
            assert np.max(np.diff(energy)) < 1e-7
    ax.set(xlabel="Czas t [s]", ylabel="Energia mechaniczna E [J]",
           title="Energia mechaniczna podczas lotu")
    ax.legend()
    save(fig, "energia")


def main():
    if not (V0 > 0 and 0 < ALPHA_DEG < 90):
        raise ValueError("Wymagane są V0 > 0 i 0 < ALPHA_DEG < 90.")
    angle = np.deg2rad(ALPHA_DEG)
    flight_time = 2 * V0 * np.sin(angle) / G
    height = V0**2 * np.sin(angle)**2 / (2 * G)
    distance = V0**2 * np.sin(2 * angle) / G
    print(f"v0 = {V0:.3f} m/s; alpha = {ALPHA_DEG:.3f} deg")
    print(f"Bez oporu: czas = {flight_time:.3f} s; "
          f"wysokość = {height:.3f} m; zasięg = {distance:.3f} m")
    for model in ("linear", "quadratic"):
        time_drag, range_drag = flight_summary(ALPHA_DEG, model)
        print(f"{model}: czas = {time_drag:.3f} s; zasięg = {range_drag:.3f} m")
    make_schematic()
    make_trajectories()
    make_velocities()
    make_range_plot()
    make_energy_plot()
    print(f"Zapisano pięć rysunków w {OUT}")


if __name__ == "__main__":
    main()