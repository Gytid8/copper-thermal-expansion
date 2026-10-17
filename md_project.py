"""Stages 1-3: environment, relaxed FCC Cu, and real NPT pilot runs.

Run with this project's .venv Python. All LAMMPS paths are relative ASCII
paths to avoid Windows C-library filename encoding problems.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))
os.chdir(ROOT)
CONFIG = json.loads((ROOT / "config.json").read_text(encoding="utf-8"))
# Microsoft redistributable extracted locally; no service/system installation.
MPI_DLL_HANDLE = (os.add_dll_directory(str(ROOT / "runtime" / "msmpi"))
                  if sys.platform == "win32" and (ROOT / "runtime" / "msmpi").is_dir() else None)


def save_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def new_lammps(log, threads=1):
    from lammps import lammps
    return lammps(cmdargs=["-screen", "none", "-log", str(log),
                          "-sf", "omp", "-pk", "omp", str(threads)])


def command(lmp, text):
    lmp.commands_string(text)


def initialization():
    return """units metal
dimension 3
boundary p p p
atom_style atomic
"""


def potential():
    return f"""mass 1 {CONFIG['mass_amu']}
pair_style eam/alloy
pair_coeff * * {CONFIG['potential']} Cu
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
"""


def environment():
    Path("outputs").mkdir(exist_ok=True)
    lmp = new_lammps("outputs/environment.log")
    try:
        info = {
            "python": sys.version, "executable": sys.executable,
            "platform": platform.platform(), "logical_cpus": os.cpu_count(),
            "lammps_version": lmp.version(),
            "eam_alloy_available": lmp.has_style("pair", "eam/alloy"),
            "openmp_available": lmp.has_package("OPENMP"),
            "potential_sha256": hashlib.sha256(Path(CONFIG['potential']).read_bytes()).hexdigest(),
            "packages": {p: importlib.metadata.version(p) for p in
                         ["lammps", "numpy", "pandas", "matplotlib", "nbformat", "ipykernel"]},
        }
        if not info["eam_alloy_available"]:
            raise RuntimeError("This LAMMPS build lacks eam/alloy")
        save_json("outputs/environment.json", info)
        print(json.dumps(info, ensure_ascii=False, indent=2), flush=True)
    finally:
        lmp.close()


def build():
    Path("outputs/model").mkdir(parents=True, exist_ok=True)
    n = CONFIG["cells"]
    text = initialization() + f"""lattice fcc {CONFIG['initial_lattice_A']}
region box block 0 {n} 0 {n} 0 {n} units lattice
create_box 1 box
create_atoms 1 box
""" + potential() + """thermo 100
thermo_style custom step atoms pe press lx ly lz vol fmax
thermo_modify format float %.12g
write_data outputs/model/initial.data
fix relax all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1.0e-12 1.0e-10 10000 100000
unfix relax
compute structure all cna/atom 3.2
variable isfcc atom c_structure==1
compute nfcc all reduce sum v_isfcc
thermo_style custom step atoms pe press lx ly lz vol fmax c_nfcc
run 0
write_data outputs/model/relaxed.data
write_dump all custom outputs/model/relaxed.lammpstrj id type x y z c_structure modify sort id
"""
    Path("outputs/model/build.in").write_text(text, encoding="ascii")
    lmp = new_lammps("outputs/model/log.lammps")
    try:
        lmp.file("outputs/model/build.in")
        info = {k: lmp.get_thermo(k) for k in ["atoms", "pe", "press", "lx", "ly", "lz", "vol", "fmax"]}
        info["fcc_atoms"] = int(round(lmp.extract_compute("nfcc", 0, 0)))
        info["lattice_A"] = info["lx"] / n
        info["energy_per_atom_eV"] = info["pe"] / info["atoms"]
        info["cna_cutoff_A"] = 3.2
        assert int(info["atoms"]) == 4 * n**3
        assert info["fcc_atoms"] == int(info["atoms"])
        assert abs(info["lx"] - info["ly"]) < 1e-7
        assert abs(info["lx"] - info["lz"]) < 1e-7
        assert abs(info["press"]) < 1.0, info
        assert 3.4 < info["lattice_A"] < 3.8
        save_json("outputs/model/model_check.json", info)
        print("MODEL", json.dumps(info), flush=True)
    finally:
        lmp.close()
    plot_model()


def plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif":
                         ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False,
                         "figure.dpi": 130, "savefig.dpi": 180})
    return plt


def plot_model():
    import numpy as np
    plt = plotting()
    data = np.loadtxt("outputs/model/relaxed.lammpstrj", skiprows=9)
    info = json.loads(Path("outputs/model/model_check.json").read_text())
    a = info["lattice_A"]
    fig = plt.figure(figsize=(11, 5))
    for index, (mask, title) in enumerate([
        (np.ones(len(data), dtype=bool), "4000原子周期性FCC单晶铜"),
        ((data[:, 2:5] <= a + 1e-6).all(axis=1), "局部常规晶胞（含边界原子）"),
    ], start=1):
        ax = fig.add_subplot(1, 2, index, projection="3d")
        points = data[mask, 2:5]
        ax.scatter(*points.T, c="#c97736", s=8 if index == 1 else 100,
                   alpha=0.38 if index == 1 else 1.0, edgecolors="none")
        ax.set(xlabel="x (Å)", ylabel="y (Å)", zlabel="z (Å)", title=title)
        ax.set_box_aspect([1, 1, 1])
    fig.suptitle(f"静态弛豫后：a = {a:.6f} Å；CNA识别4000/4000为FCC")
    fig.tight_layout()
    fig.savefig("outputs/model/model_preview.png", bbox_inches="tight")
    plt.close(fig)


def input_trial(temperature, folder, equil_ps, sample_ps, seed):
    c = CONFIG
    dt = c["timestep_ps"]
    eq = round(equil_ps / dt)
    count = round(sample_ps / dt)
    every = round(c["sample_interval_ps"] / dt)
    dump_every = round(c["dump_interval_ps"] / dt)
    return initialization() + "read_data outputs/model/relaxed.data\n" + potential() + f"""timestep {dt}
reset_timestep 0
velocity all create {temperature} {seed} mom yes rot no dist gaussian
fix integrator all npt temp {temperature} {temperature} {c['temperature_damping_ps']} iso 0 0 {c['pressure_damping_ps']}
variable time_ps equal time
variable temperature_K equal temp
variable pressure_bar equal press
variable energy_eV equal pe
variable length_x equal lx
variable length_y equal ly
variable length_z equal lz
variable volume_A3 equal vol
thermo 1000
thermo_style custom step time atoms temp press pe etotal lx ly lz vol
thermo_modify format float %.12g flush yes
fix scalar all ave/time {every} 1 {every} v_time_ps v_temperature_K v_pressure_bar v_energy_eV v_length_x v_length_y v_length_z v_volume_A3 file {folder}/thermo.dat
dump trajectory all custom {dump_every} {folder}/trajectory.lammpstrj id type x y z
dump_modify trajectory sort id
run {eq}
write_restart {folder}/equilibrated.restart
run {count}
write_data {folder}/final.data
write_restart {folder}/final.restart
"""


def trial(temperature, threads, benchmark=False):
    if not Path("outputs/model/relaxed.data").is_file():
        raise RuntimeError("Run build first")
    seed = CONFIG["trial_seed"]
    folder = Path("outputs") / (f"benchmark_t{threads}" if benchmark else "trials")
    if not benchmark:
        folder = folder / f"T{temperature}_seed{seed}"
    if (folder / "metadata.json").exists():
        raise FileExistsError(f"Preserving existing run: {folder}")
    folder.mkdir(parents=True, exist_ok=True)
    eq, sample = ((1, 2) if benchmark else
                  (CONFIG["trial_equilibration_ps"], CONFIG["trial_sampling_ps"]))
    text = input_trial(temperature, folder.as_posix(), eq, sample, seed)
    input_file = folder / "input.in"
    input_file.write_text(text, encoding="ascii")
    metadata = {"status": "running", "temperature_K": temperature, "seed": seed,
                "threads": threads, "equilibration_ps": eq, "sampling_ps": sample,
                "config": CONFIG, "started_local": time.strftime("%Y-%m-%d %H:%M:%S")}
    save_json(folder / "metadata.json", metadata)
    start = time.perf_counter()
    lmp = new_lammps((folder / "log.lammps").as_posix(), threads)
    try:
        lmp.file(input_file.as_posix())
        assert int(lmp.get_natoms()) == 4 * CONFIG["cells"]**3
        metadata["status"] = "completed"
        metadata["final_atoms"] = int(lmp.get_natoms())
    except Exception as exc:
        metadata.update(status="failed", error=str(exc))
        raise
    finally:
        lmp.close()
        metadata["elapsed_seconds"] = time.perf_counter() - start
        metadata["steps_per_second"] = round((eq + sample) / CONFIG["timestep_ps"] / metadata["elapsed_seconds"], 2)
        save_json(folder / "metadata.json", metadata)
    print("COMPLETED", folder, metadata["elapsed_seconds"], "seconds", flush=True)


def analyze():
    import numpy as np
    import pandas as pd
    plt = plotting()
    columns = ["step", "time_ps", "temperature_K", "pressure_bar", "potential_eV",
               "lx_A", "ly_A", "lz_A", "volume_A3"]
    rows = []
    for folder in sorted(Path("outputs/trials").glob("T*_seed*")):
        meta = json.loads((folder / "metadata.json").read_text())
        if meta["status"] != "completed":
            continue
        frame = pd.read_csv(folder / "thermo.dat", sep=r"\s+", comment="#", names=columns)
        assert np.isfinite(frame.to_numpy()).all()
        # ave/time with one sample also writes the initial step zero.
        expected_intervals = round((meta["equilibration_ps"] + meta["sampling_ps"]) / CONFIG["sample_interval_ps"])
        assert len(frame) == expected_intervals + 1
        step_interval = round(CONFIG["sample_interval_ps"] / CONFIG["timestep_ps"])
        assert np.array_equal(frame.step.to_numpy(), np.arange(expected_intervals + 1) * step_interval)
        production = frame[frame.time_ps > meta["equilibration_ps"]].copy()
        assert len(production) == round(meta["sampling_ps"] / CONFIG["sample_interval_ps"])
        frame.to_csv(folder / "thermo.csv", index=False)
        a = production.lx_A.to_numpy() / CONFIG["cells"]
        n_half = len(a) // 2
        row = {"temperature_target_K": meta["temperature_K"], "seed": meta["seed"],
               "samples": len(a), "temperature_mean_K": production.temperature_K.mean(),
               "pressure_mean_bar": production.pressure_bar.mean(),
               "lattice_mean_A": a.mean(), "volume_mean_A3": production.volume_A3.mean(),
               "first_half_a_A": a[:n_half].mean(), "second_half_a_A": a[n_half:].mean(),
               "halves_difference_A": a[n_half:].mean() - a[:n_half].mean(),
               "wall_seconds": meta["elapsed_seconds"],
               "status": "pilot_only_not_final_convergence"}
        blockrows = []
        for block_ps in [1, 2, 5, 10, 20]:
            blocklen = round(block_ps / CONFIG["sample_interval_ps"])
            nblocks = len(a) // blocklen
            means = a[:nblocks * blocklen].reshape(nblocks, blocklen).mean(axis=1)
            sem = means.std(ddof=1) / np.sqrt(nblocks)
            blockrows.append({"block_ps": block_ps, "n_blocks": nblocks, "a_sem_A": sem})
            if block_ps == 10:
                row["a_sem_10ps_blocks_A"] = sem
                half_block_count = len(means) // 2
                row["halves_difference_sem_A"] = float(np.sqrt(
                    means[:half_block_count].var(ddof=1) / half_block_count +
                    means[half_block_count:].var(ddof=1) / (len(means) - half_block_count)))
                for key in ["temperature_K", "pressure_bar"]:
                    values = production[key].to_numpy()
                    bmeans = values[:nblocks * blocklen].reshape(nblocks, blocklen).mean(axis=1)
                    row[key + "_sem_10ps"] = bmeans.std(ddof=1) / np.sqrt(nblocks)
        pd.DataFrame(blockrows).to_csv(folder / "block_statistics.csv", index=False)
        fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)
        for ax, key, label, scale in zip(axes,
            ["temperature_K", "pressure_bar", "potential_eV", "volume_A3"],
            ["温度 (K)", "压力 (bar)", "势能 (eV/atom)", "体积 (Å³)"],
            [1, 1, 1 / 4000, 1]):
            ax.plot(frame.time_ps, frame[key] * scale, color="#b6c7d1", lw=0.45, label="原始记录")
            ax.plot(frame.time_ps, (frame[key] * scale).rolling(50, min_periods=1).mean(),
                    color="#b7662d", lw=1, label="5 ps滑动平均（仅展示）")
            ax.axvline(meta["equilibration_ps"], color="#444444", ls="--", lw=1)
            ax.axvspan(0, meta["equilibration_ps"], color="#dddddd", alpha=.25)
            ax.set_ylabel(label)
            ax.grid(alpha=.2)
        axes[0].axhline(meta["temperature_K"], color="black", ls=":", lw=.8)
        axes[1].axhline(0, color="black", ls=":", lw=.8)
        axes[0].legend(loc="upper right", fontsize=8)
        axes[-1].set_xlabel("时间 (ps)；虚线后为采样阶段")
        fig.suptitle(f"{meta['temperature_K']} K试跑：100 ps平衡＋200 ps采样；单一种子")
        fig.tight_layout()
        fig.savefig(folder / "equilibration.png", bbox_inches="tight")
        plt.close(fig)
        rows.append(row)
    if len(rows) != 2:
        raise RuntimeError(f"Expected two completed pilots, found {len(rows)}")
    result = pd.DataFrame(rows).sort_values("temperature_target_K")
    result.to_csv("outputs/trials/pilot_summary.csv", index=False)
    print(result.to_string(index=False), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["environment", "build", "benchmark", "trial", "analyze"])
    parser.add_argument("--temperature", type=int, choices=[300, 800], default=300)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.action == "environment": environment()
    elif args.action == "build": build()
    elif args.action == "benchmark": trial(args.temperature, args.threads, benchmark=True)
    elif args.action == "trial": trial(args.temperature, args.threads)
    else: analyze()
