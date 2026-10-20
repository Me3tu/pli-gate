#!/usr/bin/env python3
"""O4 re-run driver, corrected — sequential batched GATE for the provenance re-runs.

WHY THIS EXISTS (read before changing anything)
-----------------------------------------------
The 2026-08-13 overnight attempt (`overnight_driver.sh`) ran each phantom as a
SINGLE batch of 16 shards with a 5 GB per-shard cgroup cap. That is 80 GB of
cap on a 78 GB machine. GATE's per-shard RSS grows with wall time, so a long
batch walks the whole machine into the wall: at 2026-08-13 23:37 the kernel OOM
killer fired system-wide, destroyed the user session (dbus, pipewire,
gnome-keyring) and with it the driver, ~9 h of the block, and a 12 GB partial
Sim-1 run at 27 % completion.

Three separate faults are corrected here:

  1. BATCHING. Every phantom now runs as N sequential batches sized to the
     wall-time envelope its ORIGINAL production run used. Peak memory is
     16 x 4 GB = 64 GB, the same envelope as the VIP-Man 1.4B production run
     (32 batches x ~32 min) and the 2026-08-12 water re-run (8 batches x 12
     min). This is the pattern the project's own scaling lesson mandates.

  2. MODE. `overnight_driver.sh` invoked steinberger_disks and brain_tumor in
     --total-primaries mode. Neither macro carries a
     /gate/application/setTotalNumberOfPrimaries line, so the runner would have
     hard-exited (rc 5) on both. Both originals ran in time mode. The recipes
     below are transcribed from the original run.log of each production run,
     not invented.

  3. DELETION GATE. ROOT was deleted after checking only that the triples h5
     existed. Sim-2 extracted zero rows (Sc-44 extractor on an I-124 phantom)
     and its ROOT was deleted anyway. Deletion is now gated on a row count.

Not handled here: VIP-Man, which has its own validated per-batch driver
(reconstruction/utils/run_phase3q_per_batch.py, 32 batches). Sim-2 is handled
by this driver as well (phantom `sim2_mercolli`); o4_sim2.sh is superseded.

Resume: a batch whose triples h5 already has rows is skipped, so the driver can
be restarted at any point without redoing work.
"""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import h5py
import numpy as np

REPO = Path(os.environ.get("PLI_REPO") or Path(__file__).resolve().parents[3])
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
OUT = RES / "2026-08-13_overnight"
RERUN = OUT / "O4_rerun"
RUNNER = REPO / "scripts/run_gate_parallel_musire.sh"
LOG = OUT / "o4_batched.log"

THREADS = 16
CAP_GB = 4            # 16 x 4 = 64 GB peak, the production-proven envelope
DISK_FLOOR_GB = 40
MEM_START_FLOOR_GB = 25   # refuse to start a batch without this much free
MEM_ABORT_FLOOR_GB = 12   # kill a RUNNING batch below this (see watchdog)
WATCHDOG_PERIOD_S = 15

P2B = RES / "2026-04-29_2320_phase2b_validation_sims/macros"
STB = RES / "2026-06-19_steinberger_validation/macros"

# Each recipe's `mode`, `total` and `n_batches` reproduce the ORIGINAL run's
# per-shard workload, split into batches short enough to stay inside the memory
# envelope. `per_batch` is what is passed to --total-primaries / --total-time-sec.
PHANTOMS = [
    # phase2b Sim-2 arm — cheapest run in the set (~12 min) and the one with a
    # matched pre-correction counterpart (200 M, 2026-04-30) for the yield test.
    dict(name="sim2_mercolli",
         macro=P2B / "sim2_i124_xad4_phantom.mac",
         mode="primaries", total=50_000_000, n_batches=1,
         extractor="extract_triples_i124.py", pli_log=1,
         note="I-124: MUST use extract_triples_i124.py (511: 460-545, prompt: "
              "568-639 keV). The Sc-44 extractor's [1100,1200] keV prompt window "
              "returns zero triples on this phantom."),
    # Sim-2, the arm that actually backs the manuscript. O1's provenance table
    # points the Sim-2 closure at the ARCHIVED `mercolli_te_v2` run, which used
    # mercolli_te_poisson.mac in time mode -- NOT the phase2b
    # sim2_i124_xad4_phantom.mac that o4_sim2.sh re-ran (that macro was
    # inherited from overnight_driver.sh without being checked against the
    # original). Recipe transcribed from the archived sim/batch_*/run.log:
    # 4 batches x 2.5201612903225805 s, 16 threads, 4 GB.
    dict(name="sim2_mercolli_te",
         macro=OUT / "macros/mercolli_te_poisson.mac",
         mode="time", total=10.080645161290322, n_batches=4,
         extractor="extract_triples_i124.py", pli_log=1,
         note="like-for-like supersede of the archived mercolli_te_v2 TE arm."),
    dict(name="sim1_calibration",
         macro=P2B / "sim1_sc44_validation_phantom.mac",
         mode="primaries", total=250_000_000, n_batches=10,
         extractor="extract_triples_chunked.py", pli_log=1,
         note="original: 1 x 250M primaries. Measured 2,660 primaries/s/shard "
              "=> 25M/batch is ~10 min/batch."),
    dict(name="steinberger_disks",
         macro=STB / "steinberger_disk_phantom.mac",
         mode="time", total=60.0, n_batches=12,
         extractor="extract_triples_chunked.py", pli_log=0, singles=True,
         note="original 2026-06-22: time mode, 60 s total, 16 threads, 3 GB cap "
              "-- and it FAILED (shard 2 killed at the cap; only 14/16 shards "
              "survived; see findings). MEASURED 2026-08-14: this phantom is by "
              "far the most memory-hungry per unit of simulated time -- dense "
              "glass at high activity, ~28 min of shard wall per simulated "
              "second, 0.156 GB/min/shard. A 4-batch (15 s) attempt ran 25.8 min "
              "and died at the 4 GB cap having taken the machine to 67.5 GB with "
              "10.8 GB left. 12 batches x 5 s keeps each shard ~8.6 min and peak "
              "RSS ~25 GB."),
    dict(name="steinberger_water",
         macro=RES / "2026-08-12_steinberger_allprompts/water_7MBq.mac",
         mode="time", total=23.44, n_batches=8,
         extractor="extract_triples_chunked.py", pli_log=0, singles=True,
         note="matches the 2026-08-12 S3 re-run exactly: 8 batches x 2.93 s."),
    dict(name="brain_tumor",
         macro=OUT / "macros/brain_tumor_te_poisson_1B_truetime.mac",
         mode="time", total=979.2, n_batches=16,
         extractor="extract_triples_chunked.py", pli_log=0, singles=True,
         note="matches the 2026-06-10 original: 16 batches x 61.2 s, ~16 min each. "
              "BATCH COUNT CORRECTED: run_summary.json records batch_ids 0-15, i.e. 16 "
              "batches. The earlier recipe said 17, taken from `ls -d batch_logs/*/ | wc -l` "
              "which over-counts by one because ls is aliased to a table-printing tool. "
              "total is 16 x 61.2 = 979.2 s so per_batch stays 61.2 s. "
              "singles=True: M4_recovery_curve and M3_clinical_extractor_randoms both read "
              "--singles, and the first pass kept only triples. "
              "PATCHED COPY: the original macro's line 60 /control/execute pointed at "
              "2026-05-08_0010_phase3s3_brain_tumor_1B/brain_tumor_sources_poisson.mac, "
              "a directory ARCHIVED to /data on 2026-08-12 — so GATE aborted instantly "
              "with 'Can not open a macro file'. Re-pointed at the byte-identical "
              "surviving copy in 2026-06-10_brain_tumor_1B_truetime/ (both sha256 "
              "c29860a704bad160...), so the physics is unchanged."),
]


def say(msg: str):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    with open(LOG, "a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def free_gb() -> float:
    st = os.statvfs("/home")
    return st.f_bavail * st.f_frsize / 1e9


def mem_available_gb() -> float:
    with open("/proc/meminfo") as fh:
        for ln in fh:
            if ln.startswith("MemAvailable"):
                return int(ln.split()[1]) / 1048576
    return 0.0


def rows(h5path: Path) -> int:
    """Row count of the first dataset under /detector; 0 if missing/unreadable."""
    try:
        with h5py.File(h5path, "r") as h:
            ks = list(h["detector"].keys())
            return int(h["detector"][ks[0]].shape[0]) if ks else 0
    except Exception:
        return 0


def kill_gate() -> int:
    """Kill every running Gate shard. Matches on `comm`, never on the full
    command line -- a pattern matched against argv can also match this process
    or the shell that launched it, which is how an earlier script in this block
    killed itself."""
    out = subprocess.run(["ps", "-eo", "pid,comm", "--no-headers"],
                         capture_output=True, text=True).stdout
    n = 0
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1] == "Gate":
            try:
                os.kill(int(parts[0]), 15)
                n += 1
            except (ProcessLookupError, PermissionError):
                pass
    return n


class MemoryWatchdog(threading.Thread):
    """Abort a running batch if free memory collapses.

    The pre-batch MemAvailable check is not sufficient on its own: GATE's RSS
    grows with wall time, so a batch that starts with 75 GB free can still walk
    the machine to the wall. The 2026-08-13 OOM did exactly that, and on
    2026-08-14 the steinberger_disks batch reached 67.5 GB with 10.8 GB left
    before its shards hit their cgroup caps. Per-shard cgroup caps bound each
    process but 16 x 4 GB still exceeds what is safely available, so a
    system-level guard is required as well.

    Killing the shards makes `parallel --halt now,fail=1` fail, the runner
    returns non-zero, and run_batch() records a clean batch failure.
    """

    def __init__(self, label: str):
        super().__init__(daemon=True)
        self.label = label
        self.stop = threading.Event()
        self.tripped = False
        self.low_water = float("inf")

    def run(self):
        strikes = 0
        while not self.stop.wait(WATCHDOG_PERIOD_S):
            avail = mem_available_gb()
            self.low_water = min(self.low_water, avail)
            if avail < MEM_ABORT_FLOOR_GB:
                strikes += 1
                say(f"WATCHDOG {self.label}: MemAvailable {avail:.1f} GB "
                    f"< {MEM_ABORT_FLOOR_GB} GB (strike {strikes}/2)")
                if strikes >= 2:
                    self.tripped = True
                    n = kill_gate()
                    say(f"WATCHDOG {self.label}: ABORTED batch, killed {n} Gate shard(s)")
                    return
            else:
                strikes = 0


# ---------------------------------------------------------------------------
# PRODUCT REGISTRY
#
# Which product a phantom must retain is a property of the ANALYSES it feeds,
# not of the phantom. Nothing in this driver recorded that mapping, so it was
# rediscovered the hard way: the Steinberger and brain-tumor analyses all take
# `--singles`, but the driver deletes the ROOT once triples are extracted. That
# silently destroyed the only usable input for three phantoms — caught for disks
# and water, missed for brain-tumor and Sim-1, each costing a full re-run.
#
# The mapping is now explicit and checked at preflight. A phantom whose recipe
# cannot feed an analysis it is registered for does not launch.
# ---------------------------------------------------------------------------
ANALYSIS_REQUIRES = {
    "analyze_per_tube_tau":           {"triples"},   # per-tube closure
    "filter_triples_tof_pann":        {"triples"},   # locked-recipe r-cut
    "correct_travel_time":            {"triples"},
    "simple_recon":                   {"triples"},
    "steinberger_extract_ap":         {"singles"},   # VOI oracle + clinical
    "M3_tubes":                       {"singles"},   # per-ROI clinical, randoms
    "M3_clinical_extractor_randoms":  {"singles"},
    "M4_crt_sweep":                   {"singles"},
    "M4_recovery_curve":              {"singles"},   # r_min, detection sigma
    "action_B_range_share":           {"oracle_log"},# migration budget, boundary frac
    "score_blind_products":           {"oracle_log"},# truth-scored blind evaluation
}

PHANTOM_ANALYSES = {
    "sim2_mercolli":     ["analyze_per_tube_tau", "filter_triples_tof_pann"],
    "sim2_mercolli_te":  ["analyze_per_tube_tau", "filter_triples_tof_pann"],
    # Sim-1 also feeds the M3 clinical per-ROI table (detection sigma 4.2-11.5)
    # and the M4 CRT sweep, both of which read SINGLES.
    "sim1_calibration":  ["analyze_per_tube_tau", "filter_triples_tof_pann",
                          "M3_tubes", "M4_crt_sweep"],
    "steinberger_disks": ["steinberger_extract_ap"],
    "steinberger_water": ["steinberger_extract_ap"],
    # brain-tumor feeds r_min / detection sigma via the recovery curve.
    "brain_tumor":       ["M4_recovery_curve", "M3_clinical_extractor_randoms"],
    # VIP-Man also feeds the migration budget and the truth-scored blind
    # evaluation, both of which read the PLI_LIFETIME oracle log. Its recipe
    # CANNOT currently satisfy that -- run_phase3q_per_batch's cleanup_batch
    # deletes the whole sim/ directory -- so this registration deliberately
    # BLOCKS vipman until the recipe retains the log. That is the point: the
    # log was being deleted every run because nothing claimed it.
    "vipman":            ["correct_travel_time", "simple_recon",
                          "action_B_range_share", "score_blind_products"],
}


def products_of(ph: dict) -> set[str]:
    """What a recipe actually retains after cleanup."""
    p = {"triples"}
    if ph.get("singles"):
        p.add("singles")
    if ph.get("pli_log"):
        p.add("oracle_log")     # NOTE: shard logs are deleted, so this is
                                # retained only if cleanup is changed too.
    return p


def product_gap(ph: dict) -> list[str]:
    """Analyses this phantom is registered for that its recipe cannot feed."""
    have = products_of(ph)
    gaps = []
    for a in PHANTOM_ANALYSES.get(ph["name"], []):
        need = ANALYSIS_REQUIRES.get(a, set())
        missing = need - have
        if missing:
            gaps.append(f"{a} needs {sorted(missing)}")
    return gaps


RE_ABSPATH = __import__("re").compile(r"/home/[^\s\"']+")
RE_EXECUTE = __import__("re").compile(r"^/control/execute\s+(\S+)", __import__("re").M)


def preflight(ph: dict) -> list[str]:
    """Check every absolute path a macro references, following /control/execute.

    Added after brain_tumor burned 17 batches x 2 attempts in 19 seconds: its
    macro executed a sub-macro under a directory that had been archived to
    /data, so GATE aborted at once with 'Can not open a macro file'. The runner
    reports that as a generic rc=6, indistinguishable from a resource failure,
    and the driver dutifully retried it 34 times. One second of checking would
    have caught it.
    """
    missing, seen = [], set()
    stack = [Path(ph["macro"])]
    while stack:
        mac = stack.pop()
        if mac in seen:
            continue
        seen.add(mac)
        if not mac.exists():
            missing.append(str(mac))
            continue
        txt = mac.read_text()
        for p in RE_ABSPATH.findall(txt):
            if not os.path.exists(p):
                missing.append(p)
        for inc in RE_EXECUTE.findall(txt):
            stack.append(Path(inc) if Path(inc).is_absolute() else mac.parent / inc)
    return sorted(set(missing))


def run_batch(ph: dict, batch_id: int) -> bool:
    """One batch: GATE -> extract -> row gate -> delete ROOT. True on success."""
    name = ph["name"]
    bdir = RERUN / name / f"batch_{batch_id:02d}"
    det = bdir / f"triples_detector_{name}.h5"

    n_existing = rows(det)
    if n_existing > 0:
        say(f"{name} batch {batch_id}: already done ({n_existing:,} triples), skipping")
        return True

    if free_gb() < DISK_FLOOR_GB:
        say(f"{name} batch {batch_id}: ABORT, free {free_gb():.0f} GB < {DISK_FLOOR_GB} GB floor")
        return False
    if mem_available_gb() < MEM_START_FLOOR_GB:
        say(f"{name} batch {batch_id}: ABORT, MemAvailable {mem_available_gb():.0f} GB "
            f"< {MEM_START_FLOOR_GB} GB start floor")
        return False

    for attempt in (1, 2):
        if bdir.exists():
            # Keep the failed attempt's logs. The first version of this driver
            # deleted them on retry, which destroyed the diagnostics for the
            # steinberger_disks failure and forced the cause to be recovered
            # from the memory sampler instead.
            failed = bdir.with_name(bdir.name + f"_failed_attempt{attempt-1}")
            shutil.rmtree(failed, ignore_errors=True)
            for f in bdir.glob("sim-0*.root"):
                f.unlink()
            for f in bdir.glob("sim-merged.root"):
                f.unlink()
            bdir.rename(failed)
        bdir.mkdir(parents=True)
        per_batch = ph["total"] / ph["n_batches"]
        flag = "--total-primaries" if ph["mode"] == "primaries" else "--total-time-sec"
        amount = str(int(per_batch)) if ph["mode"] == "primaries" else f"{per_batch:.6f}"
        cmd = [str(RUNNER), "--macro", str(ph["macro"]), flag, amount,
               "--threads", str(THREADS), "--bounded-max-gb", str(CAP_GB),
               "--batch-id", str(batch_id), "--label", f"o4_{name}",
               "--output-dir", str(bdir)]
        env = os.environ.copy()
        env["PLI_LIFETIME_LOG"] = str(ph["pli_log"])
        say(f"{name} batch {batch_id} attempt {attempt}: {flag} {amount} "
            f"(free {free_gb():.0f} GB, mem {mem_available_gb():.0f} GB)")
        t0 = time.time()
        dog = MemoryWatchdog(f"{name}/batch{batch_id}")
        dog.start()
        try:
            with open(bdir / "run.log", "a") as rl:
                rc = subprocess.call(cmd, cwd=REPO, env=env, stdout=rl, stderr=rl)
        finally:
            dog.stop.set()
            dog.join(timeout=WATCHDOG_PERIOD_S + 5)
        merged = bdir / "sim-merged.root"
        if dog.tripped:
            say(f"{name} batch {batch_id} attempt {attempt} ABORTED BY WATCHDOG "
                f"(low water {dog.low_water:.1f} GB) — batch is too large for this "
                f"machine; reduce per-batch simulated time, do not retry as-is")
            return False
        if rc != 0 or not merged.exists():
            say(f"{name} batch {batch_id} attempt {attempt} SIM FAILED rc={rc} "
                f"(mem low water {dog.low_water:.1f} GB)")
            continue
        say(f"{name} batch {batch_id} SIM OK wall={time.time()-t0:.0f}s "
            f"mem low water {dog.low_water:.1f} GB")

        # extractors glob sim_*.root (underscore); the runner writes sim-merged.root
        (bdir / "sim_merged.root").symlink_to("sim-merged.root")
        with open(bdir / "extract.log", "a") as el:
            subprocess.call([sys.executable, f"reconstruction/{ph['extractor']}",
                             "--root-dir", str(bdir), "--out-dir", str(bdir),
                             "--tag", f"_{name}"], cwd=REPO, env=env,
                            stdout=el, stderr=el)
        n = rows(det)
        say(f"{name} batch {batch_id}: extracted {n:,} triples")
        if n <= 0:
            say(f"{name} batch {batch_id}: EXTRACTION EMPTY — ROOT KEPT for diagnosis")
            return False

        # Some phantoms are analysed from SINGLES, not triples. The Steinberger
        # VOI analysis (steinberger_extract_ap.py --singles) is one: it does its
        # own coincidence handling, including the delayed-window randoms branch
        # that the all-prompts numbers depend on. The first pass of this driver
        # deleted the ROOT after extracting triples, which left those two
        # phantoms unable to feed their manuscript numbers at all. Convert here,
        # BEFORE the ROOT goes.
        if ph.get("singles"):
            sing = bdir / f"singles_{name}.h5"
            with open(bdir / "singles.log", "a") as sl:
                rc_s = subprocess.call(
                    [sys.executable,
                     "reconstruction/results/2026-06-10_brain_tumor_1B_truetime/"
                     "scripts/gate_root_to_singles_h5.py",
                     "--root", str(merged), "--out", str(sing),
                     "--batch-id", str(batch_id)],
                    cwd=REPO, env=env, stdout=sl, stderr=sl)
            if rc_s != 0 or not sing.exists():
                say(f"{name} batch {batch_id}: SINGLES CONVERSION FAILED rc={rc_s}"
                    f" — ROOT KEPT")
                return False
            say(f"{name} batch {batch_id}: singles h5 "
                f"{sing.stat().st_size/1e9:.2f} GB")

        # ROOT + per-shard logs are regenerable; the triples h5 are the product.
        for pat in ("sim-0*.root", "sim-merged.root", "sim_merged.root", "sim-0*.log"):
            for f in bdir.glob(pat):
                f.unlink()
        say(f"{name} batch {batch_id} DONE, free={free_gb():.0f} GB")
        return True

    say(f"{name} batch {batch_id} ABANDONED after 2 attempts")
    return False


def concat(ph: dict, batch_ids: list[int]) -> Path | None:
    """Concat per-batch triples into one file per phantom (production pattern).

    Batches are concatenated at the TRIPLES level, never by hadd-ing ROOT:
    each batch restarts eventID at 0, and the extractor groups coincidences by
    eventID, so merging ROOT across batches would manufacture cross-batch
    false triples (the documented eventID-overlap pitfall).
    """
    name = ph["name"]
    pdir = RERUN / name
    out_det = pdir / f"triples_detector_{name}.h5"
    out_tru = pdir / f"triples_truth_{name}.h5"
    dets = [pdir / f"batch_{i:02d}/triples_detector_{name}.h5" for i in batch_ids]
    trus = [pdir / f"batch_{i:02d}/triples_truth_{name}.h5" for i in batch_ids]
    dets = [d for d in dets if rows(d) > 0]
    trus = [t for t in trus if t.exists()]
    if not dets:
        say(f"{name}: nothing to concat")
        return None
    counts = [rows(d) for d in dets]
    n_total = sum(counts)

    with h5py.File(dets[0], "r") as h:
        det_keys = sorted(h["detector"].keys())
        det_attrs = dict(h.attrs)
    with h5py.File(trus[0], "r") as h:
        tru_keys = sorted(h["truth"].keys())
    int64_truth = {"eventID", "source_file_index"}

    with h5py.File(out_det, "w") as fd, h5py.File(out_tru, "w") as ft:
        dg = fd.create_group("detector")
        tg = ft.create_group("truth")
        for k in det_keys:
            dg.create_dataset(k, shape=(n_total,), dtype=np.float64,
                              compression="gzip", compression_opts=4)
        for k in tru_keys:
            tg.create_dataset(k, shape=(n_total,),
                              dtype=np.int64 if k in int64_truth else np.float64,
                              compression="gzip", compression_opts=4)
        cur = 0
        for d, t, bid, n in zip(dets, trus, batch_ids, counts):
            with h5py.File(d, "r") as h:
                for k in det_keys:
                    dg[k][cur:cur + n] = h["detector"][k][:]
            with h5py.File(t, "r") as h:
                for k in tru_keys:
                    if k == "source_file_index":
                        tg[k][cur:cur + n] = bid
                    else:
                        tg[k][cur:cur + n] = h["truth"][k][:]
            cur += n
        for k, v in det_attrs.items():
            if k != "n_triples":
                fd.attrs[k] = v
        fd.attrs["n_triples"] = n_total
        fd.attrs["n_batches"] = len(dets)
        fd.attrs["note"] = f"O4 re-run concat across {len(dets)} batches (frozen 2026-08-13 binary)"
        ft.attrs["n_triples"] = n_total
        ft.attrs["n_batches"] = len(dets)
        ft.attrs["note"] = "O4 re-run truth; source_file_index = batch_id"
    say(f"{name}: concat {len(dets)} batches -> {n_total:,} triples "
        f"({out_det.stat().st_size/1e9:.2f} GB)")
    return out_det


def filter_rcut(ph: dict, det: Path):
    """Locked recipe TOF filter, R_cut = 75 mm."""
    name = ph["name"]
    pdir = RERUN / name
    tru = pdir / f"triples_truth_{name}.h5"
    with open(pdir / "filter.log", "a") as fl:
        rc = subprocess.call([sys.executable, "reconstruction/filter_triples_tof_pann.py",
                              "--triples-detector", str(det),
                              "--triples-truth", str(tru), "--r-cut", "75.0",
                              "--output-detector", str(pdir / f"triples_detector_{name}_rcut75.h5"),
                              "--output-truth", str(pdir / f"triples_truth_{name}_rcut75.h5")],
                             cwd=REPO, stdout=fl, stderr=fl)
    say(f"{name}: r-cut filter rc={rc}, "
        f"rows={rows(pdir / f'triples_detector_{name}_rcut75.h5'):,}")


# Run order. brain_tumor goes first because the two Steinberger phantoms have
# to be re-simulated from scratch (their first pass kept only triples, and the
# VOI analysis needs singles), and everything here must finish before VIP-Man's
# ~18 h run starts.
# R5 order: cheapest first, VIP-Man last (launched by o4_chain.sh after this).
RUN_ORDER = ["sim2_mercolli", "sim2_mercolli_te", "sim1_calibration",
             "steinberger_disks", "steinberger_water", "brain_tumor"]


def main():
    only = sys.argv[1:] or None
    PHANTOMS.sort(key=lambda p: RUN_ORDER.index(p["name"])
                  if p["name"] in RUN_ORDER else 99)
    say("=== O4 batched driver START ===")
    say("run order: " + ", ".join(p["name"] for p in PHANTOMS))
    say(f"threads={THREADS} cap={CAP_GB}GB/shard peak={THREADS*CAP_GB}GB "
        f"(machine 78 GB); disk floor {DISK_FLOOR_GB} GB")
    summary = {}
    for ph in PHANTOMS:
        if only and ph["name"] not in only:
            continue
        say(f"--- {ph['name']}: {ph['n_batches']} batches, {ph['mode']} mode ---")
        say(f"    {ph['note']}")
        gaps = product_gap(ph)
        if gaps:
            say(f"{ph['name']}: PRODUCT-REGISTRY FAILURE — the recipe cannot feed "
                f"{len(gaps)} registered analysis/analyses:")
            for g in gaps:
                say(f"    MISSING PRODUCT: {g}")
            say(f"{ph['name']}: SKIPPED (add the product to the recipe before launching; "
                f"deleting ROOT without it destroys the only usable input)")
            summary[ph["name"]] = dict(batches_ok=0, product_gap=gaps)
            continue
        say(f"{ph['name']}: product registry OK — provides {sorted(products_of(ph))}")
        miss = preflight(ph)
        if miss:
            say(f"{ph['name']}: PREFLIGHT FAILED — {len(miss)} referenced path(s) missing:")
            for p in miss:
                say(f"    MISSING {p}")
            say(f"{ph['name']}: SKIPPED (no batches run; fix the macro paths first)")
            summary[ph["name"]] = dict(batches_ok=0, preflight_missing=miss)
            continue
        say(f"{ph['name']}: preflight OK — all macro paths resolve")
        ok = []
        for b in range(ph["n_batches"]):
            t_batch = time.time()
            if run_batch(ph, b):
                ok.append(b)
            else:
                # A batch that dies in seconds is a configuration error, not a
                # resource problem; retrying the other 16 cannot help.
                if time.time() - t_batch < 30 and not ok:
                    say(f"{ph['name']}: batch {b} failed in "
                        f"{time.time()-t_batch:.0f}s with no prior success — "
                        f"treating as deterministic, ABANDONING this phantom")
                    break
                say(f"{ph['name']}: batch {b} failed — continuing with remaining batches")
        if not ok:
            say(f"{ph['name']}: NO batches succeeded, skipping concat")
            summary[ph["name"]] = dict(batches_ok=0)
            continue
        det = concat(ph, ok)
        if det:
            filter_rcut(ph, det)
            summary[ph["name"]] = dict(batches_ok=len(ok), batches_total=ph["n_batches"],
                                       n_triples=rows(det))
        with open(OUT / "o4_batched_summary.json", "w") as fh:
            json.dump(summary, fh, indent=1)
    say("=== O4 batched driver DONE ===")
    say(json.dumps(summary))


if __name__ == "__main__":
    main()
