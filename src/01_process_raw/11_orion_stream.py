"""Phase O2: stream X-Atlas/Orion one batch file at a time, keep control cells and per-target sums, delete the file.

usage: python3 o2_orion_stream.py HCT116|HEK293T
Resumable: a killed run (shutdown, sleep, lost network) restarts from the last checkpoint.
  data/orion/state/<line>/genes.npy        token ids of the gene set (fixed from the first file's Non-Targeting cells)
  data/orion/state/<line>/ctrl/<batch>.npz Non-Targeting cells: CSR on the gene set (u16 counts), library size; written atomically
  data/orion/state/<line>/acc_*.npy        per (target, half) sums of CP10k and log1p(CP10k), cell counts; saved every 10 files
  data/orion/state/<line>/done.json        files whose sums are inside the saved accumulators (saved after them)
Download: curl over IPv4 with resume and retries, two files fetched ahead of processing.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import DATA, SOURCE, SCORER  # noqa: E402
import os, sys, json, time, zlib, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pyarrow.parquet as pq, requests
import urllib3.util.connection as _c
_c.HAS_IPV6 = False

LINE = sys.argv[1]
REPO = "https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion"
ST = os.path.join(os.environ.get("ORION_STATE", str(DATA / "orion" / "state")), LINE); TMP = os.path.join(ST, "tmp"); CT = os.path.join(ST, "ctrl")
for d in (ST, TMP, CT):
    os.makedirs(d, exist_ok=True)
MAXP, NTOK, EVERY = 20000, 38606, int(os.environ.get("ORION_EVERY", 10))
COLS = ["gene_token_id", "gene_expression", "total_counts", "gene_target", "cell_barcode", "pass_guide_filter"]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def file_list():
    for k in range(10):
        try:
            x = requests.get(REPO.replace("/datasets/", "/api/datasets/") + "/tree/main/data", timeout=60).json()
            fs = [f["path"] for f in x if f["type"] == "file" and os.path.basename(f["path"]).startswith(LINE + "_")]
            return sorted(fs, key=lambda p: int(p.split("Batch")[1].split(".")[0]))
        except Exception as e:
            log("file list retry", e); time.sleep(30)
    raise SystemExit("cannot list files")


def fetch(path):
    out = os.path.join(TMP, os.path.basename(path))
    if os.path.exists(out):
        return out
    part = out + ".part"
    for k in range(200):
        r = subprocess.run(["curl", "-4", "-sL", "-C", "-", "--retry", "5", "--speed-limit", "1000000", "--speed-time", "30", "-o", part,
                            f"{REPO}/resolve/main/{path}"])
        if r.returncode == 0:
            try:
                pq.ParquetFile(part).metadata            # complete parquet footer = complete file
                os.replace(part, out)
                return out
            except Exception:
                pass
        time.sleep(5 if k < 10 else 60)               # a stalled slow edge is dropped, the resume reconnects
    raise IOError(f"download failed: {path}")


def read(fn):
    t = pq.ParquetFile(fn).read(columns=COLS)
    keep = np.asarray(t["pass_guide_filter"].to_numpy() == 1)
    ex, tk = t["gene_expression"].combine_chunks(), t["gene_token_id"].combine_chunks()
    off = ex.offsets.to_numpy()
    return dict(v=ex.values.to_numpy(), tok=tk.values.to_numpy(), off=off, keep=keep,
                total=t["total_counts"].to_numpy(), target=t["gene_target"].to_numpy(zero_copy_only=False),
                bc=t["cell_barcode"].to_numpy(zero_copy_only=False))


files = file_list()
if os.environ.get("ORION_FILES"):                               # test runs: a fixed list of batches
    files = [f for f in files if os.path.basename(f) in os.environ["ORION_FILES"].split(",")]
log(f"{LINE}: {len(files)} files")
gpath = os.path.join(ST, "genes.npy")
if not os.path.exists(gpath):                                  # gene set from the first file's controls
    d = read(fetch(files[0]))
    ctrl = np.where(d["keep"] & (d["target"] == "Non-Targeting"))[0]
    s = np.zeros(NTOK)
    for c in ctrl:
        a, b = d["off"][c], d["off"][c + 1]
        np.add.at(s, d["tok"][a:b], d["v"][a:b] * 1e4 / d["total"][c])
    genes = np.where(s / len(ctrl) >= 0.5)[0]
    np.save(gpath, genes)
    log(f"gene set: {len(genes)} genes from {len(ctrl)} control cells")
genes = np.load(gpath); G = len(genes)
lut = np.full(NTOK, -1, np.int64); lut[genes] = np.arange(G)

# accumulators and checkpoint
ap = {k: os.path.join(ST, f"acc_{k}.npy") for k in ("cp", "lg", "n")}
dp, tp = os.path.join(ST, "done.json"), os.path.join(ST, "targets.json")
if os.path.exists(dp):
    done, targets = set(json.load(open(dp))), json.load(open(tp))
    acc = {k: np.load(ap[k]) for k in ap}
else:
    done, targets = set(), {}
    acc = dict(cp=np.zeros((2 * MAXP, G), np.float32), lg=np.zeros((2 * MAXP, G), np.float32), n=np.zeros(2 * MAXP, np.int64))


def checkpoint():
    for k in acc:
        np.save(ap[k] + ".tmp.npy", acc[k]); os.replace(ap[k] + ".tmp.npy", ap[k])
    json.dump(targets, open(tp + ".tmp", "w")); os.replace(tp + ".tmp", tp)
    json.dump(sorted(done), open(dp + ".tmp", "w")); os.replace(dp + ".tmp", dp)      # last: marks the sums as saved


def process(fn, name):
    d = read(fn)
    n = len(d["total"])
    cell = np.repeat(np.arange(n), np.diff(d["off"]))
    col = lut[d["tok"]]
    ok = (col >= 0) & d["keep"][cell]
    cell, col, cnt = cell[ok], col[ok], d["v"][ok]
    cp = (cnt * 1e4 / d["total"][cell]).astype(np.float32)
    lg = np.log1p(cp)
    is_ctrl = d["target"] == "Non-Targeting"
    # controls: CSR on the gene set, written atomically (doubles as a per-file checkpoint)
    cc = np.where(d["keep"] & is_ctrl)[0]
    m = is_ctrl[cell]
    ptr = np.zeros(len(cc) + 1, np.int64); ptr[1:] = np.cumsum(np.bincount(cell[m], minlength=n)[cc])
    out = os.path.join(CT, name + ".npz")
    np.savez(out + ".tmp.npz", indptr=ptr, idx=col[m].astype(np.uint16), cnt=np.minimum(cnt[m], 65535).astype(np.uint16),
             lib=d["total"][cc].astype(np.float32))
    os.replace(out + ".tmp.npz", out)
    # perturbed cells: sums per (target, half)
    pc = np.where(d["keep"] & ~is_ctrl)[0]
    half = np.array([zlib.crc32(b.encode()) & 1 for b in d["bc"][pc]])
    tid = np.array([targets.setdefault(t, len(targets)) for t in d["target"][pc]])
    if len(targets) > MAXP:
        raise SystemExit("too many targets for the accumulator")
    key_of_cell = np.full(n, -1); key_of_cell[pc] = 2 * tid + half
    k = key_of_cell[cell]; m = k >= 0
    keys, kc = np.unique(k[m], return_inverse=True)
    flat = kc.astype(np.int64) * G + col[m]
    for name_, w in (("cp", cp[m]), ("lg", lg[m])):
        s = np.bincount(flat, weights=w, minlength=len(keys) * G).reshape(len(keys), G)
        acc[name_][keys] += s.astype(np.float32)
    np.add.at(acc["n"], 2 * tid + half, 1)
    return len(cc), len(pc)


todo = [f for f in files if os.path.basename(f) not in done]
log(f"{len(done)} done, {len(todo)} to go, gene set {G}")
t0, since = time.time(), 0
with ThreadPoolExecutor(4) as ex:
    futs = {f: ex.submit(fetch, f) for f in todo[:4]}
    for i, f in enumerate(todo):
        if i + 4 < len(todo):
            futs[todo[i + 4]] = ex.submit(fetch, todo[i + 4])
        fn = futs.pop(f).result()
        nc, npc = process(fn, os.path.basename(f))
        os.remove(fn)
        done.add(os.path.basename(f)); since += 1
        log(f"{os.path.basename(f)}: {nc} control, {npc} perturbed cells | {len(done)}/{len(files)} | {(time.time()-t0)/60:.0f} min")
        if since >= EVERY or i == len(todo) - 1:
            checkpoint(); since = 0
log("finished", LINE)
