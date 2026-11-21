#!/usr/bin/env bash
# usage: classify_asm.sh LABEL ASSEMBLY(.fa|.fa.gz) [THREADS]
set -euo pipefail
export PATH=/opt/conda/bin:/opt/conda/envs/sv_env/bin:/opt/conda/envs/svimasm/bin:$PATH
command -v minimap2 >/dev/null || { echo "missing minimap2" >&2; exit 1; }
D=/data/hstn_negctrl/asm_scan; L=$1; A=$2; T=${3:-4}
nice -n 19 minimap2 -x asm20 -c --secondary=no -t "$T" "$A" "$D/flanks.fa" > "$D/$L.paf" 2> "$D/$L.mm2.log"
python3 -I - "$L" "$D/$L.paf" <<'PY'
import sys
lab, paf = sys.argv[1:3]
best = {}
for line in open(paf):
    f = line.split("\t"); q = f[0]
    d = dict(strand=f[4], tname=f[5], ts=int(f[7]), te=int(f[8]), nm=int(f[9]), al=int(f[10]), mq=int(f[11]))
    if d["mq"] >= 20 and d["al"] >= 800 and (q not in best or d["nm"] > best[q]["nm"]):
        best[q] = d
if "L" not in best or "R" not in best:
    print(f"{lab}\t.\t.\t.\tNO_MAP\t.\t.\t.\t."); sys.exit()
l, r = best["L"], best["R"]
if l["tname"] != r["tname"] or l["strand"] != r["strand"]:
    print(f"{lab}\t{l['tname']}|{r['tname']}\t{l['strand']}{r['strand']}\t.\tDISCORDANT\t{l['al']}\t{r['al']}\t{l['mq']}\t{r['mq']}"); sys.exit()
gap = r["ts"] - l["te"] if l["strand"] == "+" else l["ts"] - r["te"]
cls = "REF" if -20 <= gap <= 30 else ("INS" if 700 <= gap <= 760 else "OTHER")
print(f"{lab}\t{l['tname']}\t{l['strand']}\t{gap}\t{cls}\t{l['al']}\t{r['al']}\t{l['mq']}\t{r['mq']}")
PY
