#!/usr/bin/env bash
# usage: classify_outgroup.sh LABEL FASTA   (relaxed thresholds for divergent outgroups)
set -euo pipefail
export PATH=/opt/conda/bin:/opt/conda/envs/sv_env/bin:/opt/conda/envs/svimasm/bin:$PATH
D=/data/hstn_negctrl/asm_scan; L=$1; A=$2
minimap2 -x asm20 -k15 -w5 -c --secondary=no -t 2 "$A" "$D/flanks.fa" > "$D/$L.og.paf" 2> "$D/$L.og.log"
python3 -I - "$L" "$D/$L.og.paf" <<'PY'
import sys
lab, paf = sys.argv[1:3]
best = {}
for line in open(paf):
    f = line.split("\t"); q = f[0]
    d = dict(strand=f[4], tname=f[5], ts=int(f[7]), te=int(f[8]), nm=int(f[9]), al=int(f[10]), mq=int(f[11]))
    d["id"] = round(100*d["nm"]/d["al"], 1)
    if d["mq"] >= 10 and d["al"] >= 500 and (q not in best or d["nm"] > best[q]["nm"]):
        best[q] = d
miss = [q for q in "LR" if q not in best]
if miss:
    print(f"{lab}\t.\t.\t.\tFLANK_NOT_MAPPED({''.join(miss)})\t.\t.\t.\t.\t.\t."); sys.exit()
l, r = best["L"], best["R"]
if l["tname"] != r["tname"] or l["strand"] != r["strand"]:
    print(f"{lab}\t{l['tname']}|{r['tname']}\t{l['strand']}{r['strand']}\t.\tDISCORDANT\t{l['al']}\t{r['al']}\t{l['mq']}\t{r['mq']}\t{l['id']}\t{r['id']}"); sys.exit()
gap = r["ts"] - l["te"] if l["strand"] == "+" else l["ts"] - r["te"]
cls = "REF" if -20 <= gap <= 30 else ("INS" if 700 <= gap <= 760 else "OTHER")
print(f"{lab}\t{l['tname']}\t{l['strand']}\t{gap}\t{cls}\t{l['al']}\t{r['al']}\t{l['mq']}\t{r['mq']}\t{l['id']}\t{r['id']}")
PY
