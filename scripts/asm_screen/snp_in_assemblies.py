import gzip, re, sys
D = "/data/hstn_negctrl/asm_scan"
win = "".join(l.strip() for l in open(f"{D}/ARSUCD2_window.fa") if l[0] != ">").upper()
W0 = 85462323
SNPS = [("SNP1_85463962", 85463962), ("SNP2_85464601", 85464601)]
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

def qinfo(p):
    i = p - W0
    if 1000 <= i < 2000: return "L", i - 1000
    if 2002 <= i < 3002: return "R", i - 2002
    raise ValueError(p)

def best_alignments(paf):
    best = {}
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
        d = dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], tname=f[5], ts=int(f[7]),
                 nm=int(f[9]), al=int(f[10]), mq=int(f[11]), cg=cg)
        if cg and d["mq"] >= 20 and d["al"] >= 800 and (f[0] not in best or d["nm"] > best[f[0]]["nm"]):
            best[f[0]] = d
    return best

def target_pos(a, qoff):
    if a["strand"] == "+": q, qcur = qoff, a["qs"]
    else: q, qcur = a["qlen"] - 1 - qoff, a["qlen"] - a["qe"]
    tcur = a["ts"]
    for n, op in re.findall(r"(\d+)([MID])", a["cg"]):
        n = int(n)
        if op == "M":
            if qcur <= q < qcur + n: return tcur + (q - qcur)
            qcur += n; tcur += n
        elif op == "I":
            if qcur <= q < qcur + n: return None
            qcur += n
        else: tcur += n
    return None

def load_contig(path, name):
    op = gzip.open if path.endswith(".gz") else open
    seq, on = [], False
    with op(path, "rt") as fh:
        for l in fh:
            if l[0] == ">":
                if on: break
                on = (l[1:].split()[0] == name); continue
            if on: seq.append(l.strip())
    return "".join(seq).upper()

def snp_bases(label, path):
    best = best_alignments(f"{D}/{label}.paf")
    if "L" not in best or "R" not in best or best["L"]["tname"] != best["R"]["tname"]: return ["NA", "NA"]
    seq = load_contig(path, best["L"]["tname"])
    out = []
    for _, pos in SNPS:
        arm, off = qinfo(pos)
        a = best[arm]; tp = target_pos(a, off)
        if tp is None: out.append("NA"); continue
        b = seq[tp]
        out.append(b if a["strand"] == "+" else COMP.get(b, "N"))
    return out

cls = {}
for l in open(f"{D}/results.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] != "label": cls[f[0]] = f[4]
labels = sys.argv[1:] or list(cls)
print("label\tclass\t" + "\t".join(n for n, _ in SNPS), flush=True)
for lab in labels:
    if lab == "SELFTEST":
        exp = [win[p - W0] for _, p in SNPS]
        got = snp_bases("SELFTEST", f"{D}/ARSUCD2_window.fa")
        print(f"SELFTEST\tREF\t" + "\t".join(got) + f"\t(expected {exp}; {'PASS' if got == exp else 'FAIL'})", flush=True)
        continue
    if cls.get(lab) == "NO_MAP":
        print(f"{lab}\tNO_MAP\tNA\tNA", flush=True); continue
    print(f"{lab}\t{cls[lab]}\t" + "\t".join(snp_bases(lab, f"/data/cattle_zebu/asm/{lab}.fa.gz")), flush=True)
