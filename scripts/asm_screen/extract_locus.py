import gzip, sys
for lab in sys.argv[1:]:
    best = {}
    for line in open(f"{lab}.paf"):
        f = line.split("\t")
        q, tname, ts, te, nm, al, mq = f[0], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10]), int(f[11])
        if mq >= 20 and al >= 800 and (q not in best or nm > best[q][3]):
            best[q] = (tname, ts, te, nm)
    ctg, s, e = best["L"][0], best["L"][1], best["R"][2]
    assert best["L"][0] == best["R"][0], "flanks on different contigs"
    seq, on = [], False
    for l in gzip.open(f"/data/cattle_zebu/asm/{lab}.fa.gz", "rt"):
        if l[0] == ">":
            if on: break
            on = (l[1:].split()[0] == ctg); continue
        if on: seq.append(l.strip())
    seq = "".join(seq)[s:e].upper()
    open(f"{lab}.locus.fa", "w").write(f">{lab}\n{seq}\n")
    print(lab, ctg, s, e, len(seq))
