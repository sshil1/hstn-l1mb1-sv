import csv,collections
D="metadata/"
S1={r["sample_label"] for r in csv.DictReader(open("tables/Supplementary_Table_S1.tsv"),delimiter="\t")}
rows=list(csv.DictReader(open("tables/Supplementary_Table_S3.tsv"),delimiter="\t"))
assert len(rows)==19 and {r["sample_label"] for r in rows}==S1
g=[r for r in rows if r["graph_GT"]!="NA"]
assert len(g)==15 and all(r["graph_GT"]=="1/1" for r in g), len(g)
assert sum(r["graph_FT"]=="GQ" for r in g)==5, collections.Counter(r["graph_FT"] for r in g)
assert sum(r["review_status"]!="NA" for r in rows)==5
assert all(r["SNP_85463962_GT"]=="G/G" and r["SNP_85464601_GT"]=="G/G" for r in rows)
for r in rows:   # missing presence calls stay missing, not zero
    if r["graph_GT"]=="NA": assert r["left_bp_INS_reads"]=="NA" and r["graph_GQ"]=="NA"
# breakpoint reads agree with source (30 rows = 15 samples x 2)
src=list(csv.DictReader(open(D+"presence_breakpoints.tsv"),delimiter="\t")); m={r["sample_label"]:r for r in rows}
for b in src:
    k="left" if b["breakpoint"]=="left" else "right"
    assert m[b["sample"]][f"{k}_bp_INS_reads"]==b["INS_count"] and m[b["sample"]][f"{k}_bp_REF_reads"]==b["REF_count"]
for s in ("Tharparkar_17","Gir_885"): assert m[s]["graph_GT"]=="1/1"
print("S3 OK",len(rows),"rows;",len(g),"graph calls;",sum(r["graph_FT"]=="GQ" for r in g),"GQ-flagged")
