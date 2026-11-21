import csv, os, sys
D=(sys.argv[1] if len(sys.argv)>1 else "metadata").rstrip("/")+"/"
rd=lambda f:list(csv.DictReader(open(D+f),delimiter="\t"))
integ=rd("integrated_cohort.tsv"); pres={r["sample"]:r for r in rd("presence_cohort.tsv")}
flag={r["sample"]:r for r in rd("flagged_review_summary.tsv")}
bp={}
for r in rd("presence_breakpoints.tsv"): bp.setdefault(r["sample"],{})[r["breakpoint"]]=r
na=lambda v:"NA" if v in (None,"","." ) else v
cols=["sample_label","exact_REF_junction","exact_ALT_left_728","exact_ALT_left_729like","exact_ALT_right_shared",
"SNP_85463962_GT","SNP_85463962_DP","SNP_85463962_AD","SNP_85464601_GT","SNP_85464601_DP","SNP_85464601_AD",
"graph_GT","graph_allele_classes","graph_FT","graph_GQ","graph_DP",
"left_bp_REF_reads","left_bp_INS_reads","left_bp_GQ","left_bp_filter","right_bp_REF_reads","right_bp_INS_reads","right_bp_GQ","right_bp_filter",
"strict_REF_fragments","strict_INS_left_fragments","strict_INS_right_fragments","strict_conflicts","review_status","review_note","presence_note"]
with open(sys.argv[2] if len(sys.argv)>2 else "tables/Supplementary_Table_S3.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t",lineterminator="\n"); w.writerow(cols)
    for r in integ:
        s=r["sample"]; p=pres.get(s,{}); fl=flag.get(s,{}); L=bp.get(s,{}).get("left",{}); R=bp.get(s,{}).get("right",{})
        w.writerow([s,r["REF_junction"],r["ALT_left_728"],r["ALT_left_729like"],r["ALT_right_shared"],
        r["genotype_85463962"],r["DP_85463962"],r["AD_85463962"],r["genotype_85464601"],r["DP_85464601"],r["AD_85464601"],
        na(p.get("presence_GT")),na(p.get("caller_allele_classes")),na(p.get("FT")),na(p.get("JSON_GQ")),na(p.get("DP")),
        na(L.get("REF_count")),na(L.get("INS_count")),na(L.get("JSON_GQ")),na(L.get("filters")),
        na(R.get("REF_count")),na(R.get("INS_count")),na(R.get("JSON_GQ")),na(R.get("filters")),
        na(fl.get("REF_fragments")),na(fl.get("INS_left_fragments")),na(fl.get("INS_right_fragments")),na(fl.get("conflicts")),na(fl.get("status")),
        na(r["review_note"]),na(p.get("note"))])
