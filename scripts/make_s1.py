import csv, sys
src = sys.argv[1] if len(sys.argv)>1 else "metadata/cohort_verified_metadata.tsv"
out = sys.argv[2] if len(sys.argv)>2 else "tables/Supplementary_Table_S1.tsv"
rows = list(csv.DictReader(open(src), delimiter="\t"))
def na(v): return "NA" if v in ("", ".") else v
cols = ["sample_label","breed","run_accession","study_accession","secondary_study_accession","biosample_accession",
        "secondary_sample_accession","scientific_name","library_layout","instrument_model","ena_read_count","ena_base_count",
        "locally_recovered_read_count","note"]
with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(cols)
    for r in rows:
        breed = {"DJ":"Dingjie"}.get(r["sample"][:2], r["sample"].split("_")[0])
        w.writerow([r["sample"], breed, r["run_accession"], r["study_accession"], r["secondary_study_accession"], r["sample_accession"],
                    r["secondary_sample_accession"], r["scientific_name"], r["library_layout"], r["instrument_model"],
                    r["read_count"], r["base_count"], na(r["local_verified_read_count"]), na(r["metadata_note"])])
