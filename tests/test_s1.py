import csv, re, collections
rows = list(csv.DictReader(open("tables/Supplementary_Table_S1.tsv"), delimiter="\t"))
assert len(rows) == 19, len(rows)
assert len({r["sample_label"] for r in rows}) == 19 and len({r["run_accession"] for r in rows}) == 19 and len({r["biosample_accession"] for r in rows}) == 19
assert collections.Counter(r["study_accession"] for r in rows) == {"PRJEB78595":16,"PRJEB63525":2,"PRJNA905718":1}
assert collections.Counter(r["breed"] for r in rows) == {"Gir":5,"Sahiwal":8,"Tharparkar":5,"Dingjie":1}
for r in rows:
    assert re.fullmatch(r"(ERR|SRR)\d{7,8}", r["run_accession"]), r
    assert re.fullmatch(r"SAM(EA|N)\d+", r["biosample_accession"]), r
    assert re.fullmatch(r"(ERS|SRS)\d+", r["secondary_sample_accession"]), r
    assert all("\t" not in v and "\n" not in v for v in r.values())
# cross-check 16 PRJEB78595 runs against the live ENA listing pasted by the user (run -> alias)
ena = """ERR13461127 ERR13492134 ERR13502943 ERR13461117 ERR13461118 ERR13461128 ERR13488875 ERR13488878 ERR13488881 ERR13502944 ERR13453568 ERR13459803 ERR13461125 ERR13488860 ERR13488872 ERR13488877 ERR13488883 ERR13502970 ERR13448261 ERR13459805 ERR13461126 ERR13488880 ERR13488882 ERR13492137 ERR13502968 ERR13448232 ERR13453569 ERR13459804 ERR13461119 ERR13488873 ERR13488876 ERR13488879 ERR13492154""".split()
mine = [r["run_accession"] for r in rows if r["study_accession"] == "PRJEB78595"]
missing = [x for x in mine if x not in ena]
assert not missing, missing
print("S1 OK: 19 rows, 16/2/1 split, all 16 PRJEB78595 runs present in live ENA listing")
