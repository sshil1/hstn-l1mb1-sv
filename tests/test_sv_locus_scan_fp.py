"""Regression tests for false positives seen on Leonard 2022 MOESM4 (CNV table with LOC gene ids)."""
import os, tempfile
import sv_locus_scan as m

WIN = (85459547, 85467444, 5000)


def csvf(txt):
    p = os.path.join(tempfile.mkdtemp(), "t.csv")
    open(p, "w").write(txt)
    return p


def test_gene_id_digits_are_not_coordinates():
    # LOC100336682 -> 100336682 is inside the numeric window but is an identifier, not a position
    p = csvf("gene,chr,start,end,n\nLOC100336682,3,8564173,8568747,6\n")
    assert m.scan_loose([p], *WIN) == []


def test_count_column_six_does_not_make_row_chr6():
    p = csvf("gene,chr,start,end,n\nLOC100336682,3,8564173,8568747,6\nLOC999,9,84722752,84733970,6\n")
    assert m.scan_loose([p], *WIN) == []


def test_true_spanning_sv_still_found():
    p = csvf("id,chr,start,end\nDEL1,6,85000000,85900000\n")
    r = m.scan_loose([p], *WIN)
    assert len(r) == 1 and r[0][1] == "LOOSE_CHR6"


def test_true_in_window_row_still_found():
    p = csvf("gene,chr,start,end\nHSTN,6,85464000,85464500\n")
    assert len(m.scan_loose([p], *WIN)) == 1
