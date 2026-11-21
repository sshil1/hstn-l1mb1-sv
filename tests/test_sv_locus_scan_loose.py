"""Tests for the loose scan and per-file report (added after a silent '0 hits' on supplements
whose layout is not 'chrom in column 0')."""
import os, tempfile
import openpyxl
import sv_locus_scan as m

WIN = (85459547, 85467444, 5000)


def tmp(name):
    return os.path.join(tempfile.mkdtemp(), name)


def xlsx(rows, name="s.xlsx", sheet="Sheet1"):
    p = tmp(name)
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet
    for r in rows:
        ws.append(r)
    wb.save(p)
    return p


def test_strict_misses_when_chrom_not_in_first_column():
    # documents the old blind spot: this row is a true hit but strict scan cannot see it
    p = xlsx([["SV_id", "size", "chr", "start", "end"], ["x1", 728, "6", 85464318, 85464319]])
    assert m.scan([p], *WIN) == []


def test_loose_finds_hit_with_chrom_in_other_column():
    p = xlsx([["SV_id", "size", "chr", "start", "end"], ["x1", 728, "6", 85464318, 85464319]])
    r = m.scan_loose([p], *WIN)
    assert len(r) == 1 and r[0][1] == "LOOSE_CHR6"


def test_loose_finds_id_style_coordinates():
    p = xlsx([["id"], ["6-85464318-INS-728"]])
    r = m.scan_loose([p], *WIN)
    assert len(r) == 1


def test_loose_integer_without_chrom_is_weaker_flag():
    p = xlsx([["gene", "pos"], ["HSTN", 85464000]])
    r = m.scan_loose([p], *WIN)
    assert len(r) == 1 and r[0][1] == "LOOSE_COORD_ONLY"


def test_loose_ignores_other_chrom_with_far_coordinates():
    p = xlsx([["chr", "start"], ["7", 12345678], ["6", 1000]])
    assert m.scan_loose([p], *WIN) == []


def test_loose_ignores_small_numbers_and_counts():
    p = xlsx([["n_samples", "count"], [85, 728]])
    assert m.scan_loose([p], *WIN) == []


def test_loose_works_on_csv_and_vcf_text():
    c = tmp("a.csv")
    open(c, "w").write("id,chr,pos\nq,6,85464318\n")
    assert len(m.scan_loose([c], *WIN)) == 1


def test_report_counts_rows_and_sheets():
    p = xlsx([["a", "b"], [1, 2], [3, 4]], sheet="Tab1")
    rep = m.report([p])
    assert rep[0]["file"] == p and rep[0]["rows"] == 3 and rep[0]["sheets"] == ["Tab1"]


def test_report_survives_corrupt_xlsx():
    p = tmp("bad.xlsx")
    open(p, "wb").write(b"not a zip")
    rep = m.report([p])
    assert rep[0]["error"]
