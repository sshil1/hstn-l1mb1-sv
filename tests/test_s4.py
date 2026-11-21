import pathlib
F = pathlib.Path(__file__).resolve().parents[1] / "tables" / "Supplementary_File_S4_candidate_models.fa"

def _read():
    d, n = {}, None
    for l in F.read_text().splitlines():
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = ""
        else:
            d[n] += l.strip()
    return d

def _shared(x, y):
    p = 0
    while p < min(len(x), len(y)) and x[p] == y[p]: p += 1
    s = 0
    while s < min(len(x), len(y)) - p and x[-1 - s] == y[-1 - s]: s += 1
    return p, s

def test_records_and_lengths():
    d = _read()
    assert list(d) == ["REF", "INS_728", "INS_729"]
    assert [len(d[k]) for k in d] == [2000, 2728, 2729]
    assert all(set(v) <= set("ACGT") for v in d.values())

def test_models_replace_two_reference_bases():
    d = _read()
    for k, ins in (("INS_728", 730), ("INS_729", 731)):
        p, s = _shared(d["REF"], d[k])
        assert (p, s) == (1004, 994)
        assert len(d[k]) - p - s == ins
        assert len(d["REF"]) - p - s == 2
