# hstn-l1mb1-sv

Code and supplementary tables for the manuscript *Sequence-level characterization of an L1MB1-associated intronic structural variant in bovine HSTN across indicine and taurine cattle genomes*.

## Contents
- `tables/` Supplementary Tables S1 (19 animals, accessions), S2 (39 assemblies, flank-anchored screen) and S3 (per-animal evidence).
- `metadata/` Input tables used to build S1 and S3 (accessions and derived counts only; no raw reads).
- `scripts/make_s1.py`, `scripts/make_s3.py` Rebuild S1 and S3 from `metadata/`.
- `scripts/sv_locus_scan.py` Read-only scanner that looks for chr6 coordinates inside a locus window in supplementary SV tables (no network, exec or write calls).
- `tests/` pytest suite (TDD); run `pytest -q` from the repository root.

## Reproduce
```bash
python3 scripts/make_s1.py metadata/cohort_verified_metadata.tsv tables/Supplementary_Table_S1.tsv
python3 scripts/make_s3.py metadata tables/Supplementary_Table_S3.tsv
pytest -q
```

## Scope and caveats
Coordinates are NC_037333.1 (ARS-UCD2.0). The 19 animals are not a population survey. Graph genotypes are presence-only computational predictions, with no dosage. S2 was produced on the authors' server; the alignment pipeline scripts are not yet in this repository.

## Data
Public: PRJEB78595, PRJEB63525, PRJNA905718. Cattle pangenome record: Yang et al., 2026 (Zenodo 15313439).

## License
MIT. See `LICENSE`.

## Assembly screen (S2)
`scripts/asm_screen/` holds the flank-anchored minimap2 screen used for S2 (`run_all.sh` calls `classify_asm.sh`; `classify_outgroup.sh`, `compare_models.sh`, `extract_locus.py` and `snp_in_assemblies.py` produce the outgroup, model-comparison and SNP columns). The scripts use the authors' server paths (`/data/hstn_negctrl/asm_scan`, `/data/cattle_zebu/asm`); edit these to your layout. Download the 39 assemblies from NCBI using the accessions in `inputs/asm_list.tsv`. Inputs: `inputs/ARSUCD2_window.fa` (NC_037333.1:85462323-85466324) and `inputs/flanks.fa`. Requires minimap2.
