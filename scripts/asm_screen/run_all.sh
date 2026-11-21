#!/usr/bin/env bash
cd /data/hstn_negctrl/asm_scan
R=results.tsv
[ -s $R ] || printf "label\tcontig\tstrand\tgap_bp\tclass\tL_aln\tR_aln\tL_mapq\tR_mapq\n" > $R
( grep -E '^TAU_' /data/cattle_zebu/asm_list.tsv; grep -E '^ZEB_' /data/cattle_zebu/asm_list.tsv ) | \
while IFS=$'\t' read -r lab acc || [ -n "$lab" ]; do
  grep -q "^$lab"$'\t' $R && { echo "skip $lab"; continue; }
  ./classify_asm.sh "$lab" "/data/cattle_zebu/asm/$lab.fa.gz" 4 >> $R || echo "FAIL $lab" >> failed.txt
done
echo ASMSCAN_DONE
