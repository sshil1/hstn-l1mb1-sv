#!/usr/bin/env bash
# usage: compare_models.sh LABEL [LABEL...]   (needs LABEL.locus.fa)
export PATH=/opt/conda/bin:/opt/conda/envs/sv_env/bin:/opt/conda/envs/svimasm/bin:$PATH
M728=$(find /data /home -name "INS_728.fa" 2>/dev/null | head -1)
M729=$(find /data /home -name "INS_729.fa" 2>/dev/null | head -1)
echo "models: $M728 $M729"
for lab in "$@"; do
  for m in "$M728" "$M729"; do
    minimap2 -x asm20 -c --cs --secondary=no "$m" "$lab.locus.fa" 2>/dev/null | \
      awk -v l="$lab" -v m="$(basename "$m")" '{nm="NA"; for(i=13;i<=NF;i++) if($i ~ /^NM:i:/) nm=$i; print l, m, "alnlen="$11, nm}'
  done
done
