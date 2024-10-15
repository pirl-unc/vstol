vstol subtract \
  --target-tsv-file outputs/hg002_sniffles2.tsv \
  --query-tsv-file outputs/hg002_pbsv.tsv \
  --query-tsv-file outputs/hg002_svim.tsv \
  --query-tsv-file outputs/hg002_cutesv.tsv.gz \
  --num-threads 1 \
  --max-neighbor-distance 10 \
  --match-all-breakpoints no \
  --match-variant-types no \
  --output-tsv-file outputs/hg002_sniffles2_private.tsv

vstol subtract \
  --target-tsv-file outputs/hg002_merged_variants.tsv \
  --query-tsv-file outputs/hg002_pbsv.tsv \
  --query-tsv-file outputs/hg002_svim.tsv \
  --query-tsv-file outputs/hg002_cutesv.tsv.gz \
  --num-threads 1 \
  --max-neighbor-distance 10 \
  --match-all-breakpoints no \
  --match-variant-types no \
  --output-tsv-file outputs/hg002_merged_variants_private.tsv
