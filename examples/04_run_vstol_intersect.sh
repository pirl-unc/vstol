vstol intersect \
  --tsv-file outputs/hg002_cutesv.tsv.gz \
  --tsv-file outputs/hg002_pbsv.tsv \
  --tsv-file outputs/hg002_sniffles2.tsv \
  --tsv-file outputs/hg002_svim.tsv \
  --num-threads 4 \
  --max-neighbor-distance 10 \
  --match-all-breakpoints yes \
  --match-variant-types yes \
  --output-tsv-file outputs/hg002_intersecting_variants.tsv

