vstol diff \
  --target-tsv-file ../test/data/hg002_sniffles2.tsv \
  --query-tsv-file ../test/data/hg002_pbsv.tsv \
  --query-tsv-file ../test/data/hg002_svim.tsv \
  --query-tsv-file ../test/data/hg002_cutesv.tsv \
  --num-threads 1 \
  --max-neighbor-distance 10 \
  --output-tsv-file outputs/hg002_sniffles2_private.tsv

vstol diff \
  --target-tsv-file ../test/data/hg002_merged_variants.tsv \
  --query-tsv-file ../test/data/hg002_pbsv.tsv \
  --query-tsv-file ../test/data/hg002_svim.tsv \
  --query-tsv-file ../test/data/hg002_cutesv.tsv \
  --num-threads 1 \
  --max-neighbor-distance 10 \
  --output-tsv-file outputs/hg002_merged_variants_private.tsv
