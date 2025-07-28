vstol diff \
  --target-tsv-file outputs/colo829-dna-pacbio_savana_classify_output.somatic.chr1.tsv \
  --query-tsv-file outputs/hcc1395-dna-pacbio_savana_classify_output.somatic.chr1.tsv \
  --output-tsv-file outputs/colo829-dna-pacbio_savana_classify_output.somatic.chr1.filtered.tsv \
  --match-both-positions yes \
  --max-breakpoint-distance 1000 \
  --match-operation-types yes
