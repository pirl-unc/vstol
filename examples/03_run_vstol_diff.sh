vstol diff \
  --target-tsv-file outputs/colo829-dna-pacbio_savana_classify_output.somatic.tsv \
  --query-tsv-file outputs/hcc1395-dna-pacbio_savana_classify_output.somatic.tsv \
  --output-tsv-file outputs/colo829-dna-pacbio_savana_classify_output.somatic.filtered.tsv \
  --match-both-positions yes \
  --max-breakpoint-distance 1000 \
  --match-operation-types yes
