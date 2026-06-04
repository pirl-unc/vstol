vstol vcf2tsv \
  --vcf-file ../test/data/vcf/savana/colo829-dna-pacbio_savana_classify_output.somatic.vcf.gz \
  --method savana \
  --platform pacbio \
  --source-id colo829 \
  --output-tsv-file outputs/colo829-dna-pacbio_savana_classify_output.somatic.tsv

vstol vcf2tsv \
  --vcf-file ../test/data/vcf/savana/hcc1395-dna-pacbio_savana_classify_output.somatic.vcf.gz \
  --method savana \
  --platform pacbio \
  --source-id hcc1395 \
  --output-tsv-file outputs/hcc1395-dna-pacbio_savana_classify_output.somatic.tsv
