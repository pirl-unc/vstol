vstol vcf2tsv \
  --vcf-file ../test/data/vcf/colo829-dna-pacbio_savana_classify_output.somatic.chr1.vcf.gz \
  --method savana \
  --platform pacbio \
  --source-id colo829 \
  --output-tsv-file outputs/colo829-dna-pacbio_savana_classify_output.somatic.chr1.tsv

vstol vcf2tsv \
  --vcf-file ../test/data/vcf/hcc1395-dna-pacbio_savana_classify_output.somatic.chr1.vcf.gz \
  --method savana \
  --platform pacbio \
  --source-id hcc1395 \
  --output-tsv-file outputs/hcc1395-dna-pacbio_savana_classify_output.somatic.chr1.tsv
