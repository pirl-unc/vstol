vstol annotate \
  --tsv-file outputs/hg002_sniffles2.tsv \
  --annotator ensembl \
  --output-tsv-file outputs/hg002_sniffles2_ensembl_annotated.tsv \
  --ensembl-release 95 \
  --ensembl-species human

vstol annotate \
  --tsv-file outputs/hg002_sniffles2.tsv \
  --annotator gencode \
  --output-tsv-file outputs/hg002_sniffles2_gencode_annotated.tsv \
  --gencode-gtf-file ../test/data/gencode.v41.annotations.gtf \
  --gencode-version v41 \
  --gencode-species human \
  --num-threads 4
