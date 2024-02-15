vstol annotate \
  --tsv-file ../test/data/hg002_sniffles2.tsv \
  --annotator ensembl \
  --output-tsv-file outputs/hg002_sniffles2_ensembl_annotated.tsv \
  --ensembl-release 95 \
  --ensembl-species human

vstol annotate \
  --tsv-file ../test/data/hg002_sniffles2.tsv \
  --annotator gencode \
  --output-tsv-file outputs/hg002_sniffles2_gencode_annotated.tsv \
  --gencode-gtf-file ../test/data/gencode.v41.annotations.gtf \
  --gencode-version v41 \
  --gencode-species human \
  --num-threads 4

vstol annotate \
  --tsv-file ../test/data/hg002_sniffles2.tsv \
  --annotator refseq \
  --output-tsv-file outputs/hg002_sniffles2_refseq_annotated.tsv \
  --refseq-gtf-file ../test/data/refseq_hg38_genomic.gtf \
  --refseq-assembly-report-txt-file ../test/data/refseq_hg38_assembly_report.txt \
  --refseq-version v110 \
  --refseq-species human
