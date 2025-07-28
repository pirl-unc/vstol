vstol annotate \
  --tsv-file outputs/colo829-dna-pacbio_savana_classify_output.somatic.chr1.tsv \
  --annotator gencode \
  --output-tsv-file outputs/colo829-dna-pacbio_savana_classify_output.somatic.gencode_annotated.tsv \
  --gencode-gtf-file ../test/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz \
  --gencode-version v41 \
  --gencode-species human \
  --gencode-levels 1 2 \
  --num-processes 4
