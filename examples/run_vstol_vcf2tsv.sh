# cutesv
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_cutesv.vcf \
  --variant-calling-method cutesv \
  --sequencing-platform PACBIO \
  --source-id hg002 \
  --output-tsv-file outputs/hg002_cutesv.tsv

# deepvariant
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_deepvariant.vcf \
  --variant-calling-method deepvariant \
  --sequencing-platform PACBIO \
  --source-id hg002 \
  --output-tsv-file outputs/hg002_deepvariant.tsv

# delly2
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_delly2.vcf \
  --variant-calling-method delly2-somatic \
  --sequencing-platform ILMN \
  --source-id hg002 \
  --output-tsv-file outputs/hg002_hg001_delly2.tsv

# gatk4-mutect2
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_gatk4_mutect2.vcf \
  --variant-calling-method gatk4-mutect2 \
  --sequencing-platform ILMN \
  --source-id hg002 \
  --output-tsv-file outputs/hg002_hg001_gatk4_mutect2.tsv

# lumpy
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_lumpy.vcf \
  --variant-calling-method lumpy-somatic \
  --sequencing-platform ILMN \
  --source-id hg002 \
  --output-tsv-file outputs/hg002_hg001_lumpy.tsv

# pbsv
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_pbsv.vcf \
  --variant-calling-method pbsv \
  --sequencing-platform PACBIO \
  --source-id hg002 \
  --output-tsv-file outputs/hg002_pbsv.tsv

# savana
vstol vcf2tsv \
  --vcf-file ../test/data/sample001_savana_sv_breakpoints.vcf \
  --variant-calling-method savana \
  --sequencing-platform PACBIO \
  --source-id sample001 \
  --case-id tumor \
  --control-id normal \
  --output-tsv-file outputs/sample001_savana_sv_breakpoints.tsv

# severus
vstol vcf2tsv \
  --vcf-file ../test/data/sample001_severus_somatic.vcf \
  --variant-calling-method severus \
  --sequencing-platform PACBIO \
  --source-id sample001 \
  --case-id tumor \
  --output-tsv-file outputs/sample001_severus_somatic.tsv

# sniffles2
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_sniffles2.vcf \
  --variant-calling-method sniffles2 \
  --sequencing-platform PACBIO \
  --source-id hg002 \
  --output-tsv-file outputs/hg002_sniffles2.tsv

# strelka2-indels
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_strelka2_indels.vcf \
  --variant-calling-method strelka2-somatic \
  --sequencing-platform ILMN \
  --source-id hg002 \
  --case-id hg002 \
  --control-id hg001 \
  --output-tsv-file outputs/hg002_hg001_strelka2_indels.tsv

# strelka2-snvs
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_hg001_strelka2_snvs.vcf \
  --variant-calling-method strelka2-somatic \
  --sequencing-platform ILMN \
  --source-id hg002 \
  --case-id hg002 \
  --control-id hg001 \
  --output-tsv-file outputs/hg002_hg001_strelka2_snvs.tsv

# svim
vstol vcf2tsv \
  --vcf-file ../test/data/hg002_svim.vcf \
  --variant-calling-method svim \
  --sequencing-platform PACBIO \
  --source-id hg002 \
  --output-tsv-file outputs/hg002_svim.tsv

# svisionpro
vstol vcf2tsv \
  --vcf-file ../test/data/sample001.svision_pro_v1.8.s3.vcf \
  --variant-calling-method svisionpro \
  --sequencing-platform PACBIO \
  --source-id sample001 \
  --case-id tumor \
  --control-id normal \
  --output-tsv-file outputs/sample001.svision_pro_v1.8.s3.tsv