from vstolib.vcf import *
from .data import get_data_path


def test_vcf2tsv_manta_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-illumina_manta_outputs/results/variants/somaticSV.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_manta_somatic_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='ilmn'
    )
    assert len(df_variants) == 720


def test_vcf2tsv_delly2_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-illumina-realigned_delly2_filtered.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_delly2_somatic_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='ilmn'
    )
    assert len(df_variants) == 40


def test_vcf2tsv_mutect2_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-illumina-realigned_gatk4-mutect2.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_gatk4_mutect2_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='ilmn'
    )
    assert len(df_variants) == 728042


def test_vcf2tsv_lumpy_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-illumina-realigned_lumpy.chr1.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_lumpy_somatic_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='ilmn'
    )
    assert len(df_variants) == 2304


def test_vcf2tsv_strelka2_indels_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-illumina-realigned_strelka2_indels.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='ilmn'
    )
    assert len(df_variants) == 63584


def test_vcf2tsv_strelka2_snvs_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-illumina-realigned_strelka2_snvs.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='ilmn'
    )
    assert len(df_variants) == 339166


def test_vcf2tsv_clairs_snvs_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_clairs_outputs/snv.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clairs_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 54832


def test_vcf2tsv_clairs_indels_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_clairs_outputs/indel.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clairs_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 180697


def test_vcf2tsv_cutesv_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_cutesv.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_cutesv_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 2304


def test_vcf2tsv_deepvariant_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_deepvariant.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_deepvariant_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 788517


def test_vcf2tsv_pbsv_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_pbsv.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_pbsv_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 3251


def test_vcf2tsv_savana_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_savana_classify_output.somatic.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_savana_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 104


def test_vcf2tsv_severus_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_severus_outputs/somatic_SVs/severus_somatic.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_severus_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 141


def test_vcf2tsv_sniffles2_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_sniffles2.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_sniffles2_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 2525


def test_vcf2tsv_svim_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio_svim.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svim_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 7973


def test_vcf2tsv_svisionpro_1():
    vcf_file = get_data_path(name='vcf/colo829-dna-pacbio.svision_pro_v1.8.s3.somatic_s3.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svisionpro_vcf(
        df_vcf=df_vcf,
        source_id='colo829',
        platform='pacbio'
    )
    assert len(df_variants) == 36


def test_vcf2tsv_manta_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-illumina_manta_outputs/results/variants/somaticSV.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_manta_somatic_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='ilmn'
    )
    assert len(df_variants) == 2246


def test_vcf2tsv_delly2_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-illumina-realigned_delly2_filtered.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_delly2_somatic_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='ilmn'
    )
    assert len(df_variants) == 180


def test_vcf2tsv_mutect2_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-illumina-realigned_gatk4-mutect2.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_gatk4_mutect2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='ilmn'
    )
    assert len(df_variants) == 530992


def test_vcf2tsv_lumpy_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-illumina-realigned_lumpy.chr1.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_lumpy_somatic_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='ilmn'
    )
    assert len(df_variants) == 1678


def test_vcf2tsv_strelka2_indels_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-illumina-realigned_strelka2_indels.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='ilmn'
    )
    assert len(df_variants) == 39276


def test_vcf2tsv_strelka2_snvs_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-illumina-realigned_strelka2_snvs.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_strelka2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='ilmn'
    )
    assert len(df_variants) == 226396


def test_vcf2tsv_clairs_snvs_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_clairs_outputs/snv.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clairs_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 1080782


def test_vcf2tsv_clairs_indels_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_clairs_outputs/indel.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_clairs_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 223467


def test_vcf2tsv_cutesv_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_cutesv.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_cutesv_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 2487


def test_vcf2tsv_deepvariant_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_deepvariant.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_deepvariant_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 723570


def test_vcf2tsv_pbsv_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_pbsv.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_pbsv_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 3547


def test_vcf2tsv_savana_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_savana_classify_output.somatic.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_savana_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 779


def test_vcf2tsv_severus_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_severus_outputs/somatic_SVs/severus_somatic.vcf')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_severus_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 7745


def test_vcf2tsv_sniffles2_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_sniffles2.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_sniffles2_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 2733


def test_vcf2tsv_svim_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio_svim.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svim_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 11512


def test_vcf2tsv_svisionpro_2():
    vcf_file = get_data_path(name='vcf/hcc1395-dna-pacbio.svision_pro_v1.8.s3.somatic_s3.chr1.vcf.gz')
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    df_variants = parse_svisionpro_vcf(
        df_vcf=df_vcf,
        source_id='hcc1395',
        platform='pacbio'
    )
    assert len(df_variants) == 1256

