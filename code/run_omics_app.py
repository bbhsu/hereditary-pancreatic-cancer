from os.path import abspath, dirname, exists, join
from subprocess import PIPE, run
from warnings import warn

from genotype_to_phenotype.access_g2p import read_g2p, write_g2p

from tabix import open as tabix_open

PROJECT_DIRECTORY_PATH = dirname(dirname(abspath(__file__)))

DATA_DIRECTORY_PATH = join(PROJECT_DIRECTORY_PATH, 'data')

PERSON_GENOME_VCF_GZ_FILE_PATH = join(DATA_DIRECTORY_PATH, 'person',
                                      'genome.vcf.gz')

OUTPUT_DIRECTORY_PATH = join(PROJECT_DIRECTORY_PATH, 'output')


def run_simple_omics_app():
    """
    Run Simple Omics App.
    Arguments:
    Returns:
    """

    fasta_gz_file_path = join(DATA_DIRECTORY_PATH,
                              'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
    if not exists(fasta_gz_file_path + '.gzi'):
        print('Indexing {} with samtools faidx ...'.format(fasta_gz_file_path))
        run('samtools faidx {}'.format(fasta_gz_file_path),
            shell=True,
            check=True)

    pytabix_handle = tabix_open(PERSON_GENOME_VCF_GZ_FILE_PATH)

    input_g2p = read_g2p(join(DATA_DIRECTORY_PATH, 'input.g2p'))
    matched_rows = []

    for i, row in input_g2p['table'].iterrows():
        print('Matching input.g2p row {} ...'.format(i))

        types, features, regions, states = row[:4]

        matches = []
        for type_, feature, region, state in zip(
                types.split(';'),
                features.split(';'), regions.split(';'), states.split(';')):

            chromosome, start_position_end_position = region.split(':')

            try:
                start_position, end_position = (
                    int(position)
                    for position in start_position_end_position.split('-'))

                if end_position < start_position:
                    raise ValueError(
                        'REGION end position ({}) cannot be less than the start position ({}) (input.g2p row {}).'.
                        format(end_position, start_position, i))

            except ValueError:
                raise ValueError(
                    'REGION should be writted as <chromosome>:<start_position>-<end_position> or <chromosome>:<start_positioln>-<start_position> (input.g2p row {}).'.
                    format(i))

            variants = tuple(
                pytabix_handle.query(chromosome, start_position - 1,
                                     end_position))

            if type_ == 'variant':
                if len(variants):

                    variant = variants[0]
                    _validate_variant(variant)

                    ref, alt = variant[3], variant[4]
                    alleles = (
                        ref,
                        *alt.split(','), )

                    format_, sample = variant[8], variant[9]
                    format_fields = format_.split(':')

                    if 'GT' not in format_fields:
                        raise ValueError(
                            'Cannot determining genotype because FORMAT does not have GT field ({}).'.
                            format(variant))

                    gt = sample.split(':')[format_fields.index('GT')]

                    gt_indices = (int(i)
                                  for i in gt.replace('/', '|').split('|'))

                    sample_genotype = (alleles[i] for i in gt_indices)

                else:
                    sequence = ''.join(
                        run('samtools faidx {} {}:{}-{}'.format(
                            fasta_gz_file_path, chromosome, start_position,
                            end_position),
                            shell=True,
                            check=True,
                            stdout=PIPE).stdout.decode().strip().split('\n')[
                                1:])

                    sample_genotype = (sequence, ) * 2

                try:
                    alleles = state.split('|')
                except ValueError:
                    raise ValueError(
                        'Variant STATE should be writted as <allele_1>|<allele_2> (input.g2p row {}).'.
                        format(i))

                matches.append(set(sample_genotype) == set(alleles))

            elif type_ == 'gene':

                impact_to_int = {
                    'HIGH': 3,
                    'MODERATE': 2,
                    'LOW': 1,
                    'MODIFIER': 0,
                }

                if len(variants):

                    impact_ints = []

                    for variant in variants:
                        _validate_variant(variant)

                        try:
                            anns = [
                                info[4:] for info in variant[7].split(';')
                                if info.startswith('ANN=')
                            ][0]

                            ann_0 = anns.split(',')[0]
                            impact = ann_0.split('|')[2]

                            impact_ints.append(impact_to_int[impact])

                        except IndexError:
                            warn(
                                'Cannot get gene information because this variant is not annotated and INFO is missing ANN.'
                            )

                    if len(impact_ints):
                        matches.append(
                            impact_to_int.get(state, 2) <= max(impact_ints))
                    else:
                        matches.append(False)

                else:
                    matches.append(False)

            else:
                raise ValueError(
                    'Unknown TYPE {} (input.g2p row {}).'.format(type_, i))

        if all(matches):
            matched_rows.append(i)

    output_g2p_table = input_g2p['table'].loc[matched_rows]

    write_g2p({
        'header': input_g2p['header'],
        'table': output_g2p_table,
    }, join(OUTPUT_DIRECTORY_PATH, 'output.g2p'))

    print('This Omics App output {}.'.format(
        join(OUTPUT_DIRECTORY_PATH, 'output.g2p')))


def _validate_variant(variant):
    """Validate variant."""

    if len(variant) < 10:
        raise ValueError(
            '.vcf file columns do not contain all of CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and FORMET, and at least 1 sample.'.
            format(variant))

    elif 10 < len(variant):
        raise NotImplementedError(
            'Multi-sample .vcf file is not supported yet.')


run_simple_omics_app()
