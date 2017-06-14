from setuptools import find_packages, setup

bio_console_scripts = [
    'bam-to-wiggle=bio_utils.bio_programs.bam_to_wiggle:main',
    'bedx-to-bedy=bio_utils.bio_programs.bedx_to_bedy:main',
    'bed12-to-gtf=bio_utils.bio_programs.bed12_to_gtf:main',
    'calculate-bed-overlap=bio_utils.bio_programs.calculate_bed_overlap:main',
    'convert-ccds-to-bed=bio_utils.bio_programs.convert_ccds_to_bed:main',
    'count-aligned-reads=bio_utils.bio_programs.count_aligned_reads:main',
    'count-reads=bio_utils.bio_programs.count_reads:main',
    'create-aligned-read-count-bar-chart=bio_utils.bio_programs.plotting.create_aligned_read_count_bar_chart:main',
    'create-mygene-report=bio_utils.bio_programs.create_mygene_report:main',
    'dna-to-aa=bio_utils.bio_programs.dna_to_aa:main',
    'download-srr-files=bio_utils.bio_programs.download_srr_files:main',
    'extract-bed-sequences=bio_utils.bio_programs.extract_bed_sequences:main',
    'extract-cds-coordinates=bio_utils.bio_programs.extract_cds_coordinates:main',
    'fasta-to-fastq=bio_utils.bio_programs.fasta_to_fastq:main',
    'fastq-pe-dedupe=bio_utils.bio_programs.fastq_pe_dedupe:main',
    'filter-bam-by-ids=bio_utils.bio_programs.filter_bam_by_ids:main',
    'fix-all-bed-files=bio_utils.bio_programs.fix_all_bed_files:main',
    'get-all-utrs=bio_utils.bio_programs.get_all_utrs:main',
    'get-read-length-distribution=bio_utils.bio_programs.get_read_length_distribution:main',
    'gtf-to-bed12=bio_utils.bio_programs.gtf_to_bed12:main',
    'join-long-chromosomes=bio_utils.bio_programs.join_long_chromosomes:main',
    'merge-isoforms=bio_utils.bio_programs.merge_isoforms:main',
    'parse-meme-names=bio_utils.bio_programs.parse_meme_names:main',
    'plot-read-length-distribution=bio_utils.bio_programs.plotting.plot_read_length_distribution:main',
    'remove-duplicate-bed-entries=bio_utils.bio_programs.remove_duplicate_bed_entries:main',
    'remove-duplicate-sequences=bio_utils.bio_programs.remove_duplicate_sequences:main',
    'remove-multimapping-reads=bio_utils.bio_programs.remove_multimapping_reads:main',
    'reorder-fasta=bio_utils.bio_programs.reorder_fasta:main',
    'run-bowtie=bio_utils.bio_programs.run_bowtie:main',
    'split-bed12-blocks=bio_utils.bio_programs.split_bed12_blocks:main',
    'split-long-chromosomes=bio_utils.bio_programs.split_long_chromosomes:main',
    'subtract-bed=bio_utils.bio_programs.subtract_bed:main',
]

console_scripts = bio_console_scripts

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='bio-utils',
        version='0.2.2',
        description="This repo contains python3 bioinformatics utilities I find useful.",
        long_description=readme(),
        keywords="utilities",
        url="",
        author="Brandon Malone",
        author_email="bmmalone@gmail.com",
        license='MIT',
        packages=find_packages(),
        install_requires = [
            'cython',
            'numpy',
            'scipy',
            'matplotlib',
            'pandas',
            'docopt',
            'tqdm',
            'joblib',
            'xlrd',
            'openpyxl',
            'graphviz',
            'biopython',
            'mygene'
        ],
        include_package_data=True,
        test_suite='nose.collector',
        tests_require=['nose'],
        entry_points = {
            'console_scripts': console_scripts
        },
        zip_safe=False
        )
