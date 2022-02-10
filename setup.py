from setuptools import find_packages, setup

def _safe_read_lines(f, remove_git_lines=True):
    with open(f) as in_f:
        r = in_f.readlines()
    r = [l.strip() for l in r]

    if remove_git_lines:
        r = [ l for l in r if not l.startswith("git+ssh") ]

    return r

bio_console_scripts = [
    'bam-to-wiggle=lifesci.bio_programs.bam_to_wiggle:main',
    'bedx-to-bedy=lifesci.bio_programs.bedx_to_bedy:main',
    'bed12-to-gtf=lifesci.bio_programs.bed12_to_gtf:main',
    'calculate-bed-overlap=lifesci.bio_programs.calculate_bed_overlap:main',
    'convert-ccds-to-bed=lifesci.bio_programs.convert_ccds_to_bed:main',
    'count-aligned-reads=lifesci.bio_programs.count_aligned_reads:main',
    'count-reads=lifesci.bio_programs.count_reads:main',
    'create-aligned-read-count-bar-chart=lifesci.bio_programs.plotting.create_aligned_read_count_bar_chart:main',
    'create-mygene-report=lifesci.bio_programs.create_mygene_report:main',
    'dna-to-aa=lifesci.bio_programs.dna_to_aa:main',
    'download-srr-files=lifesci.bio_programs.download_srr_files:main',
    'extract-bed-sequences=lifesci.bio_programs.extract_bed_sequences:main',
    'extract-cds-coordinates=lifesci.bio_programs.extract_cds_coordinates:main',
    'fasta-to-fastq=lifesci.bio_programs.fasta_to_fastq:main',
    'fastq-pe-dedupe=lifesci.bio_programs.fastq_pe_dedupe:main',
    'filter-bam-by-ids=lifesci.bio_programs.filter_bam_by_ids:main',
    'fix-all-bed-files=lifesci.bio_programs.fix_all_bed_files:main',
    'get-all-utrs=lifesci.bio_programs.get_all_utrs:main',
    'get-read-length-distribution=lifesci.bio_programs.get_read_length_distribution:main',
    'gtf-to-bed12=lifesci.bio_programs.gtf_to_bed12:main',
    'join-long-chromosomes=lifesci.bio_programs.join_long_chromosomes:main',
    'merge-isoforms=lifesci.bio_programs.merge_isoforms:main',
    'parse-meme-names=lifesci.bio_programs.parse_meme_names:main',
    'plot-read-length-distribution=lifesci.bio_programs.plotting.plot_read_length_distribution:main',
    'remove-duplicate-bed-entries=lifesci.bio_programs.remove_duplicate_bed_entries:main',
    'remove-duplicate-sequences=lifesci.bio_programs.remove_duplicate_sequences:main',
    'remove-multimapping-reads=lifesci.bio_programs.remove_multimapping_reads:main',
    'reorder-fasta=lifesci.bio_programs.reorder_fasta:main',
    'run-bowtie=lifesci.bio_programs.run_bowtie:main',
    'split-bed12-blocks=lifesci.bio_programs.split_bed12_blocks:main',
    'split-long-chromosomes=lifesci.bio_programs.split_long_chromosomes:main',
    'subtract-bed=lifesci.bio_programs.subtract_bed:main',
]

console_scripts = bio_console_scripts

install_requires = _safe_read_lines("./requirements.txt")

tests_require = [
    'pytest',
    'coverage',
    'pytest-cov',
    'coveralls',
    'pytest-runner',
]

docs_require = [
    'sphinx',
    'sphinx_rtd_theme',
]

pypi_requires = [
    'twine',
    'readme_renderer[md]'
]

all_requires = (
    tests_require + 
    docs_require +
    pypi_requires
)

extras = {
    'test': tests_require,
    'all': all_requires,
    'docs': docs_require,
    'pypi': pypi_requires
}

classifiers=[
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
]

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='lifesci',
    version='0.3.1',
    description="This repo contains python3 life sciences utilities.",
    long_description=readme(),
    long_description_content_type='text/markdown',
    keywords="bioinformatics utilities",
    url="https://github.com/bmmalone/pylifesci",
    author="Brandon Malone",
    author_email="bmmalone@gmail.com",
    license='MIT',
    classifiers=classifiers,
    packages=find_packages(),
    install_requires = install_requires,
    include_package_data=True,
    tests_require=tests_require,
    extras_require=extras,
    entry_points = {
        'console_scripts': console_scripts
    },
    zip_safe=False
)
