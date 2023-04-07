# !/usr/bin/env python
import datetime
import distutils.cmd
import os
import subprocess

from typing import Dict, List

from setuptools import find_packages, setup

if os.path.exists(".git"):
    setup_kwargs: Dict = dict(
        use_scm_version=True,
        setup_requires=["setuptools_scm"],
    )
else:
    setup_kwargs: Dict = dict(version="0+d" + datetime.date.today().strftime("%Y%m%d"))

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

requires: List[str] = [
    "biopython",
    "goatools",
    "joblib",
    "mygene",
    "numpy",
    "matplotlib",
    "mhcnames",
    "pandas",
    "pepdata",
    "pillow",
    "pyensembl",
    "pyllars",
    "pysam",
    "sexpdata",
    "tqdm",
    "weblogo",
]

extras = {
    "dev": [
        "codecov==2.1.12",
        "coverage==6.5.0",
        "docutils==0.19",
        "pylint==2.15.5",
        "pytest-console-scripts==1.3.1",
        "pytest-cov==4.0.0",
        "pytest-datadir==1.4.1",
        "pytest==7.2.0",
        "setuptools_scm==7.1.0",
        "sphinx-argparse==0.3.2",
        "Sphinx==6.1.3",
        'sphinx_rtd_theme',
        'coveralls',
        'pytest-runner',
        'twine',
        'readme_renderer[md]',
    ]
}


classifiers=[
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
]

class DepsCommand(distutils.cmd.Command):
    """
    Custom command to install dependencies only, this
    is intended to be used to take advantage of Docker's
    caching.
    This installs into the user home so don't run it locally, use
    "pip install .[dev]" inside a virtual environment instead.
    """

    description = "Install dependencies"
    user_options = [
        # The format is (long option, short option, description).
        ("env=", None, "Name of environment")
    ]

    def initialize_options(self) -> None:
        self.env = ""

    def finalize_options(self) -> None:
        if self.env:
            assert self.env in ["dev"], "env: {} not available".format(self.env)

    def run(self) -> None:
        command = ["pip", "install", "--user"]
        command.extend(requires)
        if self.env:
            command.extend(extras["dev"])
        if len(command) == 3:
            self.announce("No dependencies to install", level=distutils.log.INFO)
            return
        self.announce(
            "Running command: %s" % " ".join(command), level=distutils.log.INFO
        )
        subprocess.check_call(command)


def readme() -> str:
    try:
        return "\n\n".join([open("README.md").read()])
    except:
        return ""

def description() -> str:
    """ The short description for this project """
    desc = "This repo contains python3 life sciences utilities."
    return desc


setup(
    author="Brandon Malone",
    author_email="bmmalone@gmail.com",
    classifiers=classifiers,
    cmdclass={"deps": DepsCommand},
    description=description(),
    entry_points={"console_scripts": console_scripts},
    extras_require=extras,
    include_package_data=True,
    install_requires=requires,
    keywords="bioinformatics utilities",
    license='MIT',
    long_description=readme(),
    long_description_content_type='text/markdown',
    name='lifesci',
    packages=find_packages(),
    url="https://github.com/bmmalone/pylifesci",
    zip_safe=False,
    **setup_kwargs,
)
