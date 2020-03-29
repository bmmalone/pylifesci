#! /usr/bin/env python3

import argparse
import ftplib
import os
import pandas as pd
import sys

import misc.parallel as parallel
import misc.shell_utils as shell_utils
import misc.slurm as slurm

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

default_sep = '\t'
default_accession_field = "Run"
default_paired_field = "LibraryLayout"
default_paired_values = ["PAIRED"]
default_num_cpus = 1
default_num_downloads_per_connection = 1

source_choices = ["ebi", "ncbi"]
default_source = "ebi"
servers = {
    'ebi': "ftp.sra.ebi.ac.uk",
    'ncbi': "ftp-trace.ncbi.nih.gov"
}

def get_local_sra_file(accession_num, args):
    filename = "{}.sra".format(accession_num)
    local_path = os.path.join(args.outdir, filename)
    return local_path

def get_local_fastq_files(row, args):
    local_fastq_files = []
    accession_num = row[args.accession_field]

    is_paired = row[args.paired_field] in args.paired_values

    if is_paired:
        filename = "{}_1.fastq.gz".format(accession_num)
        local_path = os.path.join(args.outdir, filename)
        local_fastq_files.append(local_path)

        filename = "{}_2.fastq.gz".format(accession_num)
        local_path = os.path.join(args.outdir, filename)
        local_fastq_files.append(local_path)
    else:
        filename = "{}.fastq.gz".format(accession_num)
        local_path = os.path.join(args.outdir, filename)
        local_fastq_files.append(local_path)

    return local_fastq_files

def get_ncbi_path(run_accession):
    # the directory structure for ncbi is here:
    #     http://www.ncbi.nlm.nih.gov/books/NBK158899/
    nbci_server = "ftp-trace.ncbi.nih.gov"
    
    ncbi_path_1 = "/sra/sra-instant/reads/ByRun/sra/SRR/"
    ncbi_path_2 = run_accession[:6]
    ncbi_path_3 = run_accession
    ncbi_path = os.path.join(ncbi_path_1, ncbi_path_2, ncbi_path_3)
    
    ncbi_file = "{}.sra".format(run_accession)
    
    return nbci_server, ncbi_path, ncbi_file

def get_ebi_path(run_accession):
    # the directory structure for ebi is here:
    #     http://www.ebi.ac.uk/ena/browse/read-download

    # ftp://ftp.sra.ebi.ac.uk/vol1/fastq/<dir1>[/<dir2>]/<run accession>
    ebi_server = "ftp.sra.ebi.ac.uk"
    
    ebi_path_1 = "/vol1/srr/"

    # <dir1> is the first 6 letters and numbers of the run accession
    ebi_part_2 = run_accession[:6]

    # <dir2> does not exist if the run accession has six digits
    accenssion_len = len(run_accession) - 3
    if accenssion_len == 6:
        ebi_part_3 = ""
    elif accenssion_len == 7:
        #<dir2> is 00 + the last digit of the run accession
        ebi_part_3 = "00" + run_accession[-1:]
    elif accenssion_len == 8:
        # <dir2> is 0 + the last two digits of the run accession.
        ebi_part_3 = "0" + run_accession[-2:]
    elif accenssion_len == 9:
        # <dir2> is the last three digits of the run accession
        ebi_part_3 = run_accession[-3]
        
    ebi_path = os.path.join(ebi_path_1, ebi_part_2, ebi_part_3)

    ebi_file = run_accession
    
    return ebi_server, ebi_path, ebi_file

def download_sra_file(row, args):

    get_path = get_ebi_path
    if args.source == "ncbi":
        get_path = get_ncbi_path

    run_accession = row[args.accession_field]
    server, path, filename = get_path(run_accession)
    local_sra_file = get_local_sra_file(run_accession, args)
    local_fastq_files = get_local_fastq_files(row, args)

    # check if we have already downloaded these
    sra_file_exists = os.path.exists(local_sra_file)
    fq_files_exist = all(os.path.exists(f) for f in local_fastq_files)
    files_exist = sra_file_exists or fq_files_exist
    
    if not args.overwrite and files_exist:
        msg = "The .sra or .fastq.gz already exist. Skipping: {}".format(local_sra_file)
        logger.info(msg)
        return
    
    msg = "Downloading {} : {}".format(path, filename)
    logger.info(msg)

    try:
        remote_path = "http://" + server + os.path.join(path, filename)
        shell_utils.download_file(remote_path, local_filename=local_sra_file)
    except:
        msg = "There was a problem downloading: {}".format(filename)
        logger.warning(msg)

def extract_sra_file(row, args):

    get_path = get_ebi_path
    if args.source == "ncbi":
        get_path = get_ncbi_path

    run_accession = row[args.accession_field]
    server, path, filename = get_path(run_accession)
    local_sra_file = get_local_sra_file(run_accession, args)
    local_fastq_files = get_local_fastq_files(row, args)

    # make sure the sra file exists
    if not os.path.exists(local_sra_file):
        # maybe we have already extracted it...
        if all(os.path.exists(f) for f in local_fastq_files):
            msg = ("The .sra file has already been extracted: {}"
                .format(local_sra_file))
            logger.info(msg)
            return

        else:
            msg = "Could not find the .sra file: {}".format(local_sra_file)
            logger.warning(msg)
            return

    msg = "Extracting .sra file: {}".format(local_sra_file)
    logger.info(msg)

    # first, check if it is paired-end or not
    split_files_str = ""

    if args.paired_field is not None:
        if row[args.paired_field] in args.paired_values:
            split_files_str = "--split-files"
    


    outdir_str = "--outdir {}".format(args.outdir)
    cmd = ("fastq-dump {} {} {} --gzip --skip-technical --readids --dumpbase "
        "--clip --origfmt".format(outdir_str, split_files_str, local_sra_file))

    try:
        ret = shell_utils.check_call(cmd)
    except:
        # just set the return code to something we will catch later
        ret = -1

    if ret == 0:
        msg = "Removing .sra file: {}".format(local_sra_file)
        logger.info(msg)

        os.remove(local_sra_file)
    else:
        msg = ("Problem extracting .sra file. It will be removed: {}".
            format(local_sra_file))
        logger.warning(msg)

        os.remove(local_sra_file)

def process_files(srr, args):

    for i, (index, row) in enumerate(srr.iterrows()):

        download_sra_file(row, args)
        extract_sra_file(row, args)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script downloads short read archive runs (i.e., SRR) files "
        "over ftp. It only requires the run number. It also converts the files from "
        "the .sra format to .fastq.gz files. It then deletes the .sra file.")

    parser.add_argument('srr', help="A csv file containing the SRR accessions to "
        "download. Optionally, it can also include whether the samples are paired-"
        "end or not.")
    parser.add_argument('outdir', help="The location for the fastq.gz files")

    parser.add_argument('-a', '--accession-field', help="The name of the column "
        "containing the SRR identifiers", default=default_accession_field)

    parser.add_argument('-p', '--paired-field', help="The name of the column "
        "indicating whether the sample is paired-end", default=default_paired_field)

    parser.add_argument('-v', '--paired-values', help="The exact string values in "
        "the paired-field which indicate the sample is paired-end", nargs="*",
        default=default_paired_values)

    parser.add_argument('-s', '--source', help="The server from which the files "
        "will be downloaded", choices=source_choices, default=default_source)

    parser.add_argument('--overwrite', help="If this flag is given, then existing "
        "files will be re-downloaded. Otherwise, if either the .sra or .fastq.gz "
        "file already exists, then the sra file will not be downloaded.",
        action='store_true')

    parser.add_argument('--num-downloads-per-connection', help="The number of "
        "files to download with each open connection. Each connections will be "
        "closed and re-opened after this many files are downloaded.", type=int,
        default=default_num_downloads_per_connection)

    parser.add_argument('--sep', help="The separator in the SRR file", 
        default=default_sep)
     
    slurm.add_sbatch_options(parser)
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    programs = [
        'fastq-dump'
    ]
    shell_utils.check_programs_exist(programs)
    
    # check if we want to use slurm
    if args.use_slurm:
        msg = ("The --use-slurm option was given, so sbatch will now be used "
            "to submit to slurm.")
        logger.warning(msg)

        cmd = ' '.join(sys.argv)

        slurm.check_sbatch(cmd, args=args)

        # and quit!
        return

    msg = "Reading SRR list"
    logger.info(msg)

    srr = pd.read_csv(args.srr, sep=args.sep)

    parallel.apply_parallel_split(srr, args.num_cpus, process_files, args)

if __name__ == '__main__':
    main()
