"""
Helper functions for working with STAR. It mostly includes things like
checking if all of the index files exist and adding command line arguments.
"""
import os
import pandas as pd
import shlex
import shutil
import sys

def create_star_tmp(tmp_path:str, tmp_name:str='STAR'):
    """ Ensure the specified directory is ready for use as the temp directory
    for STAR.
    
    In particular, it ensures that tmp_path exists and that tmp_path/tmp_name 
    does not exist (so STAR can create it).

    N.B. This function *does not* use sophisticated file locking mechanisms,
         but it should suffice for "typical" use cases.

    Parameters
    ----------
    tmp_path: string
        the path to the STAR tmp directory NOT including the STAR tmp 
        directory itself

    tmp_name: string
        the name of the STAR tmp directory

    Returns
    -------
    star_temp_path: string
        the complete path to the STAR tmp directory (tmp_path/tmp_name)

    Side effects
    ------------
    tmp_path will be created, if it does not already exist.

    tmp_path/tmp_name will be removed if it already existed.
    """
    # the path to the star directory must exist
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
    
    star_tmp_dir = os.path.join(tmp_path, tmp_name)
    # remove the actual star directory, if it exists
    if os.path.exists(star_tmp_dir):
        shutil.rmtree(star_tmp_dir)

    return star_tmp_dir

def get_star_index_files(path:str, with_sjdb_files:bool=False):
    """ Find the file names necessary for a STAR index rooted at path.

    Parameters
    ----------
    path: string
        the path to the directory containing the STAR index

    with_sjdb_files: bool
        whether to include the splice junction database files in the 
        returned list

    Returns
    -------
    star_index_files: list of strings
        the complete paths to all of the files expected for a STAR index 
        rooted at path (include sjdb/exon files, if specified)
    """
    star_files = [
        "chrLength.txt",
        "chrNameLength.txt",
        "chrName.txt",
        "chrStart.txt",
        "Genome",
        "genomeParameters.txt",
        "SA",
        "SAindex"
    ]

    sjdb_files = [
        "exonGeTrInfo.tab",
        "exonInfo.tab",
        "geneInfo.tab",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab"
    ]

    if with_sjdb_files:
        sf = star_files + sjdb_files
    else:
        sf = star_files

    star_index_files = [ os.path.join(path, star_file) for star_file in sf]
    return star_index_files

def read_star_tr_file(filename:str):
    """ Read a STAR transcript info file into a data frame. 
    
    The file is assumed to contain the number of transcripts on the first row.
    Each row then contains information about each transcript. The columns are 
    assumed to have the following meanings:

        - [ID]      transcript id
        - [S]       transcript (genomic) start
        - [E]       transcript (genomic) end
        - [Emax]    transcript (genomic) end minus "max" (???)
        - [Str]     strand
        - [ExN]     number of exons
        - [ExI]     index of first exon (in STAR exon database)

    N.B. These semantics come from the STAR 2.5.1b source code. They appear 
    stable but could change for new (or old) versions of STAR. The column names
    (in square brackets) follow the internal STAR naming conventions.

    Parameters
    ----------
    filename: string
        the complete path to the STAR transcriptInfo.tab file

    Returns
    -------
    transcript_info: pd.DataFrame
        a data frame containing the transcripts info. The order of the 
        transcripts will be the same as the order in the file.
    """
    column_names = ['ID', 'S', 'E', 'Emax', 'Str', 'ExN', 'ExI']
    transcript_info = pd.read_csv(
        filename, 
        header=None, 
        names=column_names, 
        skiprows=1, 
        sep='\t'
    )

    return transcript_info

# make an attempt to guess the correct "read a gzipped file" command
default_read_files_command = "zcat"
if sys.platform.startswith("darwin"):
    default_read_files_command = "gzcat"

def add_star_options(parser, star_executable:str="STAR", 
        read_files_command:str=default_read_files_command):
    """ Add options to a cmd parser to call STAR.

    N.B. This is primarily intended for use with the rp-bp and b-tea projects.

    Parameters
    ----------
    parser: argparse.ArgumentParser
        The parser to which the options will be added

    star_executable: string
        The name of the star executable. For example, "STARlong" is typically
        used for mapping longer reads, while "STAR" is for most HTS data.

    read_files_command: string
        The system command to read gzipped files.
    """
    star_options = parser.add_argument_group("STAR options")

    star_options.add_argument('--star-executable', help="The name of the STAR "
        "executable", default=star_executable)

    star_options.add_argument('--star-read-files-command', help="The system "
        "command to read gzipped files", default=default_read_files_command)

def get_star_options_string(args):
    """ Extract the flags and options specified for STAR added with 
    add_star_options.

    Parameters
    ---------
    args: argparse.Namespace
        The parsed arguments

    Returns
    -------
    star_options: string
        a string containing STAR options suitable to pass to another command
    """
    args_dict = vars(args)

    star_options = ['star_executable', 'star_read_files_command']

    # create a new dictionary mapping from the flag to the value
    star_options = {'--{}'.format(o.replace('_', '-')) : args_dict[o] 
        for o in star_options if args_dict[o] is not None}

    s = ' '.join("{} {}".format(k,shlex.quote(v)) 
                    for k,v in star_options.items())
    return s

