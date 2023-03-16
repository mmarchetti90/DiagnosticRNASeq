#!/usr/bin/env python3

"""
This script reads the source file listing sample ids and file names and outputs a list of fastq
file names and paths to be passed to alignment (automatically identifies read pairs based on common
sample id).
N.B. The source file has to be in the same folder as the fastq files or the fastq paths must be absolute
"""

### ---------------------------------------- ###

def listFastqFiles(src_dt):
    
    reads_list = ['SampleID\tFile1\tFile2']
    
    for sample_id in set(src_dt.Sample_ID):
        
        read1, *read2 = src_dt.loc[src_dt['Sample_ID'] == sample_id, 'Data_File']
        
        # Manage read2 for paired-end and single-read experiments
        if len(read2):
            
            read2 = read2[0]
            
        else:
            
            read2 = '/mock/path/to/mock.fastq'
        
        # Check that read1 and read2 exist and add full path to file name
        if exists(read1) and (exists(read2) if read2 != '/mock/path/to/mock.fastq' else True):
            
            read1 = abspath(read1)
            read2 = abspath(read2)
        
            reads_list.append(sample_id + '\t' + read1 + '\t' + read2)
    
    if len(reads_list) == 1:
        
        sys_exit('ERROR: no raw reads files found.')
        
    else:
    
        reads_list = '\n'.join(reads_list)

        with open('ReadsList.txt', 'w') as output:

            output.write(reads_list)

### ------------------MAIN------------------ ###

import pandas as pd

from os.path import abspath, exists
from sys import argv
from sys import exit as sys_exit

# Read arguments
if "--source_file" not in argv:
    
    sys_exit('ERROR: you need to specify --source_file.')

elif len(argv) <= argv.index("--source_file") + 1:
    
    sys_exit('ERROR: you need to specify --source_file.')

else:
    
    source_file_path = argv[argv.index("--source_file") + 1]

# Read sample table
if not exists(source_file_path):
    
    sys_exit(f'ERROR: source file {source_file_path} not found.')
    
else:

    source_data = pd.read_csv(source_file_path, sep='\t', header=None, skip_blank_lines=True)
    source_data.columns = ['Sample_ID', 'Data_File']

# Making a list of fastq files to be aligned/counted
listFastqFiles(source_data)
