#!/usr/bin/env python3

"""
This script reads the source file listing sample ids and file names and outputs a list of fastq
file names and paths to be passed to alignment (automatically identifies read pairs based on common
sample id).
N.B. The source file has to be in the same folder as the fastq files or the fastq paths must be absolute
"""

### ---------------------------------------- ###

def listFastqFiles(b_p, src_dt):
    
    # Assign files to reads_list
    reads_list = {}
    
    for sample_id,file in src_dt:
        
        if sample_id in reads_list.keys():
            
            reads_list[sample_id].append(file)
            
            reads_list[sample_id].sort()
    
        else:
            
            reads_list[sample_id] = [file]
    
    # Manage read2 for paired-end and single-read experiments
    for sample_id,files in reads_list.items():
        
        read1, *read2 = files
        
        if len(read2):
            
            read2 = read2[0]
            
        else:
            
            read2 = 'mock.fastq'
            
            reads_list[sample_id].append(read2)
            
            with open('mock.fastq', 'w') as mock_file:

                mock_file.write('Empty')
    
        # Check that read1 and read2 exist and add full path to file name
        if exists(read1) and (exists(read2) if read2 != 'mock.fastq' else True):
            
            read1 = abspath(read1)
            read2 = abspath(read2) if read2 != 'mock.fastq' else read2
        
            reads_list[sample_id] = [read1, read2]
    
        elif exists(f'{b_p}/{read1}') and (exists(f'{b_p}/{read2}') if read2 != 'mock.fastq' else True):
            
            read1 = abspath(f'{b_p}/{read1}')
            read2 = abspath(f'{b_p}/{read2}') if read2 != 'mock.fastq' else read2
        
            reads_list[sample_id] = [read1, read2]
    
        else:
    
            pass
    
    if not len(reads_list):
        
        sys_exit('ERROR: no raw reads files found.')
        
    else:
    
        reads_list = 'SampleID\tFile1\tFile2\n' + '\n'.join(['\t'.join([sample_id, read1, read2]) for sample_id,(read1,read2) in reads_list.items()])

        with open('ReadsList.txt', 'w') as output:

            output.write(reads_list)

### ------------------MAIN------------------ ###

from os import readlink
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

    source_data = [line.split('\t') for line in open(source_file_path, 'r').read().split('\n') if len(line)]

# Get absolute path of the source file (same place where reads should be)
try: # Checking id source file is a symlink

    base_path = abspath(readlink(source_file_path))

except:

    base_path = abspath(source_file_path)

base_path = '/'.join(base_path.split('/')[:-1])

# Making a list of fastq files to be aligned/counted
listFastqFiles(base_path, source_data)
