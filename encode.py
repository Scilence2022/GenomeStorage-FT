import sys; print('Python %s on %s' % (sys.version, sys.platform))

import random
from utils import *
from DNAfountain import DNAFountain
from test_utils import *

import copy
import time
import os.path
import sys
import getopt


chunk_size = 35
fountain_seed = 1
droplet_num = 2000
double_index = False
start_index = 1
work_dir = r''
input_file = ''
output_file = ''  # Will be set to input_file + ".fasta" if not specified
index_bytes = 4
anchor_bytes = 4
ec_bytes = 2
redundancy_rate = 0
test_num = 0


p1 = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
p2 = 'CTGACACTGATGCATCCG'  # complement seq of P2


opts,args = getopt.getopt(sys.argv[1:],'-h-i:-o:-c:-n:-s:-r:-d-l:',
                          ['help','input=','output=',  'chunk_size=', 'redundancy_rate=', 'droplet_num=', 'seed=', 'initial_index=', 'index_bytes=', 'ec_bytes=', 'test_num='])

usage = 'Usage:\n' + r'      python encode.py -i input_file -n number_of_droplets -o output.fasta [Options]'
options = 'Options:\n'
options = options + r'      -h, --help                              Show help information' + '\n'
options = options + r'      -i, --input   <input file>              Input file' + '\n'
options = options + r'      -o, --output  <output file>             Output file' + '\n'
options = options + r'      -r, --redundancy_rate  <number>         Redundancy rate, default 0 ' + '\n'
options = options + r'      -n, --droplet_num  <number>             Number of droplets, default 210,000 ' + '\n'
options = options + r'      -c, --chunk_size  <size>                Chunk size, default 35 (bytes)' + '\n'
options = options + r'      -s, --seed    <seed>                    Fountain random seed, default 1' + '\n'
options = options + r'      -l, --initial_index  <initial index>    Initial index, default 1' + '\n'
options = options + r'      --index_bytes  <number>                 Length of index and anchor codes, default 4 (bytes)' + '\n'
options = options + r'      --ec_bytes  <number>                    Length of ec codes, default 2 (bytes)' + '\n'
options = options + r'      --test_num  <number>                    Number of test runs, default 0' + '\n'


for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i','--input'):
        input_file = opt_value
        if not output_file:  # Set default output name if not already specified
            output_file = input_file + ".fasta"
    if opt_name in ('-o','--output'):
        output_file = opt_value
    if opt_name in ('-c','--chunk_size'):
        chunk_size = int(opt_value)
    if opt_name in ('-n', '--droplet_num'):
        droplet_num = int(opt_value)

    if opt_name in ('-r', '--redundancy_rate'):
        redundancy_rate = float(opt_value)

    if opt_name in ('-s', '--seed'):
        fountain_seed = int(opt_value)

    if opt_name in ('-l', '--low_index'):
        start_index = int(opt_value)
    if opt_name in ('--index_bytes', '--notmatch'):
        index_bytes = int(opt_value)
        anchor_bytes = index_bytes
    if opt_name in ('--ec_bytes', '--notmatch'):
        ec_bytes = int(opt_value)

    if opt_name in ('--test_num', '--notmatch'):
        test_num = int(opt_value)

start = time.perf_counter()

file1 = open(work_dir + input_file, 'rb')
filebytes1 = file1.read()
file1.close()


fdna1 = DNAFountain(filebytes1, chunk_size, start_index, fountain_seed, index_bytes*4, anchor_bytes*4, ec_bytes*4)
fdna1.degree_table_folds = 5
fdna1.gen_degrees()

print('Chunk number: ', end='')
print(fdna1.num_of_chunks)

if redundancy_rate > 0:
    droplet_num = int(fdna1.num_of_chunks * (1 + redundancy_rate)) + 1

print('Generating droplets' + str(droplet_num) + '............')
total_size = droplet_num #int(fdna1.num_of_chunks * 1.15)
core_size = int(total_size*0.95)

print(droplet_num)
print(core_size)
droplet_all = get_droplets(total_size, fdna1)

if test_num > 0:
    print('\nGenerating DNA data droplets and run strand dropout test ......')
    failed_gts = []


    suc_num_default = 0
    for i in range(0, test_num):
        droplet_sample = []
        # droplet_sample[i] = random.sample(range(0, total_size), core_size)
        random.seed(i + 1)
        for j in random.sample(range(0, total_size-1), core_size):
            droplet_sample.append(copy.deepcopy(droplet_all[j]))

        if test_droplets(droplet_sample, fdna1):
            suc_num_default = suc_num_default + 1

    print("With a dropout rate of 6%, " + str(suc_num_default) + " succeed in a total of " + str(test_num) + " runs.")
    print('\nGenerating DNA data droplets and run strand dropout test ......')


write_log_file(output_file, droplet_all, double_index)
write_fasta_file(output_file, droplet_all, p1, p2, double_index)


print('The DNA sequences in fasta format: ', end='')
print(output_file)

print('Details about the encoding: ', end='')
print(output_file + '.log')

end = time.perf_counter() - start
print('Encoding time: ', end ='')
print(end)
