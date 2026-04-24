import sys; print('Python %s on %s' % (sys.version, sys.platform))

from utils import *
from test_utils import *
from deBruijnGraph import DeBruijnGraph
from glass import Glass
from DNAdroplet import DNADroplet
import time
# import os.path
# import getopt

#####
# chunk_size = 35
# fountain_seed = 2222
# droplet_num = 2000
# double_index = False
# start_index = 1
# work_dir = r''
# input_file = ''
# output_file = ''  # Will be set to input_file + ".fasta" if not specified
# index_bytes = 4
# anchor_bytes = 4
# ec_bytes = 2
# redundancy_rate = 0
# test_num = 0
# test_dropout_rate = 0.05
#####

primerF = 'ATAAGAGGACCTGCCG'  # 5'-->3'
primerF = 'TATCGATGACCTCGAGGA'
primerE = 'CTCGAGGTCATCGATA'  # complement seq of P2

delta = 0.01
c_value = 0.01

input_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/20260408-NGS-Maqiang/X101SC26038058-Z01/X101SC26038058-Z01-J001/01.RawData/A/A_1.fq.gz'

ecc_file_bytes = 109953

# file_type = 'dump_kmers'
# output_file = 'output.ecc'

kmer_size = 27
kmer_cut_off = 10
output_file = r'output_files/' + str(kmer_cut_off) + '22.ECC'
chunk_size = 35

chunk_num = 3142
fountain_seed = 1
index_bytes = 4
crc_bytes = 2
# double_index = False
# both_way = True
index_l = 1
index_u = (chunk_num * 2) + index_l
index_u = 3700
# sec_key = 0b10101010101111011001010100110110
# # sec_key = 0b10101010101111011001010100110111 # wrong key
# sec_key = sec_key.to_bytes(8, byteorder ='big')

# droplet_template = droplet_tp()
def droplet_tp():
    adp = DNADroplet(bytes(chunk_size))
    adp.des = False
    # adp.sec_key = sec_key
    adp.num_of_chunks = chunk_num
    return adp


#
#
# opts,args = getopt.getopt(sys.argv[1:],'-h-i:-o:t:-k:-c:-n:-s:-a-b',
#                           ['help','input=','output=', 'file_type=', 'kmer_size=', 'chunk_size=', 'chunk_num=', 'seed=', 'anchor', 'both_way', 'min_index=',  'max_index=', 'index_bytes=', 'ec_bytes='])
#
# usage = 'Usage:\n' + r'      python decode.py -i input_file -t type_of_seqs -o outfile [Options]'
# options = 'Options:\n'
# options = options + r'      -h, --help                              Show help information' + '\n'
# options = options + r'      -i, --input   <input file>              Input file' + '\n'
# options = options + r'      -t, --file_type   <file type>           Input file type: FastQ, Fasta or Jellyfish dumped k-mers (default)' + '\n'
# options = options + r'      -o, --output  <output file>             Output file' + '\n'
# options = options + r'      -k, --kmer_size  <number>               k-mer size, default = 21 ' + '\n'
# options = options + r'      -c, --chunk_size  <size>                Chunk size, default = 30 (bytes)' + '\n'
# options = options + r'      -n, --chunk_num  <number>               Chunk number, default = 10,741 (for testing only)' + '\n'
# options = options + r'      --cut  <number>                         Cut_off for elimination of low coverage k-mers, default = 0 ' + '\n'
# options = options + r'      -s, --seed    <seed>                    Fountain random seed, default 1' + '\n'
# # options = options + r'      -a, --anchor                            Anchor codes, default On' + '\n'
# # options = options + r'      -b, --both_way                          Both-way search mode, default On' + '\n'
# options = options + r'      --min_index  <initial index>            Initial index, default = 1' + '\n'
# options = options + r'      --max_index  <max index>                Max index, default = 20000' + '\n'
# options = options + r'      --index_bytes  <number>                 Length of index and anchor codes, default = 4 (bytes)' + '\n'
# options = options + r'      --ec_bytes  <number>                    Length of ec codes, default = 4 (bytes)' + '\n'

#
# for opt_name,opt_value in opts:
#     if opt_name in ('-h','--help'):
#         print(usage)
#         print(options)
#         sys.exit()
#     if opt_name in ('-i','--input'):
#         input_file = opt_value
#     if opt_name in ('-o','--output'):
#         output_file = opt_value
#     if opt_name in ('-t','--file_type'):
#         file_type = opt_value
#     if opt_name in ('-k','--kmer_size'):
#         kmer_size = int(opt_value)
#     if opt_name in ('-c','--chunk_size'):
#         chunk_size = int(opt_value)
#     if opt_name in ('-n', '--chunk_num'):
#         chunk_num = int(opt_value)
#     if opt_name in ('-s', '--seed'):
#         f_seed = int(opt_value)
#     # if opt_name in ('-a', '--anchor'):
#     #     double_index = True
#     if opt_name in ('--min_index', '--notmatch'):
#         index_l = int(opt_value)
#
#     if opt_name in ('--max_index', '--notmatch'):
#         index_u = int(opt_value)
#
#     if opt_name in ('--index_bytes', '--notmatch'):
#         index_bytes = int(opt_value)
#
#     if opt_name in ('--crc_bytes', '--notmatch'):
#         crc_bytes = int(opt_value)


start = time.perf_counter()

deG = DeBruijnGraph()
deG.set_kmer_len(kmer_size)
deG.primerF = primerF
deG.primerE = primerE
deG.min_cov = kmer_cut_off + 1
deG.max_path_num = 1000000

print('\nCounting k-mers ......')
deG.count_file(input_file)

# deG.kmer_len = kmer_size
# deG.veri_kmers = False
# deG.max_path_num = 10000

# if file_type == 'FastQ' or file_type == 'fastq':
#     deG.open_fastq(input_file)
# else:
#     if file_type == 'fasta' or file_type == 'Fasta':
#         deG.open_fasta(input_file)
#     else:
#         deG.open_dump(input_file)

print('\nReconstructing DNA droplets ......')

degree_table = get_degrees(chunk_num, chunk_num * 3, fountain_seed, delta, c_value)

a = time.perf_counter()
cup = Glass(chunk_num)

droplet_template = droplet_tp()
for index in range(index_l, index_u+1):

    fd_droplets = deG.find_droplets(index, droplet_template)
    #print(len(fd_droplets), end=" ")
    if len(fd_droplets) == 1:
        fd_droplets[0].degree = degree_table[index % len(degree_table)]
        fd_droplets[0].update()
        cup.addDroplet(fd_droplets[0])
        if len(cup.droplets) % 100 < 1:
            print("Reconstructed droplet number: " + str(len(cup.droplets)))
    else:
        if len(fd_droplets) > 1:
            print("Multiple droplets found with index " + str(index))

# hashToFile(deG.path_sta, output_file + "cut_sta")

print('Decoding by fountain codes .........')

cup.decode()
if cup.isDone():
    cup.writeToFile(output_file, ecc_file_bytes)
else:
    print("Failed to decode the original information.")

end = time.perf_counter() - start
print('Decoding time: ', end ='')
print(end)



