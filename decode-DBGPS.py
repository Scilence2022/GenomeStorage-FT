import sys; print('Python %s on %s' % (sys.version, sys.platform))

from utils import *
from test_utils import *
from glass import Glass
from DNAdroplet import DNADroplet
import time
import os.path
import getopt




work_dir = ''
# input_file = work_dir + 'ZDNA.clean.greedy_paths.MP100.k27.c10.r5'
input_file = ''

output_file = work_dir + 'output.jpg'

chunk_size = 30

fountain_seed = 2
index_bytes = 4
crc_bytes = 2
delta = 0.01
c_value = 0.01

passwd = 0

opts,args = getopt.getopt(sys.argv[1:],'-h-i:-o:-p:-d:-n:',
                          ['help','input=','output=',  'pass=', 'chunk_size=', 'chunk_num=', 'seed=',  'index_bytes=', 'ec_bytes='])

usage = 'Usage:\n' + r'      python decode-DBGPS.py -i input_file -o outfile [Options]'
options = 'Options:\n'
options = options + r'      -h, --help                              Show help information' + '\n'
options = options + r'      -i, --input   <input file>              The decoded strands by DBGPS-greedy-path' + '\n'
options = options + r'      -o, --output  <output file>             Output file' + '\n'
options = options + r'      -p, --pass  password                    Password' + '\n'

options = options + r'      -d, --chunk_size  <size>                Chunk size, default = 35 (bytes)' + '\n'
options = options + r'      -n, --chunk_num  <number>               Chunk number, default = 210,000 (for testing only)' + '\n'
options = options + r'      --seed    <seed>                        Fountain random seed, default 1' + '\n'

options = options + r'      --min_index  <initial index>            Initial index, default = 1' + '\n'
options = options + r'      --max_index  <max index>                Max index, default = 20000' + '\n'
options = options + r'      --index_bytes  <number>                 Bytes of index codes, default = 4' + '\n'
options = options + r'      --ec_bytes  <number>                    Bytes of ec codes, default = 2' + '\n'


for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i','--input'):
        input_file = opt_value
    if opt_name in ('-o','--output'):
        output_file = opt_value
    if opt_name in ('-p','--pass'):
        passwd = opt_value
    if opt_name in ('-d','--chunk_size'):
        chunk_size = int(opt_value)
    if opt_name in ('-n', '--chunk_num'):
        chunk_num = int(opt_value)
    if opt_name in ('-s', '--seed'):
        fountain_seed = int(opt_value)
    if opt_name in ('--index_bytes', '--notmatch'):
        index_bytes = int(opt_value)
    if opt_name in ('--crc_bytes', '--notmatch'):
        crc_bytes = int(opt_value)

if not input_file:
    print(usage)
    print(options)
    sys.exit()


# pass_a = bytes(8)
# pass_b = 0b10101110001111011000010100110110.to_bytes(8, 'big')
# pass_c = 0b10011101011111101010101101110101.to_bytes(8, 'big')
# pass_d = 0b00110011010000110010111001000011.to_bytes(8, 'big')

passwd = int(passwd).to_bytes(8, 'big')

start = time.perf_counter()

print('Reading assembled strand sequences')
strands_recovered_multi = read_sim_multi(input_file)
strands_recovered = crc_strands(strands_recovered_multi, passwd)

#assign total chunk sizes
chunk_nums = [2045, 2000, 1914, 2052]
data_sizes = [61343, 59975, 57406, 61537]
chunk_num = chunk_nums[0]
total_bytes = data_sizes[0]
if int.from_bytes(passwd,'big') < 100000000:
    chunk_num = chunk_nums[0]
    total_bytes = data_sizes[0]
if int.from_bytes(passwd,'big') > 2900000000:
    chunk_num = chunk_nums[1]
    total_bytes = data_sizes[1]
if int.from_bytes(passwd,'big') > 2500000000 and int.from_bytes(passwd,'big') < 2900000000:
    chunk_num = chunk_nums[2]
    total_bytes = data_sizes[2]
if int.from_bytes(passwd,'big') < 900000000 and int.from_bytes(passwd,'big') > 100000000:
    chunk_num = chunk_nums[3]
    total_bytes = data_sizes[3]


print('Rebuilding DNA droplets........')
degree_table = get_degrees(chunk_num, chunk_num * 3, fountain_seed, delta, c_value)
cup = Glass(chunk_num)
for aid in strands_recovered:
    adp = DNADroplet(bytes(chunk_size))
    adp.des = True
    adp.sec_key = passwd
    adp.num_of_chunks = chunk_num
    adp.set_head_index(aid)
    adp.set_droplet_from_DNA_CRC(strands_recovered[aid])
    adp.degree = degree_table[aid % len(degree_table)]

    adp.update()
    cup.addDroplet(adp)

print('Decoding by fountain codes .........')
cup.decode()

zero_bytes = chunk_num * chunk_size - total_bytes
if cup.isDone():
    #Fix the file terminal issue
    cup.chunks[-1] = cup.chunks[-1][:-zero_bytes]
    cup.writeToFile(output_file)
    print('Successfully decoded by the Fountain codes. ')
else:
    print('Decoding failed')

cup.writeToFile(output_file)

end = time.perf_counter() - start
print('Decoding time: ', end ='')
print(end)


