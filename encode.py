import random
from utils import *
from DNAdroplet import DNADroplet
#from fountain import Fountain
from DNAfountain import DNAFountain
from glass import Glass
from crc16pure import *
from test_utils import *
import copy
import os


data_block_length = 30
test_num = 10
fountain_seed = 2

work_dir = r'./input_files'  + '/'

print('\nReading input files ...')

# fileA = open(work_dir + r'Babel-A.jpg', 'rb')
file_tju = open(work_dir + r'tju.jpg', 'rb')
fileA = open(work_dir + r'BCA.jpg', 'rb')
fileB = open(work_dir + r'HECHAIN.jpg', 'rb')
fileC = open(work_dir + r'ZONFF.jpg', 'rb')

filebytes_tju = file_tju.read()
file_tju.close()
filebytesA = fileA.read()
file_tju.close()
filebytesB = fileB.read()
fileB.close()
filebytesC = fileC.read()
fileC.close()

pass_default = bytes(8)
pass_a = 0b10101110001111011000010100110110
pass_b = 0b10011101011111101010101101110101
pass_c = 0b00110011010000110010111001000011


fountain_init_index = 3111111
ft_default = DNAFountain(filebytes_tju, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
ft_default.fix_bytes()
ft_default.des = True
ft_default.sec_key = pass_default


ft_a = DNAFountain(filebytesA, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
ft_a.des = True
ft_a.sec_key = pass_a.to_bytes(8, byteorder ='big')

ft_b = DNAFountain(filebytesB, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
ft_b.des = True
ft_b.sec_key = pass_b.to_bytes(8, byteorder ='big')

ft_c = DNAFountain(filebytesC, data_block_length, fountain_init_index, fountain_seed, 16, 16, 8)
ft_c.des = True
ft_c.sec_key = pass_c.to_bytes(8, byteorder ='big')


print('\nGenerating DNA data droplets and run strand dropout test ......')
failed_gts = []

total_size = 3000 #int(fdna1.num_of_chunks * 1.15)
core_size = int(total_size*0.96)

droplet_four_pics = []

droplet_all = get_droplets(total_size, ft_default)
droplet_four_pics.append(droplet_all)
suc_num_default = 0
for i in range(0, test_num):

    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)
    for j in random.sample(range(0, total_size), core_size):
        droplet_sample.append(copy.deepcopy(droplet_all[j]))

    if test_droplets(droplet_sample, ft_default):
        suc_num_default = suc_num_default + 1

print("With a dropout rate of 4%, " + str(suc_num_default) + " succeed in a total of " + str(test_num) + " runs.")
print('\nGenerating DNA data droplets and run strand dropout test ......')

droplet_all = get_droplets(total_size, ft_a)
droplet_four_pics.append(droplet_all)

suc_num_a = 0
for i in range(0, test_num):

    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)
    for j in random.sample(range(0, total_size), core_size):
        droplet_sample.append(copy.deepcopy(droplet_all[j]))

    if test_droplets(droplet_sample, ft_a):
        suc_num_a = suc_num_a + 1

print("With a dropout rate of 4%, " + str(suc_num_a) + " succeed in a total of " + str(test_num) + " runs.")
print('\nGenerating DNA data droplets and run strand dropout test ......')

droplet_all = get_droplets(total_size, ft_b)
droplet_four_pics.append(droplet_all)
suc_num_b = 0
for i in range(0, test_num):

    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)
    for j in random.sample(range(0, total_size), core_size):
        droplet_sample.append(copy.deepcopy(droplet_all[j]))

    if test_droplets(droplet_sample, ft_b):
        suc_num_b = suc_num_b + 1

print("With a dropout rate of 4%, " + str(suc_num_b) + " succeed in a total of " + str(test_num) + " runs.")
print('\nGenerating DNA data droplets and run strand dropout test ......')

droplet_all = get_droplets(total_size, ft_c)
droplet_four_pics.append(droplet_all)
suc_num_c = 0
for i in range(0, test_num):

    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)
    for j in random.sample(range(0, total_size), core_size):
        droplet_sample.append(copy.deepcopy(droplet_all[j]))

    if test_droplets(droplet_sample, ft_c):
        suc_num_c = suc_num_c + 1

print("With a dropout rate of 4%, " + str(suc_num_c) + " succeed in a total of " + str(test_num) + " runs.")


print('\nOutputting strand sequences and details ......')

p1 = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
p2 = 'CTGACACTGATGCATCCG'  # complement seq of P2

for i in range(0, 4):
    file2 = open(work_dir + r'babel_v2.DNAs.tab.sim.' + str(i), 'tw')
    #file2.write('Head Index\tData\tDNA\tDNA-Primers\tDegree\tChunk Nums\tTail Index\n')
    for dps in droplet_four_pics[i]:
        file2.write(str(dps.head_index))
        file2.write('\t')
        file2.write(p1)
        file2.write(dps.to_DNA_CRC_sIndex())
        file2.write(p2)
        file2.write('\n')
    file2.close()


for i in range(0, 4):
    file2 = open(work_dir + r'babel_v2.DNAs.tab.rich.' + str(i), 'tw')
    file2.write('Head Index\tData\tDNA\tDNA-Primers\tDegree\tChunk Nums\tTail Index\n')
    for dps in droplet_four_pics[i]:
        file2.write(str(dps.head_index))
        file2.write('\t')
        file2.write(str(dps.data))
        file2.write('\t')
        file2.write(dps.to_DNA_CRC_sIndex())
        file2.write('\t')
        file2.write(p1)
        file2.write(dps.to_DNA_CRC_sIndex())
        file2.write(p2)
        file2.write('\t')
        file2.write(str(dps.degree))
        file2.write('\t')
        file2.write(str(dps.get_chunk_nums()))
        file2.write('\t')
        file2.write(str(dps.tail_index))
        file2.write('\n')
    file2.close()




os.system('mv input_files/babel_v2.DNAs.tab.sim.0 input_files/tju.jpg.strands')
os.system('mv input_files/babel_v2.DNAs.tab.rich.0 input_files/tju.jpg.details')

os.system('mv input_files/babel_v2.DNAs.tab.sim.1 input_files/BCA.jpg.strands')
os.system('mv input_files/babel_v2.DNAs.tab.rich.1 input_files/BCA.jpg.details')

os.system('mv input_files/babel_v2.DNAs.tab.sim.2 input_files/HECHAIN.jpg.strands')
os.system('mv input_files/babel_v2.DNAs.tab.rich.2 input_files/HECHAIN.jpg.details')

os.system('mv input_files/babel_v2.DNAs.tab.sim.3 input_files/ZONFF.jpg.strands')
os.system('mv input_files/babel_v2.DNAs.tab.rich.3 input_files/ZONFF.jpg.details')













