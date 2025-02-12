import random
import math
import numpy as np
from random import choices
from pyDes import des, CBC, PAD_PKCS5, ECB


def	charN(str, N):
	if N < len(str):
		return str[N]
	return 'X'


def oxor(str1, str2):
	length = max(len(str1), len(str2))
	return ''.join(chr(ord(charN(str1, i)) ^ ord(charN(str2, i))) for i in range(length))

def xor(str1, str2):
	length = max(len(str1), len(str2))
	if len(str1) > len(str2):
		str2 = str2 + bytes(len(str1) - len(str2))
	if len(str2) > len(str1):
		str1 = str1 + bytes(len(str2) - len(str1))
	allBytes = b''
	i = 0
	while i < length:
		allBytes = allBytes + bytes([str1[i] ^ str2[i]])
		i =  i + 1
	return allBytes


def randChunkNums(degree, num_chunks):
	size = random.randint(1,max(5, int(num_chunks*3/5)))
	return random.sample(range(num_chunks), size)


def get_degrees(N, k, symbol_index=-1, robust_failure_probability=0.1, c_value=0.2, distribution_name='robust'):
	#print(distribution_name)
	if distribution_name == "ideal":
		probabilities = ideal_distribution(N)
	elif distribution_name == "robust":
		probabilities = robust_distribution(N, c_value, robust_failure_probability)
	else:
		probabilities = None


	print("get_degrees() probabilities: " )
	print(probabilities)
	print(len(probabilities))

	population = list(range(0, N))
	if symbol_index > 0:
		random.seed(symbol_index)
	#Not initialized yet
	else:
		return -1
	return choices(population, probabilities, k=k)


def generate_chunk_nums(blocks_quantity, degree, symbol_index):
	random.seed(symbol_index)
	indexes = random.sample(range(blocks_quantity), degree[0])
	return indexes




def ideal_distribution(K, epsilon=1e-6):

	probabilities = [0, 10/ K]
	probabilities += [1 / (k * (k - 1)) for k in range(2, K)]
	probabilities_sum = np.sum(probabilities)
	probabilities /= probabilities_sum

	probabilities_sum = np.sum(probabilities)
	assert probabilities_sum >= 1 - epsilon and probabilities_sum <= 1 + epsilon, "The ideal distribution should be standardized"
	return probabilities

def robust_distribution(K, c_value = 0.5, robust_failure_probability=0.1, epsilon=1e-6):
    """
    Returns a robust distribution used in Luby Fountain Codes, incorporating a robustness factor and spike.

    Parameters:
    - K: The total number of symbols
    - c_value: A scaling constant affecting the robustness
    - robust_failure_probability: Probability of failure used in robustness calculation
    - epsilon: Precision for checking sum normalization
    """
    S = c_value * math.log(K / robust_failure_probability) * math.sqrt(K)
    #print(S)
    # M = round(K/S)
    M = round(K / S)
    #print(M)
    extra_proba = [0] + [1/(i * M) for i in range(1, M-1)]
    extra_proba += [S * math.log(S / robust_failure_probability) / K]  # Spike at M
    extra_proba += [0 for k in range(M, K)]

    probabilities = np.add(extra_proba, ideal_distribution(K))
    #print(np.sum(probabilities))
    probabilities /= np.sum(probabilities)
    probabilities_sum = np.sum(probabilities)
    assert probabilities_sum >= 1 - epsilon and probabilities_sum <= 1 + epsilon, "The robust distribution should be standardized"
    return probabilities


def bytesToDNA(manyBytes):
	dnastr = ''
	for aByte in manyBytes:
		dnastr = dnastr + byteToDNA(aByte)
	return dnastr

def DNAToBytes(dnaStr):
	i = 0
	nBytes=b''
	while i < len(dnaStr):

		nBytes = nBytes + simDNAToByte(dnaStr[i:i+4])
		i = i + 4
	return nBytes

def byteToDNA(abyte):
	#define the converter of four bits to DNA chars
	converter = ('AT', 'AG', 'AC', 'AA',  'TA', 'TC', 'TG', 'TT', 'GG', 'GA', 'GT', 'GC',   'CC', 'CT', 'CA', 'CG')
	octNum = abyte
	#octNum = ord(abyte)
	dnastr = converter[octNum % 16]
	octNum = octNum//16
	dnastr = converter[octNum % 16] + dnastr
	return dnastr



def simDNAToByte(dnaStr):
	converter = ('AT', 'AG', 'AC', 'AA', 'TA', 'TC', 'TG', 'TT', 'GG', 'GA', 'GT', 'GC', 'CC', 'CT', 'CA', 'CG')
	bytesConverter = {}
	i = 0
	while i < 16:
		j = 0
		while j < 16:
			fourBases = converter[i] + converter[j]
			bytesConverter[fourBases] = bytes([i*16 + j])
			j = j + 1
		i = i + 1
	return bytesConverter[dnaStr[0:2] + dnaStr[2:4]]

def DNAToByte(dnaStr):
	converter = ('AT', 'AG', 'AC', 'AA', 'TA', 'TC', 'TG', 'TT', 'GG', 'GA', 'GT', 'GC', 'CC', 'CT', 'CA', 'CG')
	bytesConverter = {}
	i = 0
	while i < 16:
		j = 0
		while j < 16:
			fourBases = converter[i] + converter[j]
			bytesConverter[fourBases] = bytes([i*16 + j])
			j = j + 1
		i = i + 1
	return bytesConverter[dnaStr[0:2] + dnaStr[3:5]]




def compressDNA(dnaStr):
	#return dnaStr
	converter = ('A', 'T', 'G', 'C')
	bytesConverter = {}
	i = 0
	while i < 4:
		j = 0
		while j < 4:
			k = 0
			while k < 4:
				m = 0
				while m < 4:
					fourBases = converter[i] + converter[j] + converter[k] + converter[m]
					bytesConverter[fourBases] = bytes([i*4*4*4 + j*4*4 + k*4 + m])
					m = m + 1
				k = k + 1
			j = j + 1
		i = i + 1
	i = 0
	nBytes = b''
	while i < len(dnaStr):
		nBytes = nBytes + bytesConverter[dnaStr[i:i+4]]
		i = i + 4
	return dnaStr

def depressDNA(manyBytes):
	#define the converter of four bits to DNA chars
	converter = ('A', 'T', 'G', 'C')
	dnastr = ''
	for aByte in manyBytes:
		fourBase = ''
		fourBase = converter[aByte % 4] + fourBase
		aByte = aByte // 4
		fourBase = converter[aByte % 4] + fourBase
		aByte = aByte // 4
		fourBase = converter[aByte % 4] + fourBase
		aByte = aByte // 4
		fourBase = converter[aByte % 4] + fourBase
		dnastr = dnastr + fourBase
	return dnastr



#2019-5-14
def calc_gc(dnastr):
	dnastr = dnastr.upper()
	gc_num = dnastr.count("G") + dnastr.count("C")
	return gc_num/len(dnastr)

def max_homo_len(dnastr):
	dnastr = dnastr.upper()
	max_len = 0
	pre_len = 0
	last_c = ''
	for c in dnastr:
		if c == last_c:
			pre_len = pre_len + 1
		else:
			if pre_len > max_len:
				max_len = pre_len
			pre_len = 1
		last_c = c
	if pre_len > max_len:
		max_len = pre_len
	return max_len

#2019-05-15
#2020-03-08
#def check_dna(dnastr, min_gc=0.45, max_gc=0.55, max_homo=5):
def check_dna(dnastr, min_gc=0.4, max_gc=0.6, max_homo=7):
	# print(dnastr)
	gc_rate = calc_gc(dnastr)
	if gc_rate > max_gc:
		return False
	if gc_rate < min_gc:
		return False
	homo_poly_len = max_homo_len(dnastr)
	if homo_poly_len > max_homo:
		return False
	return True

#2019-05-16
def rev_seq(dnastr):
	complement = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
	dnastr = dnastr.upper()
	rev_dna = ""
	for i in dnastr:
		rev_dna += complement[i]
	rev_dna = rev_dna[::-1]
	return rev_dna

#20190625
def num_randint(a,b,num):
	i = 0
	int_nums = []
	while i < num:
		int_nums.append(random.randint(a, b))
		i += 1
	return int_nums


def randomATGC():
	a = "ATGC"
	ar = random.randint(0, 3)
	return a[ar:ar + 1]


def randomDNA(len):
	i = 1
	dna = ""
	while(i <=len):
		dna= dna + randomATGC()
		i = i + 1

	return dna



def DNA_complement(sequence):
	sequence = sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	return sequence.upper()

def DNA_reverse(sequence):
	sequence = sequence.upper()
	return sequence[::-1]

def DNA_rev_complement(sequence):
	sequence = DNA_complement(sequence)
	return DNA_reverse(sequence)

def expntl(L):
	"""
    negative exponential distribution
    return a double random number, L is the mean value
    """
	u = random.randint(0,2**L)
	return int(np.math.log(u)/math.log(2))



def file_to_array(file):
	f = open(file)
	arr = []
	line = f.readline()
	while line.strip():
		arr.append(line.strip())
		line = f.readline()
	return arr


def kmers_of_str(str, kmer_len=21,step_len=1):
	kmers = {}
	i = 0
	if len(str) >= kmer_len:
		i = 0
		kmstr = ''
		while i <= len(str) - kmer_len:
			kmstr = str[i:i + kmer_len]
			kmers[kmstr] = 1
			i = i + step_len
	return kmers.keys()

#2020-05-16
def kmers_of_position(str, kmer_len,pos=0):
	if len(str)-pos >= kmer_len:
		return str[pos:pos+kmer_len]




#2020-05-03
def kmers_in_dict(kmers, dict):
	for kmer in kmers:
		if kmer not in dict:
			return False
	return True

def any_kmers_in_dict(kmers, dict):
	for kmer in kmers:
		if kmer in dict:
			return True
	return False

def read_file(file):
	file1 = open(file, 'rb')
	filebytes = file1.read()
	file1.close()
	return filebytes


def des_en(s, sec_key):
	iv = sec_key
	k = des(sec_key, CBC, iv, pad=None, padmode=PAD_PKCS5)
	return k.encrypt(s)

def des_de(s, sec_key):
	iv = sec_key
	k = des(sec_key, CBC, iv, pad=None, padmode=PAD_PKCS5)
	return k.decrypt(s)

def write_log_file(file_path, droplet_all, double_index=False):
    """
    Write detailed droplet information to a log file
    
    Args:
        file_path: Path to the log file
        droplet_all: List of droplets to write
        double_index: Boolean flag for index format
    """
    with open(file_path + '.log', 'tw') as file2:
        for dps in droplet_all:
            file2.write(str(dps.head_index))
            file2.write('\t')
            file2.write(str(dps.data))
            file2.write('\t')
            if double_index:
                file2.write(dps.to_DNA_CRC())
            else:
                file2.write(dps.to_DNA_CRC_sIndex())
            file2.write('\t')
            file2.write(str(dps.degree))
            file2.write('\t')
            file2.write(str(dps.get_chunk_nums()))
            file2.write('\t')
            file2.write(str(dps.tail_index))
            file2.write('\n')

def write_fasta_file(file_path, droplet_all, p1, p2, double_index=False):
    """
    Write droplets to a FASTA format file
    
    Args:
        file_path: Path to the fasta file
        droplet_all: List of droplets to write
        p1: First primer sequence
        p2: Second primer sequence
        double_index: Boolean flag for index format
    """
    with open(file_path, 'tw') as file2:
        for dps in droplet_all:
            file2.write(">")
            file2.write(str(dps.head_index))
            file2.write("\n")
            file2.write(p1)
            if double_index:
                file2.write(dps.to_DNA_CRC())
            else:
                file2.write(dps.to_DNA_CRC_sIndex())
            file2.write(p2)
            file2.write('\n')

def reversible_hash_32(x: int) -> int:
    """
    A reversible 32-bit hash function using xor and bit rotation
    
    Args:
        x: 32-bit integer input
    Returns:
        32-bit hashed integer
    """
    x = ((x >> 16) ^ x) * 0x45d9f3b
    x = ((x >> 16) ^ x) * 0x45d9f3b
    x = (x >> 16) ^ x
    return x & 0xFFFFFFFF

def inverse_hash_32(y: int) -> int:
    """
    Inverse function of reversible_hash_32
    
    Args:
        y: 32-bit hashed integer
    Returns:
        Original 32-bit integer
    """
    y = y & 0xFFFFFFFF
    y = ((y >> 16) ^ y) * 0x119de1f3
    y = ((y >> 16) ^ y) * 0x119de1f3
    y = (y >> 16) ^ y
    return y & 0xFFFFFFFF


def jenkins_hash(anum, bit_length=32):
    data = int.to_bytes(anum, int(bit_length / 8), 'big', signed=False)
    if bit_length <= 0:
        raise ValueError("Invalid bit length. Must be greater than 0")
    hash_value = 0
    for byte in data:
        hash_value += byte
        hash_value += (hash_value << 10)
        hash_value ^= (hash_value >> 6)
        hash_value += (hash_value << 3)
        hash_value ^= (hash_value >> 11)
        hash_value += (hash_value << 15)
    hash_mask = (1 << bit_length) - 1
    return hash_value & hash_mask