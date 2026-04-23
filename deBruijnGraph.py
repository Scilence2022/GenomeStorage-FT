from utils import *
from DNAdroplet import DNADroplet
import khmer


class DeBruijnGraph:
    def __init__(self, kmer_len=21):
        self.kmer_len = kmer_len
        self.counts = khmer.Counttable(kmer_len, 1e8, 4)

        self.p_kmer = {}
        self.cRule = {}
        self.next_bases = {}
        self.min_cov = 5

        self.paths = []
        self.crc_paths = []

        self.index_byte_num = 4
        self.path_len = 0
        # self.dataEncodingDnaLength = 150
        self.ratio_tolerance = 200

        self.set_max_kmer_num = False
        self.max_path_num = 1000
        self.max_num_of_repeat_kmer = 4
        self.path_sta = {}

        # Just for test only, will be removed later
        self.primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
        self.primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2


    # 20200807 modified to fix a bug with quanlity scores
    def count_file(self, file):
        self.counts.consume_seqfile(file)

    def set_kmer_len(self, km_len, size=1e8, table_num=4):
        self.kmer_len = km_len
        self.counts = khmer.Counttable(km_len, size, table_num)

    def find_droplets(self, index, droplet_template):
        path_len = 0
        if droplet_template.des:
            path_len = droplet_template.des_data_len + droplet_template.crc_len
        else:
            path_len = len(droplet_template.data) * 4 + droplet_template.crc_len

        self.find_paths(index, path_len)
        print('find_droplets path len:')
        print(path_len)


        all_droplets = []
        for adna in self.obtained_paths:
            adrop = DNADroplet(bytes(len(droplet_template.data)))
            adrop.data = bytes(len(droplet_template.data))
            # adrop.des = droplet_template.des
            # adrop.sec_key = droplet_template.sec_key
            adrop.num_of_chunks = droplet_template.num_of_chunks
            adrop.set_index(index)
            if adrop.set_droplet_from_DNA_CRC(adna[len(self.primerF):]):
                all_droplets.append(adrop)

        return all_droplets


    def find_paths(self, index, path_len=36*4):
        index_dna = bytesToDNA(index.to_bytes(self.index_byte_num, 'big', signed=False))
        print(index,end='\t')
        print(index_dna)
        self.obtained_paths = []
        self.paths = [self.primerF + index_dna]
        self.path_len = len(index_dna) + path_len

        i = 0
        while i < path_len and len(self.paths) < self.max_path_num:
            tmpPaths = []
            for path in self.paths:
                self.p_kmer = path[len(path) - self.kmer_len:len(path)]
                maxScore = self.score_next_bases()
                for block in self.next_bases.keys():
                    if self.next_bases[block] >= maxScore / self.ratio_tolerance and self.next_bases[block] > self.min_cov:
                        tmpPaths.append(path + block)

            if len(tmpPaths) > 0:
                self.paths = tmpPaths
            else:
                return 'No path found'
            i = i + 1
        if len(self.paths) > 0 and len(self.paths[0]) == self.path_len + len(self.primerF):
        #print(len(self.pathsA[0]))
            self.obtained_paths = self.paths
        self.path_sta[index] = len(self.paths)

    def score_next_bases(self):
        self.next_bases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        max_score = 0
        for m in self.next_bases.keys():
            self.next_bases[m] = self.score_next_base(m)
            if self.next_bases[m] > max_score:
                max_score = self.next_bases[m]
        return max_score

    def score_next_base(self, aBase):
        p_kmer = self.p_kmer
        next_kmer = p_kmer[1:len(p_kmer)] + aBase[0:1]
        return self.counts.get(next_kmer)











