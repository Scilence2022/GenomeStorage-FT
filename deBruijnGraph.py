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

    def find_droplets(self, index, droplet_template:DNADroplet):
        path_len = 0
        if droplet_template.des:
            path_len = droplet_template.des_data_len + droplet_template.crc_len
        else:
            path_len = len(droplet_template.data) * 4 + droplet_template.crc_len

        droplet_template.set_head_index(index)
        index_dna = droplet_template.get_head_index_dna()

        self.find_paths(index_dna, path_len)
        print(index, end='\t')
        # print('find_droplets path len:', end = '\t')
        # print(path_len)
        print('find_droplets path num:', end='\t')
        print(len(self.obtained_paths))

        all_droplets = []
        for adna in self.obtained_paths:
            adrop = DNADroplet(bytes(len(droplet_template.data)))
            adrop.data = bytes(len(droplet_template.data))
            # adrop.des = droplet_template.des
            # adrop.sec_key = droplet_template.sec_key
            adrop.num_of_chunks = droplet_template.num_of_chunks
            adrop.set_head_index(index)
            if adrop.set_droplet_from_DNA_CRC(adna[len(self.primerF):]):
                all_droplets.append(adrop)

        return all_droplets


    def extend_one_base(self, paths):
        """
        Extend all paths by one base pair.
        
        Args:
            paths: List of current paths to extend
            
        Returns:
            List of extended paths, or empty list if no valid extensions found
        """
        tmp_paths = []
        for path in self.paths:
            self.p_kmer = path[len(path) - self.kmer_len:len(path)]
            p_kmer_score = self.counts.get(self.p_kmer)
            max_score = self.score_next_bases()
            for block in self.next_bases.keys():
                if self.get_ratio(p_kmer_score, self.next_bases[block]) <=  self.ratio_tolerance and self.next_bases[block] > self.min_cov:
                    tmp_paths.append(path + block)

        if len(tmp_paths) == 0:
            return []
        else:
            return tmp_paths

    def get_ratio(self, a, b):
        if a == 0 or b == 0: return 0
        #
        return max(a, b) / min(a, b)

    def extend_paths(self, extend_len):
        """
        Extend existing paths by a specified length.
        
        Args:
            extend_len: Number of bases to extend the paths
            
        Returns:
            String message if no paths found, otherwise None
        """
        i = 0
        while i < extend_len and len(self.paths) < self.max_path_num:
            tmp_paths = self.extend_one_base(self.paths)
            
            if len(tmp_paths) > 0:
                self.paths = tmp_paths
            else:
                return 'No path found'
            i = i + 1
        
        return None

    def get_paths(self):
        return self.paths

    def init_paths(self, dna_seq):
        self.paths = [self.primerF + dna_seq]
        self.obtained_paths = []


    def find_paths(self, index_dna, path_len=36*4):
        """
        Initialize paths with primer and index DNA, then extend to target length.
        
        Args:
            index_dna: DNA sequence representing the index
            path_len: Target length for the paths (excluding primer)
        """
        print(index_dna)
        self.init_paths(index_dna)
        # self.obtained_paths = []
        # self.paths = [self.primerF + index_dna]
        self.path_len = len(index_dna) + path_len
        
        # Extend paths by the required length
        result = self.extend_paths(path_len)
        
        # Check if paths were successfully extended to target length
        if result is None and len(self.paths) > 0 and len(self.paths[0]) == self.path_len + len(self.primerF):
            self.obtained_paths = self.paths
        
        self.path_sta[index_dna] = len(self.paths)


    def score_next_bases(self):
        self.next_bases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        max_score = 0
        for m in self.next_bases.keys():
            self.next_bases[m] = self.score_next_base(m)
            if self.next_bases[m] > max_score:
                max_score = self.next_bases[m]
        # print(self.next_bases)
        return max_score

    def score_next_base(self, aBase):
        p_kmer = self.p_kmer
        next_kmer = p_kmer[1:len(p_kmer)] + aBase[0:1]
        return self.counts.get(next_kmer)











