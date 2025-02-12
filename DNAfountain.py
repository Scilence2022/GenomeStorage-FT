from DNAdroplet import DNADroplet
from fountain import Fountain
from utils import inverse_hash_32, reversible_hash_32, xor, randChunkNums, bytesToDNA, DNAToBytes
import random


class DNAFountain(Fountain):
	def __init__(self, data=b'', chunk_size=35, index=1, seed=1, droplet_head_length=16, droplet_tail_length=16, droplet_crc_length=8):
		self.droplet_head_length = droplet_head_length
		self.droplet_tail_length = droplet_tail_length
		self.droplet_crc_length = droplet_crc_length

		self.des = False
		self.sec_key = ""
		self.primer1 = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
		self.primer2 = 'CTGACACTGATGCATCCG'  # complement seq of P2
		self.note = "FT"
		self.note_len = 2
		super(DNAFountain, self).__init__(data, chunk_size, index, seed)


	def DNAdroplet(self):
		self.updateIndex()
		# print("DNAFountain index: " + str(self.index))
		aDroplet = DNADroplet(b'', self.index, self.num_of_chunks, self.droplet_head_length, self.droplet_tail_length, self.seed, self.droplet_crc_length)
		
		if self.des:
			aDroplet.des = True
			aDroplet.sec_key = self.sec_key

			# print("DNAFountain sec_key: " + self.sec_key)

		aDroplet.data_len = self.chunk_size * 4
		aDroplet.degree = self.degrees[self.index % len(self.degrees)]
		chunk_nums = aDroplet.get_chunk_nums()
		data = None
		for num in chunk_nums:
			if data is None:
				data = self.chunk(num)
			else:
				data = xor(data, self.chunk(num))
		aDroplet.data = data
		# print(aDroplet.)
	

		if self.des:
				aDroplet.encry_data()
		return aDroplet

	def fix_bytes(self):
		blen = len(self.data)
		l_bytes = blen % self.chunk_size
		if l_bytes > 0:
			self.data = self.data + bytes(self.chunk_size - l_bytes)

	def pack_metadata(self) -> bytes:
		"""Pack number_of_chunks into bytes with reversible hash"""
		# Convert number_of_chunks to 4 bytes (32-bit integer)
		hashed_chunks = reversible_hash_32(self.num_of_chunks)
		return hashed_chunks.to_bytes(4, byteorder='big')

	def unpack_metadata(self, metadata: bytes) -> None:
		"""Unpack bytes to restore number_of_chunks using inverse hash"""
		# Get hashed value from 4 bytes
		hashed_value = int.from_bytes(metadata[0:4], byteorder='big')
		# Recover original number using inverse hash
		self.num_of_chunks = inverse_hash_32(hashed_value)

	def primers(self) -> bytes:
		"""Return the primers as bytes"""
		self.update_primers()
		return [self.primer1, self.primer2]
	
	def update_primers(self) -> None:
		"""Set the primers"""
		self.primer1 = bytesToDNA(self.pack_metadata()) 
		self.primer2 = bytesToDNA(self.pack_metadata()) 