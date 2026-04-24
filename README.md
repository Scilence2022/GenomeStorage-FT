# GenomeStorage-FT

A robust DNA-based data storage system using Luby Transform (LT) fountain codes and de Bruijn graph assembly. GenomeStorage-FT encodes arbitrary digital files into synthetic DNA sequences and recovers them from sequencing data, even in the presence of significant dropout and sequencing errors.

## Overview

GenomeStorage-FT translates binary data into DNA sequences using fountain codes, which provide inherent redundancy and rateless properties. The system is designed to tolerate strand dropouts and sequencing errors through:

- **Fountain Code Encoding**: Robust soliton distribution-based LT codes for near-optimal redundancy.
- **DNA Sequence Constraints**: GC-content balancing and homopolymer filtering to ensure synthesizability.
- **CRC Error Detection**: Per-droplet CRC-16/XMODEM checksums for integrity validation.
- **De Bruijn Graph Decoding**: k-mer-based assembly to reconstruct droplets from noisy sequencing reads.
- **Optional DES Encryption**: Data confidentiality via DES-CBC encryption before DNA encoding.

## Features

- **Rateless Encoding**: Generate as many DNA droplets as needed; only a subset is required for successful decoding.
- **Robust to Dropouts**: Handles strand loss during synthesis, storage, or sequencing.
- **Sequencing-Aware Decoding**: Uses de Bruijn graph traversal to recover original droplets from raw reads.
- **Flexible Parameters**: Configurable chunk sizes, redundancy rates, primer sequences, and index lengths.
- **Multiple Output Formats**: FASTA, tab-delimited, and detailed log files.
- **Strand Dropout Simulation**: Built-in testing mode to evaluate robustness under simulated loss conditions.

## Installation

### Prerequisites

- Python 3.6+
- pip

### Dependencies

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Core dependencies include:

| Package | Version | Purpose |
|---------|---------|---------|
| `numpy` | >=1.19.0 | Numerical computations, distribution calculations |
| `pyDes` | >=2.0.1 | DES encryption for optional data security |
| `khmer` | >=2.1.1 | Efficient k-mer counting for de Bruijn graph construction |

## Project Structure

```
GenomeStorage-FT/
├── encode.py              # Main encoding script (file -> DNA sequences)
├── decode.py              # Main decoding script (sequencing reads -> file)
├── decode-DBGPS.py        # Decoder for pre-assembled DBGPS strand data
├── DNAfountain.py         # DNAFountain class: fountain code encoder with DNA constraints
├── DNAdroplet.py          # DNADroplet class: individual DNA-encoded droplet
├── fountain.py            # Fountain base class: LT code generation and chunking
├── droplet.py             # Droplet base class: symbol representation
├── deBruijnGraph.py       # DeBruijnGraph class: k-mer graph for sequence assembly
├── glass.py               # Glass class: belief-propagation decoder for fountain codes
├── seqFountain.py         # SeqFountain class: sequencing read handling
├── utils.py               # Utility functions (DNA/byte conversion, distributions, hashing)
├── test_utils.py          # Testing utilities and simulation helpers
├── crc16pure.py           # Pure Python CRC-16/XMODEM implementation
├── requirements.txt       # Python dependency list
├── input_files/           # Example input data files
└── output_files/          # Example output data files
```

## Usage

### Encoding

Convert a binary file into DNA sequences:

```bash
python encode.py -i <input_file> -n <number_of_droplets> -o <output.fasta> [options]
```

#### Encoding Options

| Option | Long Form | Description | Default |
|--------|-----------|-------------|---------|
| `-h` | `--help` | Show help message | - |
| `-i` | `--input` | Input file path | **Required** |
| `-o` | `--output` | Output FASTA file path | `<input>.fasta` |
| `-n` | `--droplet_num` | Number of droplets to generate | 3650 |
| `-r` | `--redundancy_rate` | Redundancy rate (auto-calculates droplet count) | 0 |
| `-c` | `--chunk_size` | Chunk size in bytes | 35 |
| `-s` | `--seed` | Fountain random seed | 1 |
| `-l` | `--initial_index` | Initial droplet index | 1 |
| | `--index_bytes` | Index and anchor code length in bytes | 4 |
| | `--ec_bytes` | Error correction (CRC) length in bytes | 2 |
| | `--test_num` | Number of dropout simulation runs | 0 |
| | `--test_dropout_rate` | Simulated strand dropout rate | 0.05 |
| | `--p1` | 5'->3' primer sequence | `CCTGCAGAGTAGCATGTC` |
| | `--p2` | Complement primer sequence | `CTGACACTGATGCATCCG` |

#### Encoding Example

```bash
python encode.py -i input_files/R.MG1655.fasta -n 5000 -o output_files/MG1655_encoded.fasta -c 35 -s 1
```

This will produce:
- `output_files/MG1655_encoded.fasta` — DNA sequences in FASTA format
- `output_files/MG1655_encoded.fasta.log` — Encoding details (index, data, DNA sequence, degree, chunk list, tail index)
- `output_files/MG1655_encoded.fasta.tab` — Tab-delimited sequence list

### Decoding

Recover the original file from sequencing reads or assembled DNA sequences:

```bash
python decode.py -i <input_file> -o <output_file> [options]
```

#### Decoding Options

| Option | Long Form | Description | Default |
|--------|-----------|-------------|---------|
| `-h` | `--help` | Show help message | - |
| `-i` | `--input` | Input sequencing file (FastQ, FastA, or Jellyfish dump) | **Required** |
| `-o` | `--output` | Output file path | `output_files/<cutoff>.decoded.ECC` |
| `-t` | `--file_type` | Input file type: FastQ, Fasta, or Jellyfish dumped k-mers | k-mer dump |
| `-k` | `--kmer_size` | k-mer size for de Bruijn graph | 21 |
| | `--cut` | Coverage cutoff for low-frequency k-mer elimination | 10 |
| `-c` | `--chunk_size` | Chunk size in bytes (must match encoding) | 35 |
| `-n` | `--chunk_num` | Number of chunks (must match encoding) | 3142 |
| `-s` | `--seed` | Fountain random seed (must match encoding) | 1 |
| | `--min_index` | Minimum index to search | 1 |
| | `--max_index` | Maximum index to search | 3700 |
| | `--index_bytes` | Index code length in bytes | 4 |
| | `--ec_bytes` | CRC length in bytes | 2 |
| | `--total_bytes` | Total bytes of original file | 109953 |
| | `--primer` | Forward primer sequence | `TATCGATGACCTCGAGGA` |

#### Decoding Example

```bash
python decode.py -i input_files/Ecoli-MG1655-LR881938.1.ECC.fasta -o output_files/decoded.bin -c 35 -n 3142 -s 1 --total_bytes 109953
```

### Decoding Pre-Assembled Strands (DBGPS)

If you have strands already assembled by DBGPS (de Bruijn Graph Path Search):

```bash
python decode-DBGPS.py -i <assembled_strands_file> -o <output_file> -p <password> [options]
```

## Core Concepts

### Fountain Codes

GenomeStorage-FT uses Luby Transform (LT) codes with a robust soliton distribution. The source file is split into fixed-size chunks. Each droplet is formed by XORing a randomly selected subset of chunks. The random subset (degree and chunk indices) is deterministically derived from the droplet's index and a shared seed, allowing the decoder to reconstruct the same mapping without additional metadata.

### DNA Encoding Scheme

Each droplet is serialized as:

```
[Head Index Hash] + [Data Payload] + [CRC Checksum]
```

- **Head Index Hash**: Jenkins hash of the droplet index (used as a pseudorandom seed)
- **Data Payload**: XOR-combined chunk data (optionally DES-encrypted)
- **CRC Checksum**: CRC-16/XMODEM for error detection

This byte sequence is converted to DNA using a fixed 2-base-pair encoding per byte, producing sequences that are flanked by PCR primer sequences for amplification.

### DNA Constraints

Droplets are validated against synthesizability constraints:
- **GC Content**: 40% to 60%
- **Homopolymer Run Length**: Maximum 7 consecutive identical bases

Droplets failing these constraints are discarded and regenerated.

### De Bruijn Graph Assembly

During decoding, sequencing reads are decomposed into k-mers. The decoder initiates path searches from known index DNA sequences through the de Bruijn graph, extending paths base-by-base using k-mer coverage information. Successfully reconstructed droplet sequences are validated by CRC before being fed into the fountain code decoder.

### Belief Propagation Decoding

The `Glass` class implements the standard LT code peeling decoder:
1. Identify droplets of degree 1 (containing a single unknown chunk).
2. Recover that chunk directly.
3. XOR the recovered chunk out of all other droplets containing it.
4. Repeat until all chunks are recovered or no more degree-1 droplets exist.

## API Reference

### DNAFountain

```python
from DNAfountain import DNAFountain

fdna = DNAFountain(data, chunk_size=35, index=1, seed=1,
                   droplet_head_length=16, droplet_tail_length=16,
                   droplet_crc_length=8)
fdna.gen_degrees()
droplet = fdna.DNAdroplet()  # Generate next droplet
```

### DeBruijnGraph

```python
from deBruijnGraph import DeBruijnGraph

deG = DeBruijnGraph(kmer_len=21)
deG.set_kmer_len(kmer_len)
deG.count_file(sequencing_file)
droplets = deG.find_droplets(index, droplet_template)
```

### Glass

```python
from glass import Glass

cup = Glass(num_chunks)
cup.addDroplet(droplet)
cup.decode()
if cup.isDone():
    cup.writeToFile(output_path, total_bytes)
```

## Input/Output File Formats

### FASTA Output (`*.fasta`)

```
>1
CCTGCAGAGTAGCATGTC<encoded DNA>CTGACACTGATGCATCCG
>2
CCTGCAGAGTAGCATGTC<encoded DNA>CTGACACTGATGCATCCG
...
```

### Log Output (`*.log`)

Tab-delimited file containing:
- Head index
- Raw data bytes
- DNA sequence
- Degree
- Chunk number list
- Tail index

### Tab Output (`*.tab`)

```
head_index	sequence
1	CCTGCAGAGTAGCATGTC...CTGACACTGATGCATCCG
2	CCTGCAGAGTAGCATGTC...CTGACACTGATGCATCCG
```

## Example Dataset

The `input_files/` directory contains sample datasets:

| File | Description |
|------|-------------|
| `R.MG1655.fasta` | Reference genome sequence for E. coli MG1655 |
| `Ecoli-MG1655-LR881938.1.ECC` | Encoded data for E. coli MG1655 |
| `Corg-13032-Fixed-nxm.ECC` | Encoded data for C. glutamicum |

## Testing

Run built-in dropout simulation during encoding:

```bash
python encode.py -i myfile.bin -n 10000 --test_num 100 --test_dropout_rate 0.1
```

This simulates 100 decoding attempts with a 10% strand dropout rate and reports the success rate.

## License

This project is licensed under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.

## Citation

If you use GenomeStorage-FT in your research, please cite:

> GenomeStorage-FT: A Fountain Code-Based DNA Data Storage System with de Bruijn Graph Assembly.

## Acknowledgments

- LT fountain codes based on Luby's original design.
- k-mer counting powered by the [khmer](https://github.com/dib-lab/khmer) library.
- DES encryption provided by [pyDes](https://github.com/twhiteman/pyDes).