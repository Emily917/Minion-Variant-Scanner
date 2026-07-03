# Minion-Variant-Scanner

**A C++ tool for FASTQ parsing and K-mer frequency computation in genomic variant detection.**

[![C++](https://img.shields.io/badge/C++-11-blue.svg)](https://en.cppreference.com/w/cpp/11)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## 🎯 What is this?

This project is a lightweight, high-performance **FASTQ file parser and K-mer counter** written in C++. It is designed to handle DNA sequencing data from third-generation sequencers (e.g., Oxford Nanopore MinION) and extract meaningful features for variant detection.

**Key features:**
- Fast FASTQ parsing using C++ standard library
- DNA sequence extraction from raw sequencing data
- K-mer frequency computation with `unordered_map`
- Modular design for easy extension

## 📊 Sample Output

```bash
=== Minion Variant Scanner ===
Step 1: Counting reads in FASTQ file...
Total reads in SRR390728_1.fastq: 1000

Step 2: Extracting sequences...
Extracted 1000 sequences.

Step 3: Computing 11-mer frequencies...
Total distinct 11-mers: 17057
Most frequent 11-mer: TTTTTTTTTTT (appears 44 times)
```

## 🛠️ Build & Run

### Prerequisites
- C++11 compatible compiler (g++ or clang)
- macOS / Linux / WSL

### Compile
```bash
g++ -std=c++11 *.cpp -o scanner
```

### Run
```bash
./scanner
```

### With custom FASTQ file
```bash
# Modify the filename in main.cpp, or use command-line args (coming soon)
```

## 📁 Project Structure

```
Minion-Variant-Scanner/
├── main.cpp              # Entry point
├── FastqParser.hpp       # Header file
├── FastqParser.cpp       # Implementation
├── SRR390728_1.fastq     # Test data (forward reads)
├── SRR390728_2.fastq     # Test data (reverse reads)
└── README.md             # This file
```

## 🧬 Biological Context

This tool is a **"starter version"** of variant callers like **Clair** (developed by Prof. Ruibang Luo's group at HKU). It addresses the high error rate of **third-generation sequencing (long-read)** data by converting raw reads into K-mer features, which can later be fed into machine learning models for SNP and indel detection.

## 🚀 Next Steps (Roadmap)

- [ ] Command-line argument support (input file, K-mer length)
- [ ] Performance optimization with integer-encoded K-mers
- [ ] Multi-threaded parsing
- [ ] Integration with lightweight ML models (SVM / PyTorch)

## 📚 References

- [Clair3](https://github.com/HKU-BAL/Clair3) – Prof. Ruibang Luo's state-of-the-art variant caller
- [SRA Toolkit](https://github.com/ncbi/sra-tools) – For downloading FASTQ data

## 📬 Contact

Emily917 – [13372579039@163.com](mailto:13372579039@163.com)

GitHub: [github.com/Emily917/Minion-Variant-Scanner](https://github.com/Emily917/Minion-Variant-Scanner)
