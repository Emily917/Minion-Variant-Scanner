# Minion-Variant-Scanner

**A C++ tool for FASTQ parsing, K-mer frequency analysis, and statistical anomaly detection in genomic sequencing data.**

[![C++](https://img.shields.io/badge/C++-11-blue.svg)](https://en.cppreference.com/w/cpp/11)
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)


## 🎯 What is this?

This project is a **lightweight, high-performance tool** for processing FASTQ files from third-generation sequencing (e.g., Oxford Nanopore MinION). It extracts DNA sequences, computes **K-mer frequencies**, and uses **statistical anomaly detection** (mean ± 2σ) to identify over-represented or under-represented sequence motifs.

**Key features:**

- Fast FASTQ parsing using C++ standard library
- Interactive command-line interface
- Multi-file merging support (analyze multiple FASTQ files together)
- K-mer frequency computation with `unordered_map`
- Statistical anomaly detection (mean ± 2 standard deviations)
- Automatic CSV export of Top N results
- Python-based visualization with automatic chart generation
- Charts saved to desktop `Minion_Results` folder with dynamic filenames

---

## 📊 Example Output
=== Minion Variant Scanner (Interactive Mode) ===

FASTQ file name（using space blank if there's more than 1 file, eg: SRR390728_1.fastq SRR390728_2.fastq）: SRR390728_1.fastq
Input K-mer length（eg: 11）: 11
Input Top N high/low frequency（eg: 10）: 10

Analyzing...
Number of files: 1

SRR390728_1.fastq
K-mer length: 11
Top N: 10

[1/1] read files: SRR390728_1.fastq
Sequence number: 1000
Successfully extracted: 1000 sequences
Cumulative total sequence: 1000

Calculate 11-mer frequencies...

=== k-mer Distribution Summary (k=11) ===
Total distinct 11-mers: 17057
Average frequency per k-mer: 1.52
Standard deviation: 1.09
Most frequent 11-mer: TTTTTTTTTTT (appears 44 times)

=== Anomaly Detection (threshold: mean ± 2σ) ===
High frequency threshold: 3.70
Low frequency threshold: -0.65

Top 10 most frequent k-mers (high anomalies):

TTTTTTTTTTT: 44 times

GGGGGGGGGGG: 15 times

AAAAAAAAAAA: 11 times

TCTCTAGAACG: 10 times

ATCTCTAGAAC: 10 times

GCATCCTTTTT: 10 times

CTCTAGAACGT: 10 times

AGAGCGCATCC: 9 times

AGAATCCAACT: 9 times

GGAGAATCCAA: 9 times

Top 10 least frequent k-mers (low anomalies):
(No k-mers below low threshold)

✅ Top 10 data saved to SRR390728_1_k11_top10.csv

✅ All analyses complete！
Merge file: SRR390728_1
CSV file: SRR390728_1_k11_top10.csv
input --- python plot_kmer.py --- generate the chart.

<img width="983" height="491" alt="Screenshot 2026-07-13 at 02 29 32" src="https://github.com/user-attachments/assets/53337481-8035-4de3-b449-3c6b8c2cf1f1" />

## 🛠️ Build & Run

**Prerequisites**

- C++11 compatible compiler (g++ or clang)
- Python 3.8+ with `matplotlib` installed
- macOS / Linux / WSL

**Install Python Dependencies**
g++ *.cpp -o scanner

## 📁 Project Structure
<img width="689" height="177" alt="Screenshot 2026-07-13 at 02 10 53" src="https://github.com/user-attachments/assets/8fbe0d90-f6b8-4fd8-ac9c-52498bed03a5" />

## 🧬 Biological Context
This tool is inspired by state-of-the-art variant callers like **Clair** (developed by Prof. Ruibang Luo's group at HKU). It addresses the high error rate of **third-generation sequencing (long-read)** data by:

1. Converting raw reads into **K-mer features**
2. Using **statistical anomaly detection** to flag over-represented sequences
3. These flagged regions often correspond to **repetitive elements**, **regulatory motifs**, or **sequencing artifacts**

Future work will integrate machine learning models (SVM / PyTorch) for SNP and indel detection.

## 📚 References

- [Clair3](https://github.com/HKU-BAL/Clair3) – State-of-the-art variant caller by Prof. Ruibang Luo
- [SRA Toolkit](https://github.com/ncbi/sra-tools) – For downloading FASTQ data

## 📬 Contact

**Emily917** – [13372579039@163.com](mailto:13372579039@163.com)

GitHub: [github.com/Emily917/Minion-Variant-Scanner](https://github.com/Emily917/Minion-Variant-Scanner)

## 📄 License

This project is licensed under the MIT License.
