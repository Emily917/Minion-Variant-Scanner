#include <iostream>
#include "FastqParser.hpp"

using namespace std;

//g++ *.cpp -o scanner && ./scanner 开始输入终端
//python plot_kmer.py 生成图像

int main() {
    cout << "=== Minion Variant Scanner (Interactive Mode) ===" << endl;
    cout << endl;

    // ---- 交互式输入 ----
    string inputLine;
    cout << "FASTQ file name（using space blank if there's more than 1 file, eg: SRR390728_1.fastq SRR390728_2.fastq）: ";
    getline(cin, inputLine);

    // 用空格分割文件名
    vector<string> filenames;
    string currentFile;
    for (size_t i = 0; i <= inputLine.length(); ++i) {
        if (i == inputLine.length() || inputLine[i] == ' ') {
            if (!currentFile.empty()) {
                filenames.push_back(currentFile);
                currentFile.clear();
            }
        } else {
            currentFile += inputLine[i];
        }
    }

    if (filenames.empty()) {
        cout << "error: didn't recognize File name" << endl;
        return 1;
    }

    int k;
    cout << "Input K-mer length（eg: 11）: ";
    cin >> k;

    int topN;
    cout << "Input Top N high/low frequency（eg: 10）: ";
    cin >> topN;

    cout << endl;
    cout << "Analyzing..." << endl;
    cout << "  Number of files: " << filenames.size() << endl;
    for (size_t i = 0; i < filenames.size(); ++i) {
        cout << "    " << (i + 1) << ". " << filenames[i] << endl;
    }
    cout << "  K-mer length: " << k << endl;
    cout << "  Top N: " << topN << endl;
    cout << endl;

    // ---- 1. 统计并提取所有文件的序列 ----
    vector<string> allSequences;
    int totalReads = 0;

    for (size_t i = 0; i < filenames.size(); ++i) {
        string filename = filenames[i];
        cout << "[" << (i + 1) << "/" << filenames.size() << "] read files: " << filename << endl;

        int readCount = countReads(filename);
        if (readCount < 0) {
            cout << "  Jump（can't open the file）" << endl;
            continue;
        }
        totalReads += readCount;
        cout << "  Sequence number: " << readCount << endl;

        vector<string> sequences = extractSequences(filename);
        cout << "  Successfully extracted: " << sequences.size() << " sequences" << endl;

        // 合并到总容器
        allSequences.insert(allSequences.end(), sequences.begin(), sequences.end());
        cout << "  Cumulative total sequence: " << allSequences.size() << endl;
        cout << endl;
    }

    if (allSequences.empty()) {
        cout << "Error: Failed to read any sequences." << endl;
        return 1;
    }

    cout << "========================================" << endl;
    cout << "Total after consolidation: " << allSequences.size() << " sequences" << endl;
    cout << "========================================" << endl;

    // ---- 2. 计算 K-mer 频率 ----
    cout << "Calculate " << k << "-mer frequencies..." << endl;
    unordered_map<string, int> kmerFreq = computeKmerFrequencies(allSequences, k);

    // ---- 3. 生成合并后的文件名（用于 CSV 和 PNG 命名） ----
    string mergedName;
    for (size_t i = 0; i < filenames.size(); ++i) {
        // 去掉 .fastq 或 .fq 后缀
        string baseName = filenames[i];
        size_t dotPos = baseName.find_last_of('.');
        if (dotPos != string::npos) {
            baseName = baseName.substr(0, dotPos);
        }
        if (i > 0) {
            mergedName += "_and_";
        }
        mergedName += baseName;
    }

    // ---- 4. 分析并导出 CSV（传入合并后的文件名） ----
    analyzeKmerDistribution(kmerFreq, k, topN, mergedName);

    cout << "\n✅ All analyses complete！" << endl;
    cout << "   Merge file: " << mergedName << endl;
    cout << "   CSV file: " << mergedName << "_k" << k << "_top" << topN << ".csv" << endl;
    cout << "   input --- python plot_kmer.py --- generate the chart." << endl;

    return 0;
}