#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdlib>  // for system()
#include "FastqParser.hpp"

using namespace std;

/*
cd /Users/e./Desktop/Minion-Variant-Scanner
 source .venv/bin/activate 
./scanner 
（运行指令）

source .venv/bin/activate 
python plot_kmer.py 
（绘图指令）
*/


int main() {
    cout << "=== Minion Variant Scanner v2.0 ===" << endl;
    cout << "Support: 3-line & 4-line FASTQ format" << endl;
    cout << endl;

    while (true) {  // 支持多次分析
        // 交互式输入
        string inputLine;
        cout << "FASTQ file name (space separated, or 'quit' to exit): ";
        getline(cin, inputLine);

        if (inputLine == "quit" || inputLine == "q" || inputLine == "exit") {
            cout << "Goodbye!" << endl;
            break;
        }

        // 分割文件名
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
            cout << "Error: No filename provided" << endl;
            continue;
        }

        int k;
        cout << "K-mer length (default 11): ";
        string kInput;
        getline(cin, kInput);
        if (!kInput.empty()) k = stoi(kInput);
        else k = 11;

        int topN;
        cout << "Top N (default 10): ";
        string nInput;
        getline(cin, nInput);
        if (!nInput.empty()) topN = stoi(nInput);
        else topN = 10;

        cout << endl;
        cout << "Analyzing..." << endl;
        cout << "  Files: " << filenames.size() << endl;
        for (size_t i = 0; i < filenames.size(); ++i) {
            cout << "    " << (i+1) << ". " << filenames[i] << endl;
        }
        cout << "  K-mer length: " << k << endl;
        cout << "  Top N: " << topN << endl;
        cout << endl;

        auto startTime = chrono::high_resolution_clock::now();

        // 创建解析器
        FastqParser parser;

        // 读取所有序列
        vector<string> allSequences;
        allSequences.reserve(1000000);

        int totalReads = 0;
        for (size_t i = 0; i < filenames.size(); ++i) {
            cout << "[" << (i+1) << "/" << filenames.size() << "] Reading: " << filenames[i] << endl;

            int readCount = parser.countReads(filenames[i]);
            if (readCount < 0) {
                cout << "  ⚠️  Skipping file" << endl;
                continue;
            }
            totalReads += readCount;
            cout << "  Total reads: " << readCount << endl;

            vector<string> sequences = parser.extractSequences(filenames[i], 0);
            allSequences.insert(allSequences.end(), sequences.begin(), sequences.end());

            cout << "  Cumulative: " << allSequences.size() << " sequences" << endl;
            cout << endl;
        }

        if (allSequences.empty()) {
            cout << "Error: No valid sequences read" << endl;
            continue;
        }

        cout << "========================================" << endl;
        cout << "Total sequences: " << allSequences.size() << endl;
        cout << "========================================" << endl;

        // 计算K-mer频率
        cout << "\nCalculating " << k << "-mer frequencies..." << endl;
        auto kmerStart = chrono::high_resolution_clock::now();

        unordered_map<uint64_t, int> kmerFreq = parser.computeKmerFrequencies(allSequences, k);

        auto kmerEnd = chrono::high_resolution_clock::now();
        chrono::duration<double> kmerTime = kmerEnd - kmerStart;
        cout << "  K-mer computation time: " << kmerTime.count() << " seconds" << endl;
        cout << "  Distinct k-mers: " << kmerFreq.size() << endl;

        // 生成输出文件名
        string mergedName;
        for (size_t i = 0; i < filenames.size(); ++i) {
            string baseName = filenames[i];
            size_t dotPos = baseName.find_last_of('.');
            if (dotPos != string::npos) {
                baseName = baseName.substr(0, dotPos);
            }
            if (i > 0) mergedName += "_and_";
            mergedName += baseName;
        }

        // 分析并导出
        parser.analyzeAndExport(kmerFreq, k, topN, mergedName);

        auto endTime = chrono::high_resolution_clock::now();
        chrono::duration<double> totalTime = endTime - startTime;
        cout << "\n⏱️  Total execution time: " << totalTime.count() << " seconds" << endl;

        // ===== 关键改进：自动运行 Python 绘图 =====
        cout << "\n📊 Generating plot from CSV..." << endl;
        int result = system("python plot_kmer.py");
        if (result == 0) {
            cout << "✅ Plot generated successfully!" << endl;
            cout << "   Check: ~/Desktop/Minion_Results/" << endl;
        } else {
            cout << "⚠️  Could not run plot_kmer.py" << endl;
            cout << "   You can manually run: python plot_kmer.py" << endl;
        }

        cout << "\n========================================" << endl;
        cout << "✅ Analysis and plot complete!" << endl;
        cout << "========================================" << endl;
        cout << "\nPress Enter to analyze another file, or type 'quit' at the prompt..." << endl;
    }

    return 0;
}