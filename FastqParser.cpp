#include "FastqParser.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

// ============================================================
// 辅助排序函数（从大到小）
// ============================================================
void sortHigh(vector<pair<string, int> >& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = i + 1; j < vec.size(); ++j) {
            if (vec[j].second > vec[i].second) {
                swap(vec[i], vec[j]);
            }
        }
    }
}

// ============================================================
// 辅助排序函数（从小到大）
// ============================================================
void sortLow(vector<pair<string, int> >& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = i + 1; j < vec.size(); ++j) {
            if (vec[j].second < vec[i].second) {
                swap(vec[i], vec[j]);
            }
        }
    }
}

// ============================================================
// 1. 统计序列条数
// ============================================================
int countReads(const string& filename) {
    ifstream file(filename.c_str());
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return -1;
    }

    int readCount = 0;
    string line;
    while (getline(file, line)) {
        if (!line.empty() && line[0] == '@') {
            readCount++;
        }
    }
    file.close();
    return readCount;
}

// ============================================================
// 2. 提取所有序列
// ============================================================
vector<string> extractSequences(const string& filename) {
    vector<string> sequences;
    ifstream file(filename.c_str());
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return sequences;
    }

    string line;
    int lineCount = 0;

    while (getline(file, line)) {
        if (lineCount % 4 == 1) {
            sequences.push_back(line);
        }
        lineCount++;
    }
    file.close();
    return sequences;
}

// ============================================================
// 3. 计算 K-mer 频率
// ============================================================
unordered_map<string, int> computeKmerFrequencies(const vector<string>& sequences, int k) {
    unordered_map<string, int> kmerCount;

    for (size_t i = 0; i < sequences.size(); ++i) {
        const string& seq = sequences[i];
        if (seq.length() < (size_t)k) continue;

        for (size_t j = 0; j <= seq.length() - (size_t)k; ++j) {
            string kmer = seq.substr(j, k);
            unordered_map<string, int>::iterator it = kmerCount.find(kmer);
            if (it == kmerCount.end()) {
                kmerCount.insert(pair<string, int>(kmer, 1));
            } else {
                it->second++;
            }
        }
    }
    return kmerCount;
}

// ============================================================
// 4. 分析 k-mer 分布 + 异常检测 + CSV 导出
// ============================================================
void analyzeKmerDistribution(const unordered_map<string, int>& kmerFreq, int k, int topN, const string& filename) {
    if (kmerFreq.empty()) {
        cout << "No k-mers to analyze." << endl;
        return;
    }

    cout << "\n=== k-mer Distribution Summary (k=" << k << ") ===" << endl;

    size_t totalDistinct = kmerFreq.size();
    cout << "Total distinct " << k << "-mers: " << totalDistinct << endl;

    long long sum = 0;
    int maxCount = 0;
    string mostFreqKmer;

    for (unordered_map<string, int>::const_iterator it = kmerFreq.begin(); it != kmerFreq.end(); ++it) {
        int count = it->second;
        sum += count;
        if (count > maxCount) {
            maxCount = count;
            mostFreqKmer = it->first;
        }
    }

    double mean = static_cast<double>(sum) / totalDistinct;
    cout << "Average frequency per k-mer: " << mean << endl;

    double varianceSum = 0.0;
    for (unordered_map<string, int>::const_iterator it = kmerFreq.begin(); it != kmerFreq.end(); ++it) {
        double diff = it->second - mean;
        varianceSum += diff * diff;
    }
    double stddev = sqrt(varianceSum / totalDistinct);
    cout << "Standard deviation: " << stddev << endl;

    cout << "Most frequent " << k << "-mer: " << mostFreqKmer
         << " (appears " << maxCount << " times)" << endl;

    double thresholdHigh = mean + 2 * stddev;
    double thresholdLow = mean - 2 * stddev;

    cout << "\n=== Anomaly Detection (threshold: mean ± 2σ) ===" << endl;
    cout << "High frequency threshold: " << thresholdHigh << endl;
    cout << "Low frequency threshold: " << thresholdLow << endl;

    vector<pair<string, int> > highAnomalies;
    vector<pair<string, int> > lowAnomalies;

    for (unordered_map<string, int>::const_iterator it = kmerFreq.begin(); it != kmerFreq.end(); ++it) {
        if (it->second > thresholdHigh) {
            highAnomalies.push_back(*it);
        } else if (it->second < thresholdLow) {
            lowAnomalies.push_back(*it);
        }
    }

    sortHigh(highAnomalies);
    sortLow(lowAnomalies);

    cout << "\nTop " << topN << " most frequent k-mers (high anomalies):" << endl;
    if (highAnomalies.empty()) {
        cout << "  (No k-mers above high threshold)" << endl;
    } else {
        size_t limit = (highAnomalies.size() < (size_t)topN) ? highAnomalies.size() : (size_t)topN;
        for (size_t i = 0; i < limit; ++i) {
            cout << "  " << (i + 1) << ". " << highAnomalies[i].first
                 << ": " << highAnomalies[i].second << " times" << endl;
        }
    }

    cout << "\nTop " << topN << " least frequent k-mers (low anomalies):" << endl;
    if (lowAnomalies.empty()) {
        cout << "  (No k-mers below low threshold)" << endl;
    } else {
        size_t limit = (lowAnomalies.size() < (size_t)topN) ? lowAnomalies.size() : (size_t)topN;
        for (size_t i = 0; i < limit; ++i) {
            cout << "  " << (i + 1) << ". " << lowAnomalies[i].first
                 << ": " << lowAnomalies[i].second << " times" << endl;
        }
    }

    // ---- 导出 CSV（使用传入的 filename，它已经是合并后的名称） ----
    string csvFilename = filename + "_k" + to_string(k) + "_top" + to_string(topN) + ".csv";

    ofstream outFile(csvFilename.c_str());
    if (outFile.is_open()) {
        outFile << "# metadata\n";
        outFile << "# filename: " << filename << "\n";
        outFile << "# k: " << k << "\n";
        outFile << "# topN: " << topN << "\n";
        outFile << "kmer,frequency\n";
        size_t limit = (highAnomalies.size() < (size_t)topN) ? highAnomalies.size() : (size_t)topN;
        for (size_t i = 0; i < limit; ++i) {
            outFile << highAnomalies[i].first << "," << highAnomalies[i].second << "\n";
        }
        outFile.close();
        cout << "\n✅ Top " << topN << " data saved to " << csvFilename << endl;
    } else {
        cerr << "Warning: Could not write " << csvFilename << endl;
    }
}