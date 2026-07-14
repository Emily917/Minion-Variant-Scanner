#include "FastqParser.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <future>
#include <cctype>

using namespace std;

// ============ 线程池 ============
class FastqParser::ThreadPool {
public:
    ThreadPool(size_t numThreads = thread::hardware_concurrency()) 
        : stop(false) {
        for (size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    function<void()> task;
                    {
                        unique_lock<mutex> lock(queueMutex);
                        condition.wait(lock, [this] { 
                            return stop || !tasks.empty(); 
                        });
                        if (stop && tasks.empty()) return;
                        task = move(tasks.front());
                        tasks.pop();
                    }
                    task();
                }
            });
        }
    }
    
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> future<typename result_of<F(Args...)>::type> {
        using returnType = typename result_of<F(Args...)>::type;
        auto task = make_shared<packaged_task<returnType()>>(
            bind(forward<F>(f), forward<Args>(args)...)
        );
        future<returnType> result = task->get_future();
        {
            unique_lock<mutex> lock(queueMutex);
            if (stop) throw runtime_error("enqueue on stopped ThreadPool");
            tasks.emplace([task]() { (*task)(); });
        }
        condition.notify_one();
        return result;
    }
    
    ~ThreadPool() {
        {
            unique_lock<mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for (thread& worker : workers) {
            worker.join();
        }
    }
    
private:
    vector<thread> workers;
    queue<function<void()>> tasks;
    mutex queueMutex;
    condition_variable condition;
    bool stop;
};

// ============ FastqParser 实现 ============

FastqParser::FastqParser() : pool_(new ThreadPool()) {}
FastqParser::~FastqParser() { delete pool_; }

// 检测FASTQ格式：3行（无质量）还是4行（有质量）
bool FastqParser::detectFormat(const string& filename) {
    ifstream file(filename.c_str());
    if (!file.is_open()) return false;
    
    vector<string> lines;
    string line;
    while (getline(file, line) && lines.size() < 10) {
        if (!line.empty()) lines.push_back(line);
    }
    file.close();
    
    if (lines.size() < 4) return false;  // 默认4行格式
    
    // 如果第4行（索引3）以@开头，说明是3行格式（没有质量行）
    // 因为3行格式: @行, 序列, +行, @行...
    if (lines[3][0] == '@') {
        return true;  // 3行格式
    }
    
    return false;  // 4行格式
}

// 检查是否为有效DNA序列
bool FastqParser::isValidDNA(const string& seq) {
    for (char c : seq) {
        char uc = toupper(c);
        if (uc != 'A' && uc != 'C' && uc != 'G' && uc != 'T' && uc != 'N') {
            return false;
        }
    }
    return !seq.empty();
}

// 编码K-mer (A=0, C=1, G=2, T=3)
uint64_t FastqParser::encodeKmer(const string& kmer) {
    uint64_t encoded = 0;
    for (char c : kmer) {
        encoded <<= 2;
        switch (toupper(c)) {
            case 'A': encoded |= 0; break;
            case 'C': encoded |= 1; break;
            case 'G': encoded |= 2; break;
            case 'T': encoded |= 3; break;
            default:  encoded |= 0; break;  // N 当作 A
        }
    }
    return encoded;
}

// 解码K-mer
string FastqParser::decodeKmer(uint64_t encoded, int k) {
    string kmer;
    kmer.resize(k);
    for (int i = k - 1; i >= 0; --i) {
        int code = encoded & 3;
        switch (code) {
            case 0: kmer[i] = 'A'; break;
            case 1: kmer[i] = 'C'; break;
            case 2: kmer[i] = 'G'; break;
            case 3: kmer[i] = 'T'; break;
        }
        encoded >>= 2;
    }
    return kmer;
}

// 1. 统计read数量（自动检测格式）
int FastqParser::countReads(const string& filename) {
    ifstream file(filename.c_str());
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return -1;
    }
    
    bool isThreeLine = detectFormat(filename);
    if (isThreeLine) {
        cout << "  Detected 3-line FASTQ format (no quality scores)" << endl;
    } else {
        cout << "  Detected 4-line FASTQ format" << endl;
    }
    
    // 重新打开文件
    file.clear();
    file.seekg(0, ios::beg);
    
    int readCount = 0;
    string line;
    int lineNum = 0;
    
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        if (isThreeLine) {
            // 3行格式：每3行一个read，@行在每组的第1行
            if (lineNum % 3 == 0 && line[0] == '@') {
                readCount++;
            }
        } else {
            // 4行格式：每4行一个read，@行在每组的第1行
            if (lineNum % 4 == 0 && line[0] == '@') {
                readCount++;
            }
        }
        lineNum++;
    }
    file.close();
    return readCount;
}

// 2. 提取序列（支持3行格式）
vector<string> FastqParser::extractSequences(const string& filename, size_t maxReads) {
    vector<string> sequences;
    if (maxReads == 0) maxReads = SIZE_MAX;
    
    ifstream file(filename.c_str());
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return sequences;
    }
    
    bool isThreeLine = detectFormat(filename);
    
    // 重新打开文件
    file.clear();
    file.seekg(0, ios::beg);
    
    string line;
    int lineNum = 0;
    sequences.reserve(min(maxReads, size_t(1000000)));
    
    while (getline(file, line) && sequences.size() < maxReads) {
        if (line.empty()) continue;
        
        if (isThreeLine) {
            // 3行格式：@行(0), 序列(1), +行(2)
            if (lineNum % 3 == 1) {  // 序列行
                if (isValidDNA(line)) {
                    sequences.push_back(line);
                } else {
                    // 尝试清理序列
                    string cleanSeq;
                    cleanSeq.reserve(line.size());
                    for (char c : line) {
                        char uc = toupper(c);
                        if (uc == 'A' || uc == 'C' || uc == 'G' || uc == 'T' || uc == 'N') {
                            cleanSeq += uc;
                        }
                    }
                    if (!cleanSeq.empty()) {
                        sequences.push_back(cleanSeq);
                    }
                }
            }
        } else {
            // 4行格式：@行(0), 序列(1), +行(2), 质量(3)
            if (lineNum % 4 == 1) {  // 序列行
                if (isValidDNA(line)) {
                    sequences.push_back(line);
                }
            }
        }
        lineNum++;
    }
    file.close();
    
    cout << "  Extracted " << sequences.size() << " valid sequences" << endl;
    return sequences;
}

// 3. 并行计算K-mer频率（使用整数编码 + 滑动窗口）
unordered_map<uint64_t, int> FastqParser::computeKmerFrequencies(
    const vector<string>& sequences, int k) {
    
    unordered_map<uint64_t, int> kmerCount;
    if (sequences.empty() || k <= 0) return kmerCount;
    
    kmerCount.reserve(sequences.size() * 2);
    
    // 使用线程池并行处理
    mutex mtx;
    vector<future<void>> futures;
    size_t numThreads = min(
    static_cast<size_t>(thread::hardware_concurrency()), 
    sequences.size()
);
   
    if (numThreads == 0) numThreads = 1;
    size_t chunkSize = max(size_t(1), sequences.size() / numThreads);
    
    cout << "  Using " << numThreads << " threads..." << endl;
    
    for (size_t start = 0; start < sequences.size(); start += chunkSize) {
        size_t end = min(start + chunkSize, sequences.size());
        futures.push_back(pool_->enqueue([&, start, end, k] {
            unordered_map<uint64_t, int> localMap;
            localMap.reserve((end - start) * 2);
            
            for (size_t i = start; i < end; ++i) {
                const string& seq = sequences[i];
                if (seq.length() < static_cast<size_t>(k)) continue;
                
                // 初始编码前k个碱基
                uint64_t kmer = 0;
                for (int j = 0; j < k; ++j) {
                    kmer <<= 2;
                    char c = toupper(seq[j]);
                    switch(c) {
                        case 'A': kmer |= 0; break;
                        case 'C': kmer |= 1; break;
                        case 'G': kmer |= 2; break;
                        case 'T': kmer |= 3; break;
                        default: kmer |= 0; break;
                    }
                }
                
                // 滑动窗口
                uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2*k)) - 1);
                for (size_t j = 0; j <= seq.length() - static_cast<size_t>(k); ++j) {
                    if (j > 0) {
                        // 左移2位，加入新碱基
                        kmer <<= 2;
                        char c = toupper(seq[j + k - 1]);
                        switch(c) {
                            case 'A': kmer |= 0; break;
                            case 'C': kmer |= 1; break;
                            case 'G': kmer |= 2; break;
                            case 'T': kmer |= 3; break;
                            default: kmer |= 0; break;
                        }
                        kmer &= mask;
                    }
                    
                    auto it = localMap.find(kmer);
                    if (it == localMap.end()) {
                        localMap.emplace(kmer, 1);
                    } else {
                        it->second++;
                    }
                }
            }
            
            // 合并到全局map
            lock_guard<mutex> lock(mtx);
            for (const auto& pair : localMap) {
                auto it = kmerCount.find(pair.first);
                if (it == kmerCount.end()) {
                    kmerCount.emplace(pair.first, pair.second);
                } else {
                    it->second += pair.second;
                }
            }
        }));
    }
    
    // 等待所有线程完成
    for (auto& f : futures) {
        f.get();
    }
    
    return kmerCount;
}

// 4. 获取统计信息
KmerStatistics FastqParser::getStatistics(const unordered_map<uint64_t, int>& kmerFreq, int k) {
    KmerStatistics stats;
    stats.totalDistinct = kmerFreq.size();
    stats.totalOccurrences = 0;
    stats.maxCount = 0;
    stats.mostFreqKmer = "";
    
    if (kmerFreq.empty()) {
        stats.mean = 0;
        stats.stddev = 0;
        stats.highThreshold = 0;
        stats.lowThreshold = 0;
        return stats;
    }
    
    for (const auto& pair : kmerFreq) {
        int count = pair.second;
        stats.totalOccurrences += count;
        if (count > stats.maxCount) {
            stats.maxCount = count;
            stats.mostFreqKmer = decodeKmer(pair.first, k);
        }
    }
    
    stats.mean = static_cast<double>(stats.totalOccurrences) / stats.totalDistinct;
    
    // 计算标准差
    double varianceSum = 0.0;
    for (const auto& pair : kmerFreq) {
        double diff = pair.second - stats.mean;
        varianceSum += diff * diff;
    }
    stats.stddev = sqrt(varianceSum / stats.totalDistinct);
    
    stats.highThreshold = stats.mean + 2 * stats.stddev;
    stats.lowThreshold = stats.mean - 2 * stats.stddev;
    
    return stats;
}

// 5. 分析并导出CSV
void FastqParser::analyzeAndExport(const unordered_map<uint64_t, int>& kmerFreq, 
                                   int k, int topN, 
                                   const string& filename) {
    if (kmerFreq.empty()) {
        cout << "No k-mers to analyze." << endl;
        return;
    }
    
    // 获取统计信息
    KmerStatistics stats = getStatistics(kmerFreq, k);
    
    cout << "\n=== K-mer Distribution Summary (k=" << k << ") ===" << endl;
    cout << "  Total distinct k-mers: " << stats.totalDistinct << endl;
    cout << "  Total occurrences: " << stats.totalOccurrences << endl;
    cout << "  Average frequency: " << stats.mean << endl;
    cout << "  Standard deviation: " << stats.stddev << endl;
    cout << "  Most frequent: " << stats.mostFreqKmer 
         << " (appears " << stats.maxCount << " times)" << endl;
    
    cout << "\n=== Anomaly Detection (threshold: mean ± 2σ) ===" << endl;
    cout << "  High threshold: " << stats.highThreshold << endl;
    cout << "  Low threshold: " << stats.lowThreshold << endl;
    
    // 收集异常
    vector<pair<uint64_t, int>> highAnomalies;
    vector<pair<uint64_t, int>> lowAnomalies;
    highAnomalies.reserve(kmerFreq.size() / 10);
    lowAnomalies.reserve(kmerFreq.size() / 10);
    
    for (const auto& pair : kmerFreq) {
        if (pair.second > stats.highThreshold) {
            highAnomalies.push_back(pair);
        } else if (pair.second < stats.lowThreshold) {
            lowAnomalies.push_back(pair);
        }
    }
    
    // 使用 std::sort (O(n log n))
    sort(highAnomalies.begin(), highAnomalies.end(),
         [](const auto& a, const auto& b) { return a.second > b.second; });
    sort(lowAnomalies.begin(), lowAnomalies.end(),
         [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // 打印高频异常
    cout << "\nTop " << topN << " most frequent k-mers:" << endl;
    size_t limit = min(highAnomalies.size(), static_cast<size_t>(topN));
    if (limit == 0) {
        cout << "  (No k-mers above high threshold)" << endl;
    } else {
        for (size_t i = 0; i < limit; ++i) {
            cout << "  " << (i+1) << ". " << decodeKmer(highAnomalies[i].first, k)
                 << ": " << highAnomalies[i].second << " times" << endl;
        }
    }
    
    // 打印低频异常
    cout << "\nTop " << topN << " least frequent k-mers:" << endl;
    limit = min(lowAnomalies.size(), static_cast<size_t>(topN));
    if (limit == 0) {
        cout << "  (No k-mers below low threshold)" << endl;
    } else {
        for (size_t i = 0; i < limit; ++i) {
            cout << "  " << (i+1) << ". " << decodeKmer(lowAnomalies[i].first, k)
                 << ": " << lowAnomalies[i].second << " times" << endl;
        }
    }
    
    // 导出CSV（包含高频和低频）
    string csvFilename = filename + "_k" + to_string(k) + "_top" + to_string(topN) + ".csv";
    ofstream outFile(csvFilename.c_str());
    
    if (outFile.is_open()) {
        outFile << "# metadata\n";
        outFile << "# filename: " << filename << "\n";
        outFile << "# k: " << k << "\n";
        outFile << "# topN: " << topN << "\n";
        outFile << "# total_distinct_kmers: " << stats.totalDistinct << "\n";
        outFile << "# total_occurrences: " << stats.totalOccurrences << "\n";
        outFile << "# mean: " << stats.mean << "\n";
        outFile << "# stddev: " << stats.stddev << "\n";
        outFile << "# high_threshold: " << stats.highThreshold << "\n";
        outFile << "# low_threshold: " << stats.lowThreshold << "\n";
        outFile << "# most_frequent_kmer: " << stats.mostFreqKmer << "\n";
        outFile << "category,kmer,frequency\n";
        
        // 写入高频异常
        limit = min(highAnomalies.size(), static_cast<size_t>(topN));
        for (size_t i = 0; i < limit; ++i) {
            outFile << "high," << decodeKmer(highAnomalies[i].first, k) 
                    << "," << highAnomalies[i].second << "\n";
        }
        
        // 写入低频异常
        limit = min(lowAnomalies.size(), static_cast<size_t>(topN));
        for (size_t i = 0; i < limit; ++i) {
            outFile << "low," << decodeKmer(lowAnomalies[i].first, k) 
                    << "," << lowAnomalies[i].second << "\n";
        }
        
        outFile.close();
        cout << "\n✅ Data saved to " << csvFilename << endl;
    } else {
        cerr << "Warning: Could not write " << csvFilename << endl;
    }
}