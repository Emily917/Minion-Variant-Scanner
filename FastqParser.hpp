#ifndef FASTQ_PARSER_HPP
#define FASTQ_PARSER_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>

// 统计信息结构体
struct KmerStatistics {
    size_t totalDistinct;
    uint64_t totalOccurrences;
    double mean;
    double stddev;
    uint64_t maxCount;
    std::string mostFreqKmer;
    double highThreshold;
    double lowThreshold;
};

class FastqParser {
public:
    FastqParser();
    ~FastqParser();

    // 统计read数量（自动检测3行/4行格式）
    int countReads(const std::string& filename);
    
    // 提取序列（支持3行格式）
    std::vector<std::string> extractSequences(const std::string& filename, size_t maxReads = 0);
    
    // 计算K-mer频率（使用整数编码 + 并行处理）
    std::unordered_map<uint64_t, int> computeKmerFrequencies(
        const std::vector<std::string>& sequences, int k);
    
    // 解码K-mer
    std::string decodeKmer(uint64_t encoded, int k);
    
    // 分析并导出CSV
    void analyzeAndExport(const std::unordered_map<uint64_t, int>& kmerFreq, 
                          int k, int topN, 
                          const std::string& filename);
    
    // 获取统计信息
    KmerStatistics getStatistics(const std::unordered_map<uint64_t, int>& kmerFreq, int k);

private:
    // 编码K-mer为64位整数
    uint64_t encodeKmer(const std::string& kmer);
    
    // 检查是否为有效的DNA序列
    bool isValidDNA(const std::string& seq);
    
    // 检测FASTQ格式（3行还是4行）
    bool detectFormat(const std::string& filename);
    
    // 线程池
    class ThreadPool;
    ThreadPool* pool_;
};

#endif