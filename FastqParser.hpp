#ifndef FASTQ_PARSER_HPP
#define FASTQ_PARSER_HPP
// 这两行是“防护措施”，防止这个头文件被不小心重复包含导致编译出错。是 C++ 的标准写法。

#include <string>
#include <unordered_map> // 新增：用于 K-mer 计数（哈希表）
#include <vector>        // 新增：用于存储序列（动态数组）

// 函数声明：统计 FASTQ 文件中的序列条数
// 声明一个叫 countReads 的函数：
// 它接收一个 const std::string& 类型的参数（就是文件名，而且保证不会修改它）
// 它返回一个 int（整数）
int countReads(const std::string& filename);

// 新增函数声明：从 FASTQ 文件中提取所有序列
// 它接收一个 const std::string& 类型的参数（文件名）
// 它返回一个 vector<string>（装着所有序列的“动态数组”）
std::vector<std::string> extractSequences(const std::string& filename);

// 新增函数声明：计算 K-mer 频率
// 它接收一个 const std::vector<std::string>&（所有序列的容器，而且保证不会修改它）
// 和一个 int（K-mer 的长度 k）
// 它返回一个 unordered_map<string, int>（哈希表，key=K-mer, value=出现次数）
std::unordered_map<std::string, int> computeKmerFrequencies(const std::vector<std::string>& sequences, int k);

#endif // 结束防护措施