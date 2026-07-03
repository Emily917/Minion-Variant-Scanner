#include <iostream>
#include "FastqParser.hpp" // 调用 FastqParser.hpp 文件
using namespace std;

int main() {
    cout << "=== Minion Variant Scanner ===" << endl; // 三代测序变异扫描器
    cout << "Step 1: Counting reads in FASTQ file..." << endl;

    string filename = "SRR390728_1.fastq"; // 创建一个叫 filename 的“盒子”，里面装着我们要读取的文件名。

    // 1. 统计序列条数（原有功能）
    int totalReads = countReads(filename);
    // 调用 countReads 这个函数（在 FastqParser.cpp 里定义的），把文件名传给它。
    // 它会返回一个整数（读到的序列条数），我们把它存到 totalReads 这个“盒子”里。

    if (totalReads >= 0) { // 判断一下：如果 totalReads 大于等于 0（说明读取没出错），就执行下面的代码。
        cout << "Total reads in " << filename << ": " << totalReads << endl;
    } else {
        return 1; // 如果出错了（返回 -1），就退出程序，返回 1 表示有错误
    }

    // 2. 提取所有序列
    cout << "\nStep 2: Extracting sequences..." << endl;
    vector<string> sequences = extractSequences(filename);
    // 调用 extractSequences 函数，把文件名传给它。
    // 它会返回一个装着所有序列的“容器”（vector），我们把它存到 sequences 这个“盒子”里。

    cout << "Extracted " << sequences.size() << " sequences." << endl;
    // sequences.size() 表示这个容器里有多少条序列，打印出来给用户看。

    // 3. 计算 K-mer 频率（以 k=11 为例）
    cout << "\nStep 3: Computing 11-mer frequencies..." << endl;
    int k = 11; // 设置 K-mer 的长度为 11
    unordered_map<string, int> kmerFreq = computeKmerFrequencies(sequences, k);
    // 调用 computeKmerFrequencies 函数，把 sequences（所有序列）和 k（长度）传给它。
    // 它会返回一个“哈希表”，里面记录了每个 K-mer 出现了几次。

    // 4. 展示一些统计结果
    cout << "Total distinct 11-mers: " << kmerFreq.size() << endl;
    // kmerFreq.size() 表示有多少种不同的 K-mer（去重后的数量）。

    // 找出出现次数最多的 K-mer
    string mostFreqKmer; // 用来存放出现次数最多的那个 K-mer
    int maxCount = 0; // 用来存放出现次数的最大值

    for (const auto& pair : kmerFreq) {
        // 遍历 kmerFreq 这个哈希表里的每一对（key = K-mer, value = 出现次数）
        // const auto& pair 是一种简写，意思是“每次取出一对，并且不修改它”
        if (pair.second > maxCount) {
            // pair.second 是“值”（出现次数），如果它比当前记录的最大值还大
            maxCount = pair.second; // 更新最大值
            mostFreqKmer = pair.first; // 记录下这个 K-mer（pair.first 是“键”）
        }
    }

    if (!mostFreqKmer.empty()) {
        // 如果 mostFreqKmer 不是空的（说明我们确实找到了出现次数最多的那个）
        cout << "Most frequent 11-mer: " << mostFreqKmer
             << " (appears " << maxCount << " times)" << endl;
    }

    return 0; // 程序正常结束，返回 0 给系统（表示一切顺利）
}