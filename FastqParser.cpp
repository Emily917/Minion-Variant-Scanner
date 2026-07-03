#include "FastqParser.hpp"
#include <fstream>   // 用于文件输入输出
#include <iostream>  // 用于输入输出（cerr, endl）

using namespace std;

// --- 原有的 countReads 函数（保持不变） ---
int countReads(const string& filename) {
    // int: 函数执行完会返回一个整数（序列数量）。
    // const string& filename: 这个函数接收一个参数，类型是 string（字符串），名字叫 filename。
    // const: 保证函数内部不会修改这个字符串。
    // &: 引用传递，传的是文件的“地址”，而不是复制一份整个文件名，这样更高效。

    ifstream file(filename);
    // 创建一个名叫 file 的“输入文件流”对象，并用 filename 这个文件名去打开它。
    // 在 C++ 程序和磁盘上的 .fastq 文件之间，建立一条“读取管道”。

    if (!file.is_open()) { // 如果没有成功打开 file
        cerr << "Error: Cannot open file " << filename << endl;
        // 在标准错误流（cerr）中输出一行错误信息，方便排查问题。
        // endl：输出一个换行符。
        return -1; // 报错，返回 -1（因为正常的序列条数不可能是负数）
    }

    int readCount = 0; // 创建计数器，初始值为0
    string line; // 代表读取时的每一行内容，内容可被新行覆盖

    while (getline(file, line)) { // 开始读取每一行内容，一直循环直到最后一行
        // getline: 读取一整行
        // (file, line): 从 file 文件中读取一整行，存在 line 里

        if (!line.empty() && line[0] == '@') {
            // !line.empty(): 这一行不是空的（防止读到空行时访问 line[0] 导致错误）
            // line[0] == '@': 这一行的第一个字符是 '@'（这个是基于标准的 FASTQ 文件标识行开头格式）
            readCount++; // 计数器加1
        }
    }

    file.close(); // 关闭之前打开的文件，释放文件资源
    return readCount; // 返回计数器的最终值（总序列条数）
}


// --- 新增：提取所有序列 ---
vector<string> extractSequences(const string& filename) {
    // 这个函数的作用：打开一个 FASTQ 文件，把每条序列的碱基行提取出来，装进一个容器里返回。

    vector<string> sequences; // 用来存放所有序列的“容器”（动态数组）
    // 创建一个空的 vector，叫 sequences，里面装的都是 string。

    ifstream file(filename); // 创建一个“文件流”对象，尝试打开文件

    if (!file.is_open()) { // 如果文件没打开
        cerr << "Error: Cannot open file " << filename << endl; // 输出错误信息
        return sequences; // 返回空容器（因为没法读数据了）
    }

    string line; // 临时存放每一行内容
    int lineCount = 0; // 记录当前是第几行（从 0 开始计数）

    while (getline(file, line)) { // 逐行读取文件内容
        // FASTQ 格式：每 4 行一组，第 1 行是标识行（以 @ 开头），第 2 行是序列（碱基行），
        // 第 3 行是分隔行（以 + 开头），第 4 行是质量行。
        // 所以，序列行在每组里的位置是：行号 % 4 == 1（0-based 计数）
        // 也就是第 0 行是 @...，第 1 行是 ATGC...，第 2 行是 +...，第 3 行是质量值。

        if (lineCount % 4 == 1) {
            // 如果当前行号除以 4 的余数是 1，说明这一行是序列行（碱基行）
            sequences.push_back(line);
            // push_back 意思是“把它放到容器的末尾”
            // 这一行代码把当前读到的序列（line）添加到 sequences 这个容器的最后面
        }
        lineCount++; // 行号加1，继续读下一行
    }

    file.close(); // 关闭文件
    return sequences; // 把装满了所有序列的容器返回给调用者
}


// --- 新增：计算 K-mer 频率 ---
unordered_map<string, int> computeKmerFrequencies(const vector<string>& sequences, int k) {
    // 这个函数的作用：把所有序列里的 K-mer 都切出来，统计每个 K-mer 出现了几次。
    // 参数：sequences（所有序列的容器），k（K-mer 的长度）
    // 返回：一个 unordered_map（哈希表），key=K-mer 字符串，value=出现次数

    unordered_map<string, int> kmerCount;
    // 创建一个空的哈希表，叫 kmerCount
    // key：K-mer 的序列（比如 "ATCG"），value：这个 K-mer 出现的次数（整数）

    for (const string& seq : sequences) {
        // 遍历 sequences 这个容器里的每一条序列
        // const string& seq 意思是“每次从容器里取出一条序列，放到 seq 里，并且不修改它”
        // 就像从盒子里一个一个拿出小球来处理

        if (seq.length() < k) {
            // 如果这条序列的长度比 k 还要短，那就没法切出 K-mer，直接跳过
            continue; // 跳过这一条，处理下一条
        }

        // 滑动窗口：从位置 0 开始，一直滑动到 seq.length() - k 的位置
        // 比如序列 "ATCGATCG"，k=3，那么窗口会从位置 0,1,2,3,4,5 开始切
        // 分别切出 "ATC", "TCG", "CGA", "GAT", "ATC", "TCG"
        for (size_t i = 0; i <= seq.length() - k; ++i) {
            // size_t 是一种无符号整数类型，专门用来表示“大小”和“位置”
            // i <= seq.length() - k：确保窗口不超出序列范围

            string kmer = seq.substr(i, k);
            // substr(i, k) 是 C++ 字符串的“截取”方法
            // 意思是从位置 i 开始，截取长度为 k 的子串
            // 比如 seq = "ATCGATCG"，i=2，k=3，那么 kmer = "CGA"

            kmerCount[kmer]++;
            // 这一行是“哈希表”的精髓！
            // 如果 kmer 这个 key 已经在哈希表里了，就把它的 value 加 1
            // 如果还没在哈希表里，就把它加进去，value 初始为 0，然后加 1 变成 1
            // 相当于：找到这个 K-mer 对应的“计数器”，然后让它加 1
        }
    }

    return kmerCount; // 把统计好的哈希表返回给调用者
}