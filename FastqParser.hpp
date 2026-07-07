#ifndef FASTQ_PARSER_HPP
#define FASTQ_PARSER_HPP

#include <string>
#include <unordered_map>
#include <vector>

int countReads(const std::string& filename);
std::vector<std::string> extractSequences(const std::string& filename);
std::unordered_map<std::string, int> computeKmerFrequencies(const std::vector<std::string>& sequences, int k);
void analyzeKmerDistribution(const std::unordered_map<std::string, int>& kmerFreq, int k, int topN, const std::string& filename);

#endif
