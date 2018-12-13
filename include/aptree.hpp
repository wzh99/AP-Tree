#pragma once

#include "sttree.hpp"
#include <map>
#include <memory>
#include <unordered_map>

class APTree : public STTree {
public:
    APTree(const std::vector<std::string> &vocab, const std::vector<Query> &queries, size_t nCuts, size_t threshold);
    std::vector<Query> Match(const STObject &obj) const override;

private:
    struct QueryNested;
    struct STObjectNested;
    struct Node;
    struct KeywordCut;
    struct KeywordPartition;
    struct SpatialPartition;

    std::unique_ptr<Node> root; // root node of the whole AP-Tree
    std::vector<std::string> dict; // stores all the vocabulary
    std::unordered_map<std::string, size_t> dictIndex; // stores index of keywords in dictionart
    const size_t nCuts; // partition of f-ary tree node
    const size_t threshold; // maximum number of queries in a query node

    void match(const STObjectNested &obj, size_t start, const Node *node, std::set<QueryNested *> &out) const;
    void build(Node *node, const std::vector<QueryNested *> &subQueries, size_t offset, bool useKeyword, bool useSpatial);
    KeywordPartition keywordHeuristic(const std::vector<QueryNested *> &subQueries, size_t offset);
    SpatialPartition spatialHeuristic(const std::vector<QueryNested *> &subQueries, const Boundf &bound);
};