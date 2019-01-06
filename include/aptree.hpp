#pragma once

#include "sttree.hpp"
#include <map>
#include <memory>
#include <unordered_map>

class APTree : public STTree {
public:
    APTree(const std::vector<std::string> &vocab, const std::vector<Query> &queries, size_t f, size_t theta_Q, double theta_KL);
    ~APTree();
    std::vector<Query> Match(const STObject &obj) const override;
	void Register(const std::vector<Query> &qry) override;

private:
    struct QueryNested;
    struct STObjectNested;
    struct Node;
    struct KeywordCut;
    struct KeywordPartition;
    struct SpatialPartition;

    Node *root; // root node of the whole AP-Tree
    std::vector<std::string> dict; // stores all the vocabulary
    std::unordered_map<std::string, size_t> dictIndex; // stores index of keywords in dictionary
    const size_t f; // partition of f-ary tree node
    const size_t theta_Q; // maximum number of queries in a query node
	const double theta_KL; // KL-Divergence threshold for reconstruction

	// Index construction methods
    void build(Node *node, const std::vector<QueryNested *> &subQry);
    KeywordPartition keywordHeuristic(const std::vector<QueryNested *> &subQry, size_t offset);
    SpatialPartition spatialHeuristic(const std::vector<QueryNested *> &subQry, const Boundf &bound);

	// Query match method
	void match(const STObjectNested &obj, size_t offset, const Node *node, std::set<QueryNested> &out) const;

	// Query registration methods
	Node * regist(Node *node, const std::vector<QueryNested *> &newQry);
	std::vector<QueryNested *> collectAndMerge(const Node *node, const std::vector<QueryNested *> &newQry) const;
	void collect(const Node *node, std::vector<QueryNested *> &out) const;
};