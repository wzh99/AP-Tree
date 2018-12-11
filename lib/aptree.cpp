#include "../include/aptree.hpp"

struct APTree::QueryNested {
    Boundf region;
    std::set<size_t> dict; // stores indexes of keywords in dictionart for algorithm efficiency
    ~QueryNested() {}
};

struct APTree::STObjectNested {
    Pointf location;
    std::set<size_t> dict; // the same as QueryNested
};

struct APTree::KeywordCut {
    size_t start, end; // both ends included
    KeywordCut(size_t start, size_t end) : start(start), end(end) {}
    bool operator == (const KeywordCut &other) const { return start == other.start && end == other.end; }
    bool operator < (const KeywordCut &other) const { 
        if (start != other.start) return start < other.start; 
        return end < other.end;
    }
};

struct APTree::KeywordPartition {
    std::vector<KeywordCut> cuts;
    float cost;
};

struct APTree::SpatialPartition {
    std::vector<std::vector<std::vector<QueryNested *>>> cells; // m * n grid of query pointers
    float cost;
};

struct APTree::Node {
    enum NodeType {
        QUERY, // stores queries
        KEYWORD, // stores keyword cuts
        SPATIAL // stores spatial quadtree cells
    };
    union { // package data for different node types for memory efficiency
        struct { // query node
            std::vector<QueryNested> queries;
        };
        struct { // keyword node
            std::map<KeywordCut, std::unique_ptr<Node>> cuts;
        };
        struct { // spatial node
            Boundf bound;
            Pointf divide;
            std::vector<std::vector<std::unique_ptr<Node>>> cells; // an m * n array of pointer to child nodes
        };
    };
    std::unique_ptr<Node> dummy = nullptr; // shared by query and spatial node, can be null
    ~Node() {}
};

// For verifying selected queries
template <class Type>
static bool isSubset(const std::set<Type> &super, const std::set<Type> sub) {
    auto superIte = super.begin(), subIte = sub.begin();
    while (subIte != sub.end()) {
        if (superIte == super.end()) return false;
        auto superVal = *superIte, subVal = *subIte;
        if (superVal > subVal) return false;
        else if (superVal == subVal) {
            superIte++;
            subIte++;
        } else 
            while (*superIte < *subIte && superIte != super.end())
                superIte++;
    }
    return true;
}

APTree::APTree(const std::vector<std::string> &vocab, const std::vector<Query> &queries, size_t nCuts, size_t threshold)
    : dict(vocab.begin(), vocab.end()), nCuts(nCuts), threshold(threshold)
{
    // Build vocabulary index map
    for (size_t i = 0; i < dict.size(); i++)
        dictIndex[dict[i]] = i;

    // Build nested queries
    std::vector<QueryNested> nestedQueries;
    for (const auto &qry : queries) {
        std::set<size_t> indexes;
        for (const auto &str : qry.keywords)
            indexes.insert(dictIndex[str]);
        nestedQueries.push_back({qry.region, std::move(indexes)});
    }

    // Build nested query pointer vector
    std::vector<QueryNested *> nstdQryPtrs;
    std::transform(nestedQueries.begin(), nestedQueries.end(), nstdQryPtrs.begin(),
                   [] (auto &qry) { return &qry; });

    // Start building tree
    build(root.get(), nstdQryPtrs, 0, true, true);
}

void APTree::build(Node *node, const std::vector<QueryNested *> &subQueries, size_t offset, bool keyword, bool spatial)
{

}
