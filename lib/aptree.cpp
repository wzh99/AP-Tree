#include "../include/aptree.hpp"
#include <limits>

struct APTree::QueryNested {
    Boundf region;
    std::set<size_t> dict; // stores indexes of keywords in dictionart for algorithm efficiency
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
    std::map<KeywordCut, std::vector<QueryNested *>> cuts;
    std::vector<QueryNested *> dummy; // queries which has less queries than current offset
    float cost;
};

struct APTree::SpatialPartition {
    Boundf bound;
    std::vector<float> divX, divY; // start position of X and Y partition, respectively
    std::vector<std::vector<QueryNested *>> cells; // row-major vector of query pointers in each cell
    std::vector<QueryNested *> dummy; // queries which covers the whole region
    float cost;
};

struct APTree::Node {
    enum NodeType {
        QUERY, // stores queries
        KEYWORD, // stores keyword cuts
        SPATIAL // stores spatial quadtree cells
    } type;
    union { // package data in an union for different node types for sake of memory efficiency
        struct { // query node
            std::vector<QueryNested> queries;
        };
        struct { // keyword node
            size_t offset; // the same with KeywordPartition
            std::map<KeywordCut, std::unique_ptr<Node>> cuts;
        };
        struct { // spatial node
            Boundf bound; // the same with SpatialPartition
            std::vector<float> divX, divY;
            std::vector<std::unique_ptr<Node>> cells;
        };
    };
    std::unique_ptr<Node> dummy = nullptr; // shared by query and spatial node, can be null
    Node() {}
    ~Node() {} // explicit destructor for std::unique_ptr
};

// For verifying selected queries
template <class Type>
static bool isSubset(const std::set<Type> &super, const std::set<Type> &sub) {
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
    build(root.get(), nstdQryPtrs, 0, true, true); // different from paper, offset is 0-indexed instead of 1-indexed
}

void APTree::build(Node *node, const std::vector<QueryNested *> &subQueries, size_t offset, bool useKeyword,
                   bool useSpatial)
{
    // Build query node
    if ((!useKeyword && !useKeyword) || subQueries.size() < threshold) {
        node->type = Node::QUERY;
        std::transform(subQueries.begin(), subQueries.end(), node->queries.begin(),
                       [] (auto ptr) { return *ptr; }); // copy queries into node
        return;
    }

    // Compute cost for both keyword and spatial partitions
    float kwCost, spCost;
    kwCost = spCost = std::numeric_limits<float>::infinity();
    KeywordPartition kwPart;
    SpatialPartition spPart;
    if (useKeyword) { // try keyword partition
        kwPart = keywordHeuristic(subQueries, offset);
        kwCost = kwPart.cost;
    }
    if (useSpatial) { // try spatial partition
        spPart = spatialHeuristic(subQueries);
        spCost = spPart.cost;
    }

    // Build keyword or partition node according to computed costs
    if (kwCost < spCost) { // keyword partition is chosen
        node->type = Node::KEYWORD;
        node->offset = offset;
        node->dummy = std::make_unique<Node>();
        build(node->dummy.get(), kwPart.dummy, offset + 1, false, useSpatial); 
            // keyword dummy node can no longer be partitioned by keyword again
        for (const auto &pair : kwPart.cuts) {
            auto ptr = new Node();
            node->cuts[pair.first] = std::unique_ptr<Node>(ptr);
            build(ptr, pair.second, offset + 1, useKeyword, useSpatial);
        }
    } else { // spatial partition is chosen
        node->type = Node::SPATIAL;
        node->dummy = std::make_unique<Node>();
        build(node->dummy.get(), spPart.dummy, offset, useKeyword, false);
            // spatial dummy node can no longer be paritioned by space again
        size_t sizeX = node->divX.size(), sizeY = node->divY.size();
        for (size_t j = 0; j < sizeY; j++)
            for (size_t i = 0; i < sizeX; i++) {
                size_t index = i + j * sizeX;
                auto ptr = new Node();
                node->cells[index] = std::unique_ptr<Node>(ptr);
                build(ptr, spPart.cells[index], offset, useKeyword, useSpatial);
            }
    }
}
