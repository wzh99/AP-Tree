#include "../include/aptree.hpp"
#include <limits>
#include <numeric>

struct APTree::QueryNested {
    Boundf region;
    std::vector<size_t> dict; // stores indexes of keywords in dictionary for algorithm efficiency
};

struct APTree::STObjectNested {
    Pointf location;
    std::vector<size_t> dict; // the same as QueryNested
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
    std::map<KeywordCut, std::vector<QueryNested *>> queries;
    std::vector<QueryNested *> dummy; // queries which has less queries than current offset
    float cost;
};

struct APTree::SpatialPartition {
    std::vector<float> partX, partY; // start position of X and Y partition, respectively
    std::vector<std::vector<QueryNested *>> cells; // row-major vector of query pointers in each cell
    std::vector<QueryNested *> dummy; // queries which covers the whole region
    float cost;
};

struct APTree::Node {
    enum NodeType { // assigned in callee code
        QUERY, // stores queries
        KEYWORD, // stores keyword cuts
        SPATIAL // stores spatial quadtree cells
    } type;

    Boundf bound; // assigned in caller code, for use of spatial partition

    // Package different type of node data in a union will cause problems
    struct QueryNode {
        std::vector<QueryNested> queries;
    };
    std::unique_ptr<QueryNode> query;

    struct KeywordNode {
        size_t offset; // the same with KeywordPartition
        std::vector<KeywordCut> cuts; // sorted keyword cuts for binary search
        std::map<KeywordCut, std::unique_ptr<Node>> children;
    };
    std::unique_ptr<KeywordNode> keyword;

    struct SpatialNode {
        std::vector<float> partX, partY;
        std::vector<std::unique_ptr<Node>> cells;
    };
    std::unique_ptr<SpatialNode> spatial;

    // Shared
    std::unique_ptr<Node> dummy = nullptr; // shared by query and spatial node, can be null

    Node() {}
    ~Node() {} // explicit destructor for std::unique_ptr
};

// For verifying selected queries, elements in two input vector must be sorted
template <class Type>
static bool isSubset(const std::vector<Type> &super, const std::vector<Type> &sub) {
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
        std::vector<size_t> indexVec;
        for (const auto &str : qry.keywords)
            indexVec.push_back(dictIndex[str]);
        std::sort(indexVec.begin(), indexVec.end()); // indexes must be sorted for verifying subset
        nestedQueries.push_back({qry.region, std::move(indexVec)});
    }

    // Build nested query pointer vector
    std::vector<QueryNested *> nstdQryPtrs;
    std::transform(nestedQueries.begin(), nestedQueries.end(), nstdQryPtrs.begin(),
                   [] (auto &qry) { return &qry; });

    // Start building tree
    root = std::make_unique<Node>();
    root->bound = Boundf({0, 0}, {1, 1}); // the whole region is a unit square
    build(root.get(), nstdQryPtrs, 0, true, true); // different from paper, offset is 0-indexed instead of 1-indexed
}

void APTree::build(Node *node, const std::vector<QueryNested *> &subQueries, size_t offset, bool useKeyword,
                   bool useSpatial)
{
    // Build query node
    if ((!useKeyword && !useKeyword) || subQueries.size() <= threshold) {
        node->type = Node::QUERY;
        node->query = std::make_unique<Node::QueryNode>();
        std::transform(subQueries.begin(), subQueries.end(), node->query->queries.begin(),
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
        spPart = spatialHeuristic(subQueries, node->bound);
        spCost = spPart.cost;
    }

    // Build keyword or partition node according to computed costs
    if (kwCost < spCost) { // keyword partition is chosen
        node->type = Node::KEYWORD;
        node->keyword = std::make_unique<Node::KeywordNode>();
        node->dummy = std::make_unique<Node>();
        node->keyword->offset = offset;
        build(node->dummy.get(), kwPart.dummy, offset + 1, false, useSpatial); 
            // keyword dummy node can no longer be partitioned by keyword again

        for (const auto &pair : kwPart.queries) {
            auto ptr = new Node();
            ptr->bound = node->bound; // keyword partition doesn't divide space
            node->keyword->cuts.push_back(pair.first);
            node->keyword->children[pair.first] = std::unique_ptr<Node>(ptr);
            build(ptr, pair.second, offset + 1, useKeyword, useSpatial);
        }

    } else { // spatial partition is chosen
        node->type = Node::SPATIAL;
        node->spatial = std::make_unique<Node::SpatialNode>();
        node->dummy = std::make_unique<Node>();
        build(node->dummy.get(), spPart.dummy, offset, useKeyword, false);
            // spatial dummy node can no longer be partitioned by space again
            
        size_t sizeX = node->spatial->partX.size(), sizeY = node->spatial->partY.size();
        const auto &partX = node->spatial->partX, &partY = node->spatial->partY; // shorter name for axis partition
        for (size_t j = 0; j < sizeY; j++)
            for (size_t i = 0; i < sizeX; i++) {
                size_t index = i + j * sizeX;
                auto ptr = new Node();
                ptr->bound = {{partX[i], (i == sizeX - 1) ? node->bound.max.x : partX[i + 1]},
                              {partY[i], (j == sizeY - 1) ? node->bound.max.y : partY[i + 1]}};
                node->spatial->cells[index] = std::unique_ptr<Node>(ptr);
                build(ptr, spPart.cells[index], offset, useKeyword, useSpatial);
            }
    }
}

APTree::KeywordPartition APTree::keywordHeuristic(const std::vector<QueryNested *> &subQueries, size_t offset)
{
    KeywordPartition partition;

    // Get statistics of queries to be partitioned
    std::map<size_t, std::vector<QueryNested *>> wordStatMap;
    std::vector<QueryNested *> dummy;
    for (const auto ptr : subQueries) {
        // Divide keyword into dummies and descendents
        if (ptr->dict.size() <= offset)
            dummy.push_back(ptr); // this query can't be partitioned by current offset
        else {
            size_t curWord = ptr->dict[offset];
            if (wordStatMap.find(curWord) == wordStatMap.end()) // new word to stats
                wordStatMap[curWord] = std::vector<QueryNested *>(); 
            wordStatMap[curWord].push_back(ptr); // count word frequency
        }
    }

    // Get total frequency of all appearing words
    size_t totalFreq = 0;
    for (const auto &pair : wordStatMap) totalFreq += pair.second.size();

    // Convert the word statistic map to vector for random access
    std::vector<std::pair<size_t, std::vector<QueryNested *>>> wordStatVec;
    wordStatVec.insert(wordStatVec.begin(), wordStatMap.begin(), wordStatMap.end());

    // Propose a partition which divides words by a fixed offset
    size_t nPart = std::min(nCuts, wordStatVec.size()); // the number of words to be partitioned may be less than nCuts
    size_t nWordPerPart = (wordStatVec.size() - 1) / nPart + 1;
    std::vector<KeywordCut> proposedCuts; 
    for (size_t i = 0; i < nPart; i++) {
        size_t cutStart = i * nWordPerPart; // cutStart and cutEnd are indexes in wordStatVec, not actual word indexes
        size_t cutEnd = std::min(wordStatVec.size(), cutStart + nWordPerPart) - 1;
        proposedCuts.push_back({wordStatVec[cutStart].first, wordStatVec[cutEnd].first});
    }

    // Find a partition which evenly partitions keywords by weight through heuristic algorithm
    for (size_t i = 1; i < nCuts; i++) {
    }

    // Prepare partition strategy to be returned
    partition.dummy = std::move(dummy); // dummy vector can no longer be used in this scope
    return partition;
}
