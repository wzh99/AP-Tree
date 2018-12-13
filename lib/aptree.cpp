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
    float cost = 0;
};

struct APTree::SpatialPartition {
    std::vector<float> partX, partY; // start position of X and Y partition, respectively
    std::vector<std::vector<QueryNested *>> cells; // row-major vector of query pointers in each cell
    std::vector<QueryNested *> dummy; // queries which covers the whole region
    float cost = 0;
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
    std::unique_ptr<QueryNode> query = nullptr;

    struct KeywordNode {
        size_t offset; // the same with KeywordPartition
        std::vector<KeywordCut> cuts; // sorted keyword cuts for binary search
        std::map<KeywordCut, std::unique_ptr<Node>> children;

        KeywordCut Search(size_t word) const;
    };
    std::unique_ptr<KeywordNode> keyword = nullptr;

    struct SpatialNode {
        std::vector<float> partX, partY;
        std::vector<std::unique_ptr<Node>> cells;
    };
    std::unique_ptr<SpatialNode> spatial = nullptr;

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
    if ((!useKeyword && !useKeyword) || subQueries.size() < threshold) {
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

    // Get statistics of query keywords to be partitioned
    std::map<size_t, std::vector<QueryNested *>> offsetWordQryMap; // queries related to keywords of current offset
    std::map<size_t, size_t> allWordFreqMap; // frequencies related to keywords of all offsets
    std::vector<QueryNested *> dummy;
    for (const auto ptr : subQueries) {
        // Count keyword frequencies in all query offsets
        for (auto word : ptr->dict)
            if (allWordFreqMap.find(word) == allWordFreqMap.end())
                allWordFreqMap[word] = 1;
            else
                allWordFreqMap[word]++;

        // Divide keyword into dummies and descendents
        if (ptr->dict.size() <= offset)
            dummy.push_back(ptr); // this query can't be partitioned by current offset
        else {
            size_t curWord = ptr->dict[offset];
            if (offsetWordQryMap.find(curWord) == offsetWordQryMap.end()) // new word to offset keyword map
                offsetWordQryMap[curWord] = std::vector<QueryNested *>(); 
            offsetWordQryMap[curWord].push_back(ptr); // count word frequency
        }
    }

    // Get total frequency of all appearing words
    size_t totalFreq = 0;
    for (const auto &pair : allWordFreqMap) totalFreq += pair.second;

    // Convert the word statistic map to vector for random access
    std::vector<std::pair<size_t, std::vector<QueryNested *>>> offsetWordQryVec;
    offsetWordQryVec.insert(offsetWordQryVec.begin(), offsetWordQryMap.begin(), offsetWordQryMap.end());

    // Propose a partition which divides words by a fixed offset
    size_t nPart = std::min(nCuts, offsetWordQryVec.size()); // number of words to be partitioned may be less than nCuts
    size_t nWordPerPart = (offsetWordQryVec.size() - 1) / nPart + 1;
    std::vector<KeywordCut> propCuts; // proposed cuts
    propCuts.reserve(nPart);
    for (size_t i = 0; i < nPart; i++) {
        size_t cutStart = i * nWordPerPart; // cutStart and cutEnd are indexes in offsetWordQryVec, not actual word indexes
        size_t cutEnd = std::min(offsetWordQryVec.size(), cutStart + nWordPerPart) - 1;
        propCuts.push_back({offsetWordQryVec[cutStart].first, offsetWordQryVec[cutEnd].first});
    }

    // Defines cost computing function for each keyword cut
    auto computeCutCost = [&] (const KeywordCut &cut) {
        float weight = 0; // number of queries associated to current cut(bucket)
        float probSum = 0; // sum of probabilities of all words in current cut
        for (size_t word = cut.start; word <= cut.end; word++) {
            if (offsetWordQryMap.find(word) == offsetWordQryMap.end()) // word not found
                continue;
            weight += offsetWordQryMap[word].size();
            probSum += float(allWordFreqMap[word]) / totalFreq; // use local query distribution to estimate the whole one
        }
        return weight * probSum;
    };

    // Find a partition which evenly partitions keywords by weight through heuristic algorithm
    for (size_t i = 1; i < nPart; i++) {
        // Compute current cost of proposed cuts
        // Only account for two involved cuts, since other cuts has no change on total partition cost
        float curCost = computeCutCost(propCuts[i - 1]) + computeCutCost(propCuts[i]);

        // Try all possible cuts
        size_t leftCutStart = propCuts[i - 1].start, rightCutEnd = propCuts[i].end;
        for (size_t leftCutEnd = propCuts[i - 1].start; leftCutEnd != propCuts[i].end; leftCutEnd++) {
            if (offsetWordQryMap.find(leftCutEnd) == offsetWordQryMap.end()) // proposed left cut ending word not found
                continue;
            size_t rightCutStart = leftCutEnd + 1;
            while (offsetWordQryMap.find(rightCutStart) == offsetWordQryMap.end() && rightCutStart < rightCutEnd)
                rightCutStart++; // set start of right cut to valid index for efficiency
            float propCost = computeCutCost({leftCutStart, leftCutEnd}) + computeCutCost({rightCutStart, rightCutEnd});
            if (propCost < curCost) { // better cut
                curCost = propCost;
                propCuts[i - 1].end = leftCutEnd; // update proposed cuts
                propCuts[i].start = rightCutStart;
            }
        }
    }

    // Prepare partition strategy to be returned
    partition.cost = 0;
    partition.cuts = std::move(propCuts);
    partition.dummy = std::move(dummy); 
    for (const auto &cut : partition.cuts) {
        partition.cost += computeCutCost(cut);
        partition.queries[cut] = std::vector<QueryNested *>();
        auto &cutQryVec = partition.queries[cut];
        for (size_t word = cut.start; word <= cut.end; word++) { 
            if (offsetWordQryMap.find(word) == offsetWordQryMap.end()) continue;
            auto &wordQryVec = offsetWordQryMap[word];
            cutQryVec.insert(cutQryVec.end(), wordQryVec.begin(), wordQryVec.end());
        }
    }
    return partition;
}
