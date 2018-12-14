#include "../include/aptree.hpp"
#include <cmath>
#include <limits>
#include <numeric>

static constexpr auto INDEX_NOT_FOUND = std::numeric_limits<size_t>::max();

// For indicating range of space or keyword
template <class Type>
static size_t rangeSearch(const std::vector<Type> &ranges, const Type target) {
    size_t start = 0, end = ranges.size() - 1, mid = (start + end) / 2;
    if (target < ranges[start] || target >= ranges[end]) 
        return INDEX_NOT_FOUND; // can't be in any range
    while (true) {
        if (end - start == 1) return start; // range is obvious
        mid = (start + end) / 2;
        auto startEle = ranges[start], midEle = ranges[mid], endEle = ranges[end];
        if (target == startEle) return start;
        if (target == midEle) return mid;
        if (target == endEle) return end;
        if (target > startEle && target < midEle) {
            end = mid;
            continue;
        } else if (target > midEle && target < endEle) {
            start = mid;
            continue;
        }
    }
    return INDEX_NOT_FOUND;
}

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

// For selecting common queries from two ranges in corresponding axis
template <class Type>
static std::vector<Type> commonElements(const std::vector<Type> &v1, const std::vector<Type> &v2) {
    std::vector<Type> common;
    if (v1.size() == 0 || v2.size() == 0) return common;
    auto ite1 = v1.begin(), ite2 = v2.begin();
    while (ite1 != v1.end() && ite2 != v2.end()) {
        if (*ite1 == *ite2) {
            common.push_back(*ite1);
            ite1++;
            ite2++;
        }
        else if (*ite1 < *ite2)
            ite1++;
        else
            ite2++;
    }
    return std::move(common);
}

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
    double cost = 0;
};

struct APTree::SpatialPartition {
    size_t nPartX, nPartY;
    std::vector<double> partX, partY; // start position of X and Y partition, respectively
    std::vector<std::vector<QueryNested *>> cells; // column-major vector of query pointers in each cell
    std::vector<QueryNested *> dummy; // queries which covers the whole region
    double cost = 0;
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
        size_t nPartX, nPartY;
        std::vector<double> partX, partY; // has one more element than nPart
        std::vector<std::unique_ptr<Node>> cells; // 1D vector of m * n cells [0 ... m-1][0 ... n-1]

        Pointu GetCellIndex(const Pointf &pt) const;
    };
    std::unique_ptr<SpatialNode> spatial = nullptr;

    // Shared
    std::unique_ptr<Node> dummy = nullptr; // shared by query and spatial node, can be null

    Node() {}
    ~Node() {} // explicit destructor for std::unique_ptr
};

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
    double kwCost, spCost;
    kwCost = spCost = std::numeric_limits<double>::infinity();
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

        node->keyword->cuts = kwPart.cuts;
        for (const auto &pair : kwPart.queries) {
            auto ptr = new Node();
            ptr->bound = node->bound; // keyword partition doesn't divide space
            node->keyword->children[pair.first] = std::unique_ptr<Node>(ptr);
            build(ptr, pair.second, offset + 1, useKeyword, useSpatial);
        }

    } else { // spatial partition is chosen
        node->type = Node::SPATIAL;
        node->spatial = std::make_unique<Node::SpatialNode>();
        node->dummy = std::make_unique<Node>();
        build(node->dummy.get(), spPart.dummy, offset, useKeyword, false);
            // spatial dummy node can no longer be partitioned by space again
            
        node->spatial->nPartX = spPart.nPartX;
        node->spatial->nPartY = spPart.nPartY;
        node->spatial->partX = spPart.partX;
        node->spatial->partY = spPart.partY;
        for (size_t i = 0; i < spPart.nPartX; i++)
            for (size_t j = 0; j < spPart.nPartY; j++) {
                size_t index = j + i * spPart.nPartY;
                auto ptr = new Node();
                ptr->bound = {{spPart.partX[i], spPart.partY[j]}, {spPart.partX[i + 1], spPart.partY[j + 1]}};
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
        double weight = 0; // number of queries associated to current cut(bucket)
        double freqSum = 0; // sum of probabilities of all words in current cut
        for (size_t word = cut.start; word <= cut.end; word++) {
            if (offsetWordQryMap.find(word) == offsetWordQryMap.end()) // word not found
                continue;
            weight += offsetWordQryMap[word].size();
            freqSum += allWordFreqMap[word]; // use local query distribution to estimate the whole one
        }
        return weight * freqSum / totalFreq;
    };

    // Find a partition which evenly partitions keywords by weight through heuristic algorithm
    for (size_t i = 1; i < nPart; i++) {
        // Compute current cost of proposed cuts
        // Only account for two involved cuts, since other cuts has no change on total partition cost
        double curCost = computeCutCost(propCuts[i - 1]) + computeCutCost(propCuts[i]);

        // Try all possible cuts
        size_t leftCutStart = propCuts[i - 1].start, rightCutEnd = propCuts[i].end;
        for (size_t leftCutEnd = leftCutStart; leftCutEnd != rightCutEnd; leftCutEnd++) {
            if (offsetWordQryMap.find(leftCutEnd) == offsetWordQryMap.end()) // proposed left cut ending word not found
                continue;
            size_t rightCutStart = leftCutEnd + 1;
            while (offsetWordQryMap.find(rightCutStart) == offsetWordQryMap.end() && rightCutStart < rightCutEnd)
                rightCutStart++; // set start of right cut to valid index for efficiency
            double propCost = computeCutCost({leftCutStart, leftCutEnd}) + computeCutCost({rightCutStart, rightCutEnd});
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

APTree::SpatialPartition APTree::spatialHeuristic(const std::vector<QueryNested *> &subQueries, const Boundf &bound)
{
    SpatialPartition partition;
    size_t nPartX = partition.nPartX = std::floor(std::sqrt(nCuts));
    size_t nPartY = partition.nPartY = nCuts / nPartX;

    // Prune dummy queries
    std::vector<QueryNested *> cellQry, dummy;
    for (const auto qry : subQueries)
        if (qry->region.Contains(bound))
            dummy.push_back(qry);
        else
            cellQry.push_back(qry);
    partition.dummy = std::move(dummy);
    
    std::map<double, std::vector<QueryNested *>> bndQryMapX, bndQryMapY;
    { // X axis
        // Get statistics of all query bounds
        std::set<double> bndStatSet;
        for (const auto qry : cellQry) {
            bndStatSet.insert(qry->region.min.x);
            bndStatSet.insert(qry->region.max.x);
        }
        bndStatSet.insert(bound.max.x); 
        std::vector<double> bndStatVec;
        bndStatVec.insert(bndStatVec.end(), bndStatSet.begin(), bndStatSet.end());

        // Get statistics of queries related to each bound
        for (size_t i = 0; i < bndStatVec.size() - 1; i++) {
            auto &qryVec = bndQryMapX[bndStatVec[i]] = std::vector<QueryNested *>();
            if (i == bndStatSet.size() - 1) continue;
            Boundf curBnd = {{bndStatVec[i], bound.min.y}, {bndStatVec[i + 1], bound.max.y}};
            for (const auto qry : subQueries)
                if (curBnd.Overlaps(qry->region))
                    qryVec.push_back(qry);
        }

        // Convert query map to vector for random access
        std::vector<std::pair<double, std::vector<QueryNested *>>> bndQryVec;
        bndQryVec.insert(bndQryVec.end(), bndQryMapX.begin(), bndQryMapX.end());

        // Propose a partition which divides queries by a fixed bound number offset
        size_t nBndPerPart = (bndStatVec.size() - 1) / nPartX;
        std::vector<size_t> partIdx; // indexes in bound statistics vector
        partIdx.reserve(nPartX + 1);
        for (size_t i = 0; i < nPartX; i++) 
            partIdx[i] = i * nBndPerPart;
        partIdx[nPartX] = bndStatVec.size() - 1;

        // Define cost function for each proposed axis partition
        auto computeCost = [&] (size_t thisIndex, size_t nextIndex) {
            size_t weight = 0;
            for (size_t i = thisIndex; i < nextIndex; i++)
                weight += bndQryVec[i].second.size();
            double length = bndStatVec[nextIndex] - bndStatVec[thisIndex];
            return weight * length / (bound.max.x - bound.min.x);
        };

        // Find a partition which evenly partitions X axis by weight through heuristic algorithm
        for (size_t i = 1; i < nPartX; i++) {
            // Compute current cost of proposed partition
            double curCost = computeCost(partIdx[i - 1], partIdx[i]) + computeCost(partIdx[i], partIdx[i + 1]);
            // Try all possible partition
            size_t leftPartIdx = partIdx[i - 1], nextPartIdx = partIdx[i + 1];
            for (size_t curPartIdx = leftPartIdx + 1; curPartIdx < nextPartIdx; curPartIdx++) {
                double newCost = computeCost(leftPartIdx, curPartIdx) + computeCost(curPartIdx, nextPartIdx);
                if (newCost < curCost) { // better partition
                    curCost = newCost;
                    partIdx[i] = curPartIdx;
                }
            }
        }

        // Convert index into real number partitions
        std::transform(partIdx.begin(), partIdx.end(),
                       partition.partX.begin(), [&](size_t i) { return bndStatVec[i]; });
    }
    { // Y axis
        // Get statistics of all query bounds
        std::set<double> bndStatSet;
        for (const auto qry : cellQry) {
            bndStatSet.insert(qry->region.min.y);
            bndStatSet.insert(qry->region.max.y);
        }
        bndStatSet.insert(bound.max.y); 
        std::vector<double> bndStatVec;
        bndStatVec.insert(bndStatVec.end(), bndStatSet.begin(), bndStatSet.end());

        // Get statistics of queries related to each bound
        for (size_t i = 0; i < bndStatVec.size() - 1; i++) {
            auto &qryVec = bndQryMapY[bndStatVec[i]] = std::vector<QueryNested *>();
            if (i == bndStatSet.size() - 1) continue;
            Boundf curBnd = {{bound.min.x, bndStatVec[i]}, {bound.max.x, bndStatVec[i + 1]}};
            for (const auto qry : subQueries)
                if (curBnd.Overlaps(qry->region))
                    qryVec.push_back(qry);
        }

        // Convert query map to vector for random access
        std::vector<std::pair<double, std::vector<QueryNested *>>> bndQryVec;
        bndQryVec.insert(bndQryVec.end(), bndQryMapY.begin(), bndQryMapY.end());

        // Propose a partition which divides queries by a fixed bound number offset
        size_t nBndPerPart = (bndStatVec.size() - 1) / nPartX;
        std::vector<size_t> partIdx; // indexes in bound statistics vector
        partIdx.reserve(nPartY + 1);
        for (size_t i = 0; i < nPartY; i++) 
            partIdx[i] = i * nBndPerPart;
        partIdx[nPartY] = bndStatVec.size() - 1;

        // Define cost function for each proposed axis partition
        auto computeCost = [&] (size_t thisIndex, size_t nextIndex) {
            size_t weight = 0;
            for (size_t i = thisIndex; i < nextIndex; i++)
                weight += bndQryVec[i].second.size();
            double length = bndStatVec[nextIndex] - bndStatVec[thisIndex];
            return weight * length / (bound.max.y - bound.min.y);
        };

        // Find a partition which evenly partitions X axis by weight through heuristic algorithm
        for (size_t i = 1; i < nPartY; i++) {
            // Compute current cost of proposed partition
            double curCost = computeCost(partIdx[i - 1], partIdx[i]) + computeCost(partIdx[i], partIdx[i + 1]);
            // Try all possible partition
            size_t leftPartIdx = partIdx[i - 1], nextPartIdx = partIdx[i + 1];
            for (size_t curPartIdx = leftPartIdx + 1; curPartIdx < nextPartIdx; curPartIdx++) {
                double newCost = computeCost(leftPartIdx, curPartIdx) + computeCost(curPartIdx, nextPartIdx);
                if (newCost < curCost) { // better partition
                    curCost = newCost;
                    partIdx[i] = curPartIdx;
                }
            }
        }

        // Convert index into real number partitions
        std::transform(partIdx.begin(), partIdx.end(),
                       partition.partY.begin(), [&](size_t i) { return bndStatVec[i]; });
    }

    // Generate cells and compute final cost
    partition.cost = 0;
    for (size_t i = 0; i < nPartX; i++) {
        for (size_t j = 0; j < nPartY; j++) {
            size_t index = j + i * nPartY;
            auto &qryVecX = bndQryMapX[partition.partX[i]], &qryVecY = bndQryMapY[partition.partY[j]];
            std::sort(qryVecX.begin(), qryVecX.end());
            std::sort(qryVecY.begin(), qryVecY.end());
            auto qryCell = commonElements(qryVecX, qryVecY);
            size_t weight = qryCell.size();
            double prob = (partition.partX[i + 1] - partition.partX[i]) * (partition.partY[i + 1] - partition.partY[i]);
            partition.cost += weight * prob;
            partition.cells[index] = std::move(qryCell);
        }
    }

    return partition;
}
