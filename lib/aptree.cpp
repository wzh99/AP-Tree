#include "../include/aptree.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

#define COUT(expr) std::cout << expr << std::endl;

static constexpr auto INDEX_NOT_FOUND = std::numeric_limits<size_t>::max();

// For indicating range of space or keyword
// Array pointer version
template <class Type>
static size_t rangeSearch(const Type *ranges, size_t size, const Type target) noexcept {
    size_t start = 0, end = size - 1, mid = (start + end) / 2;
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
        }
        else if (target > midEle && target < endEle) {
            start = mid;
            continue;
        }
    }
    return INDEX_NOT_FOUND;
}

// std::vector version
template <class Type>
static inline size_t rangeSearch(const std::vector<Type> &ranges, const Type target) noexcept {
    return rangeSearch(ranges.data(), ranges.size(), target);
}

// For verifying selected queries, elements in two input vector must be sorted
template <class Type>
static bool isSubset(const std::vector<Type> &super, const std::vector<Type> &sub) noexcept {
    auto superIte = super.begin(), subIte = sub.begin();
    while (subIte != sub.end()) {
        if (superIte == super.end()) return false;
        auto superVal = *superIte, subVal = *subIte;
        if (superVal > subVal) return false;
        else if (superVal == subVal) {
            superIte++;
            subIte++;
        } else 
            while (superIte != super.end() && *superIte < *subIte)
                superIte++;
    }
    return true;
}

// For selecting common queries from two ranges in corresponding axis
template <class Type>
static std::vector<Type> commonElements(const std::set<Type> &v1, const std::set<Type> &v2) noexcept {
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
    return common;
}

struct APTree::QueryNested {
    Boundf region;
    std::vector<size_t> keywords; // stores indexes of keywords in dictionary for algorithm efficiency
    bool operator < (const  QueryNested &qry) const noexcept { return region.Area() < qry.region.Area(); }
    bool operator == (const QueryNested &qry) const noexcept { return region == qry.region && keywords == qry.keywords; }
};

struct APTree::STObjectNested {
    Pointf location;
    std::vector<size_t> keywords; // the same as QueryNested
};

struct APTree::KeywordCut {
    size_t start, end; // both ends included
    KeywordCut() {}
    KeywordCut(size_t start, size_t end) : start(start), end(end) {}
    bool isInCut(size_t index) const noexcept { return index >= start && index <= end; }
    bool operator == (const KeywordCut &other) const noexcept { return start == other.start && end == other.end; }
    bool operator < (const KeywordCut &other) const noexcept { 
        if (start != other.start) return start < other.start; 
        return end < other.end;
    }
};

struct APTree::KeywordPartition {
    size_t nPart;
    std::unique_ptr<KeywordCut[]> cuts;
    std::unique_ptr<std::vector<QueryNested *>[]> queries;
    std::vector<QueryNested *> dummy; // queries which has less queries than current offset
    double cost = 0;
};

struct APTree::SpatialPartition {
    size_t nPartX, nPartY;
    std::unique_ptr<double[]> partX, partY; // start position of X and Y partition, respectively
    std::unique_ptr<std::vector<QueryNested *>[]> cells; // column-major vector of query pointers in each cell
    std::vector<QueryNested *> dummy; // queries which covers the whole region
    double cost = 0;
};

struct APTree::Node {
    enum NodeType { // assigned in callee code
        QUERY, // stores queries
        KEYWORD, // stores keyword cuts
        SPATIAL // stores spatial quadtree cells
    } type;

    size_t offset; // assigned in caller code, for use of keyword partition
    Boundf bound; // assigned in caller code, for use of spatial partition
    bool useKw, useSp; // assigned in caller code, instead of arguments passed in build method

    // Package different type of node data in a union will cause problems
    struct QueryNode {
        std::vector<QueryNested> queries;
    };
    std::unique_ptr<QueryNode> query = nullptr;

    struct KeywordNode {
        size_t offset; // the same with KeywordPartition
        size_t nPart;
        std::unique_ptr<KeywordCut[]> cuts; // sorted keyword cuts for binary search
        std::unique_ptr<std::unique_ptr<Node>[]> children;
        std::unique_ptr<size_t[]> nOld; // record number of descendent queries in last construction and 
        std::unique_ptr<size_t[]> nAdd; // record number of newly added queries since last construction

        // Find corresponding cut based on object keywords, return index of cut found
        // The algorithm is similar to rangeSearch(), but deals with keyword cuts instead of ranges
        size_t Search(size_t target) const noexcept {
            size_t start = 0, end = nPart - 1, mid = (start + end) / 2;
            if (target < cuts[start].start || target > cuts[end].end)
                return INDEX_NOT_FOUND; // can't be in any cut
            while (true) {
                if (end - start == 1) {
                    if (cuts[start].isInCut(target))
                        return start;
                    else 
                        return INDEX_NOT_FOUND; // between gap of two cuts
                }
                mid = (start + end) / 2;
                auto startCut = cuts[start], midCut = cuts[mid], endCut = cuts[end];
                if (startCut.isInCut(target)) return start;
                if (midCut.isInCut(target)) return mid;
                if (endCut.isInCut(target)) return end;
                if (target >= startCut.start && target < midCut.start) {
                    end = mid;
                    continue;
                } else if (target >= midCut.start && target < endCut.start) {
                    start = mid;
                    continue;
                }
            }
            return INDEX_NOT_FOUND;
        }
    };
    std::unique_ptr<KeywordNode> keyword = nullptr;

    struct SpatialNode {
        size_t nPartX, nPartY;
        std::unique_ptr<double[]> partX, partY; // has one more element than nPart
        std::unique_ptr<std::unique_ptr<Node>[]> cells; // 1D vector of m * n cells [0 ... m-1][0 ... n-1]
        std::unique_ptr<size_t[]> nOld, nAdd; // the same as 

        Pointu GetCellIndex(const Pointf &pt) const noexcept {
            size_t ix = rangeSearch(partX.get(), nPartX, pt.x);
            size_t iy = rangeSearch(partY.get(), nPartY, pt.y);
            return {ix, iy};
        }
    };
    std::unique_ptr<SpatialNode> spatial = nullptr;

    // Shared
    std::unique_ptr<Node> dummy = nullptr; // shared by query and spatial node, can be null

    Node() {}
    Node(const Node &node) 
        : offset(node.offset), bound(node.bound), useKw(node.useKw), useSp(node.useSp) {}
    Node(size_t offset, const Boundf &bd, bool useKw, bool useSp)
        : offset(offset), bound(bd), useKw(useKw), useSp(useSp) {}
    ~Node() {} // explicit destructor for std::unique_ptr
};

APTree::APTree(const std::vector<std::string> &vocab, const std::vector<Query> &queries, size_t f, size_t theta_Q, double theta_KL)
    : dict(vocab.begin(), vocab.end()), f(f), theta_Q(theta_Q), theta_KL(theta_KL)
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

    // Remove repeated queries
    std::sort(nestedQueries.begin(), nestedQueries.end());
    auto newEnd = std::unique(nestedQueries.begin(), nestedQueries.end());
    nestedQueries.erase(newEnd, nestedQueries.end());

    // Build nested query pointer vector
    std::vector<QueryNested *> nstdQryPtrs;
    auto qryArrPtr = nestedQueries.data();
    for (size_t i = 0; i < nestedQueries.size(); i++) 
        nstdQryPtrs.push_back(qryArrPtr + i);

    // Start building tree
    root = new Node(0, Boundf({ 0, 0 }, { 1, 1 }), true, true);
    build(root, nstdQryPtrs); // different from paper, offset is 0-indexed instead of 1-indexed
}

APTree::~APTree() { delete root; }

void APTree::build(Node *node, const std::vector<QueryNested *> &subQry)
{
    // Build query node
    if ((!node->useKw && !node->useSp) || subQry.size() < theta_Q) {
        node->type = Node::QUERY;
        node->query = std::make_unique<Node::QueryNode>();
        // Copy queries into node
        auto &queries = node->query->queries;
        for (const auto ptr : subQry)
            queries.push_back(*ptr);
        return;
    }

    // Compute cost for both keyword and spatial partitions
    double kwCost, spCost;
    kwCost = spCost = std::numeric_limits<double>::infinity();
    KeywordPartition kwPart;
    SpatialPartition spPart;
    if (node->useKw) { // try keyword partition
        kwPart = keywordHeuristic(subQry, node->offset);
        kwCost = kwPart.cost;
    }
    if (node->useSp) { // try spatial partition
        spPart = spatialHeuristic(subQry, node->bound);
        spCost = spPart.cost;
    }

    // Build keyword or partition node according to computed costs
    if (kwCost < spCost) { // keyword partition is chosen
        node->type = Node::KEYWORD;
        node->keyword = std::make_unique<Node::KeywordNode>();
        if (kwPart.dummy.size() > 0) {
            node->dummy = std::make_unique<Node>(node->offset + 1, node->bound, false, node->useSp);
            build(node->dummy.get(), kwPart.dummy);
            // keyword dummy node can no longer be partitioned by keyword again
        }

        size_t nPart = node->keyword->nPart = kwPart.nPart;
        node->keyword->cuts = std::move(kwPart.cuts);
        node->keyword->children = std::make_unique<std::unique_ptr<Node>[]>(nPart);
        node->keyword->nOld = std::make_unique<size_t[]>(nPart);
        node->keyword->nAdd = std::make_unique<size_t[]>(nPart);
        memset(node->keyword->nAdd.get(), 0, nPart * sizeof(size_t)); // set nOld to zeros
        for (size_t i = 0; i < nPart; i++) {
            auto ptr = new Node(node->offset + 1, node->bound, node->useKw, node->useSp);
            node->keyword->children[i] = std::unique_ptr<Node>(ptr);
            node->keyword->nOld[i] = kwPart.queries[i].size();
            build(ptr, kwPart.queries[i]);
        }

    } else { // spatial partition is chosen
        node->type = Node::SPATIAL;
        node->spatial = std::make_unique<Node::SpatialNode>();
        if (spPart.dummy.size() > 0) {
            node->dummy = std::make_unique<Node>(node->offset, node->bound, node->useKw, false);
            build(node->dummy.get(), spPart.dummy);
            // spatial dummy node can no longer be partitioned by space again
        }
            
        // Move data from partition strategy to node
        auto nPartX = node->spatial->nPartX = spPart.nPartX;
        auto nPartY = node->spatial->nPartY = spPart.nPartY;
        auto nPart = nPartX * nPartY;
        node->spatial->partX = std::move(spPart.partX);
        node->spatial->partY = std::move(spPart.partY);
        node->spatial->cells = std::make_unique<std::unique_ptr<Node>[]>(nPart);
        node->spatial->nOld = std::make_unique<size_t[]>(nPart);
        node->spatial->nAdd = std::make_unique<size_t[]>(nPart);
        memset(node->spatial->nAdd.get(), 0, nPart * sizeof(size_t)); // set nAdd to zeros
        auto partX = node->spatial->partX.get(), partY = node->spatial->partY.get();
        for (size_t i = 0; i < nPartX; i++)
            for (size_t j = 0; j < nPartY; j++) {
                size_t index = j + i * nPartY;
                auto bound = Boundf{ {partX[i], partY[j]}, {partX[i + 1], partY[j + 1]} };
                auto ptr = new Node(node->offset, bound, node->useKw, node->useSp);
                node->spatial->cells[index] = std::unique_ptr<Node>(ptr);
                node->spatial->nOld[index] = spPart.cells[index].size();
                build(ptr, spPart.cells[index]);
            }
    }
}

// Bucket problem for heuristic algorithm
// Put m queries in n buckets. Suppose buckets have the same capacity, what is the minimum capacity required?
inline static size_t bucketCapacity(size_t m, size_t n) noexcept {
    return m % n == 0 ? (m / n) : (m / n + 1);
}

APTree::KeywordPartition APTree::keywordHeuristic(const std::vector<QueryNested *> &subQry, size_t offset)
{
    KeywordPartition partition;

    // Get statistics of query keywords to be partitioned
    std::map<size_t, std::vector<QueryNested *>> offsetWordQryMap; // queries related to keywords of current offset
    std::map<size_t, size_t> allWordFreqMap; // frequencies related to keywords of all offsets
    std::vector<QueryNested *> dummy;
    double dictSize = 0;
    for (const auto ptr : subQry) {
        // Count keyword frequencies in all query offsets
        dictSize += ptr->keywords.size();
        for (auto word : ptr->keywords)
            if (allWordFreqMap.find(word) == allWordFreqMap.end())
                allWordFreqMap[word] = 1;
            else
                allWordFreqMap[word]++;

        // Divide keyword into dummies and descendents
        if (ptr->keywords.size() <= offset)
            dummy.push_back(ptr); // this query can't be partitioned by current offset
        else {
            size_t curWord = ptr->keywords[offset];
            if (offsetWordQryMap.find(curWord) == offsetWordQryMap.end()) // new word to offset keyword map
                offsetWordQryMap[curWord] = std::vector<QueryNested *>(); 
            offsetWordQryMap[curWord].push_back(ptr); // count word frequency
        }
    }
    dictSize /= subQry.size(); // get average keyword size in passed queries.

    // Get total frequency of all appearing words
    double totalFreq = 0;
    for (const auto &pair : allWordFreqMap) totalFreq += double(pair.second);
    totalFreq /= dictSize;

    // Convert the word statistic map to vector for random access
    std::vector<std::pair<size_t, std::vector<QueryNested *>>> offsetWordQryVec;
    offsetWordQryVec.insert(offsetWordQryVec.begin(), offsetWordQryMap.begin(), offsetWordQryMap.end());

    // Propose a partition which divides words by a fixed offset
    size_t nPart = std::min(f, offsetWordQryVec.size()); // number of words to be partitioned may be less than nPart
    size_t nWordPerPart = bucketCapacity(offsetWordQryVec.size(), nPart);
    std::vector<KeywordCut> propCuts; // proposed cuts
    propCuts.reserve(nPart);
    for (size_t i = 0; i < nPart; i++) {
        size_t cutStart = i * nWordPerPart; // cutStart and cutEnd are indexes in offsetWordQryVec, not actual word indexes
        if (cutStart >= offsetWordQryVec.size()) { 
            nPart = i;
            break; // break in case of empty buckets
        }
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
    partition.nPart = nPart;
    partition.cost = 0;
    partition.cuts = std::make_unique<KeywordCut[]>(nPart);
    memcpy(partition.cuts.get(), propCuts.data(), nPart * sizeof(nPart));
    partition.dummy = std::move(dummy); 
    partition.queries = std::make_unique<std::vector<QueryNested *>[]>(nPart);
    for (size_t i = 0; i < nPart; i++) {
        partition.cost += computeCutCost(propCuts[i]);
        partition.queries[i] = std::vector<QueryNested *>();
        auto &cutQryVec = partition.queries[i];
        for (size_t word = propCuts[i].start; word <= propCuts[i].end; word++) { 
            if (offsetWordQryMap.find(word) == offsetWordQryMap.end()) continue;
            auto &wordQryVec = offsetWordQryMap[word];
            cutQryVec.insert(cutQryVec.end(), wordQryVec.begin(), wordQryVec.end());
        }
    }
    return partition;
}

APTree::SpatialPartition APTree::spatialHeuristic(const std::vector<QueryNested *> &subQry, const Boundf &bound)
{
    SpatialPartition partition;
    size_t nPart = std::min(f, subQry.size());
    size_t nPartX = partition.nPartX = std::floor(std::sqrt(nPart));
    size_t nPartY = partition.nPartY = nPart / nPartX;

    // Prune dummy queries
    std::vector<QueryNested *> cellQry, dummy;
    for (const auto qry : subQry)
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
        for (size_t i = 0; i < bndStatVec.size(); i++) {
            auto &qryVec = bndQryMapX[bndStatVec[i]] = std::vector<QueryNested *>();
            if (i == bndStatSet.size() - 1) continue;
            Boundf curBnd = {{bndStatVec[i], bound.min.y}, {bndStatVec[i + 1], bound.max.y}};
            for (const auto qry : subQry)
                if (curBnd.Overlaps(qry->region))
                    qryVec.push_back(qry);
        }

        // Convert query map to vector for random access
        std::vector<std::pair<double, std::vector<QueryNested *>>> bndQryVec;
        bndQryVec.insert(bndQryVec.end(), bndQryMapX.begin(), bndQryMapX.end());

        // Propose a partition which divides queries by a fixed bound number offset
        size_t nBndPerPart = bucketCapacity(bndStatVec.size(), nPartX);
        std::unique_ptr<size_t[]> partIdx(new size_t[nPartX + 1]); // indexes in bound statistics vector
        for (size_t i = 0; i < nPartX; i++) {
            if (i * nBndPerPart >= bndStatVec.size()) {
                nPartX = i;
                break;
            }
            partIdx[i] = i * nBndPerPart;
        }
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
        partition.partX = std::make_unique<double[]>(nPartX + 1);
        for (size_t i = 0; i <= nPartX; i++)
            partition.partX[i] = bndStatVec[partIdx[i]];
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
        for (size_t i = 0; i < bndStatVec.size(); i++) {
            auto &qryVec = bndQryMapY[bndStatVec[i]] = std::vector<QueryNested *>();
            if (i == bndStatSet.size() - 1) continue;
            Boundf curBnd = {{bound.min.x, bndStatVec[i]}, {bound.max.x, bndStatVec[i + 1]}};
            for (const auto qry : subQry)
                if (curBnd.Overlaps(qry->region))
                    qryVec.push_back(qry);
        }

        // Convert query map to vector for random access
        std::vector<std::pair<double, std::vector<QueryNested *>>> bndQryVec;
        bndQryVec.insert(bndQryVec.end(), bndQryMapY.begin(), bndQryMapY.end());

        // Propose a partition which divides queries by a fixed bound number offset
        size_t nBndPerPart = bucketCapacity(bndStatVec.size(), nPartY);
        std::unique_ptr<size_t[]> partIdx(new size_t[nPartY + 1]);
        for (size_t i = 0; i < nPartY; i++) {
            if (i * nBndPerPart >= bndStatVec.size()) {
                nPartY = i;
                break;
            }
            partIdx[i] = i * nBndPerPart;
        }
        partIdx[nPartY] = bndStatVec.size() - 1;

        // Define cost function for each proposed axis partition
        auto computeCost = [&] (size_t thisIndex, size_t nextIndex) {
            size_t weight = 0;
            for (size_t i = thisIndex; i < nextIndex; i++)
                weight += bndQryVec[i].second.size();
            double length = bndStatVec[nextIndex] - bndStatVec[thisIndex];
            return weight * length / (bound.max.y - bound.min.y);
        };

        // Find a partition which evenly partitions Y axis by weight through heuristic algorithm
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
        partition.partY = std::make_unique<double[]>(nPartY + 1);
        for (size_t i = 0; i <= nPartY; i++)
            partition.partY[i] = bndStatVec[partIdx[i]];
    }

    // Generate cells and compute final cost
    partition.cost = 0;
    partition.cells = std::unique_ptr<std::vector<QueryNested *>[]>(new std::vector<QueryNested *>[nPartX * nPartY]);
    for (size_t i = 0; i < nPartX; i++) {
        for (size_t j = 0; j < nPartY; j++) {
            size_t index = j + i * nPartY;
            // Pick up queries between bounds
            double thisBndX = partition.partX[i], nextBndX = partition.partX[i + 1];
            auto bndIteX = bndQryMapX.find(thisBndX);
            std::set<QueryNested *> qrySetX;
            while (bndIteX->first != nextBndX) {
                qrySetX.insert(bndIteX->second.begin(), bndIteX->second.end());
                bndIteX++;
            }
            double thisBndY = partition.partY[j], nextBndY = partition.partY[j + 1];
            auto bndIteY = bndQryMapY.find(thisBndY);
            std::set<QueryNested *> qrySetY;
            while (bndIteY->first != nextBndY) {
                qrySetY.insert(bndIteY->second.begin(), bndIteY->second.end());
                bndIteY++;
            }
            // Find common queries of two axes
            auto qryCell = commonElements(qrySetX, qrySetY);
            partition.cells[index] = std::vector<QueryNested *>();
            if (qryCell.size() > 0) {
                size_t weight = qryCell.size();
                double prob = (partition.partX[i + 1] - partition.partX[i]) * 
                              (partition.partY[i + 1] - partition.partY[i]);
                prob /= bound.Area();
                partition.cost += weight * prob;
                partition.cells[index] = std::move(qryCell);
            }
        }
    }

    return partition;
}

std::vector<Query> APTree::Match(const STObject &obj) const {
    // Convert input spatial-textual object to nested one
    STObjectNested stObjN{obj.location, std::vector<size_t>()};
    for (const auto &kw : obj.keywords)
        stObjN.keywords.push_back(dictIndex.find(kw)->second);
    std::sort(stObjN.keywords.begin(), stObjN.keywords.end());

    // Start match recursion
    std::set<QueryNested> queries;
    match(stObjN, 0, root, queries);

    // Convert nested queries to output struct
    std::vector<Query> output;
    for (const auto &qry : queries) {
        Query oQry{qry.region, std::set<std::string>()};
        for (const size_t i : qry.keywords) 
            oQry.keywords.insert(dict[i]);
        output.push_back(std::move(oQry));
    }
    return output;
}

void APTree::match(const STObjectNested &obj, size_t offset, const Node *node, std::set<QueryNested> &out) const
{
    if (node->type == Node::QUERY) {
        // Verify queries in query storage and insert the matched ones to output set
        for (const auto &qry : node->query->queries)
            if (qry.region.Contains(obj.location) && isSubset(obj.keywords, qry.keywords))
                out.insert(qry);

    } else if (node->type == Node::KEYWORD) {
        // Record if certain cut has been visited before
        auto nPart = node->keyword->nPart;
        auto visited = std::make_unique<bool[]>(nPart);
        for (size_t i = 0; i < nPart; i++)
            visited[i] = false;
        // Find the corresponding cut based on ith word in object keywords
        for (size_t i = offset; i < obj.keywords.size(); i++) {
            auto cutIdx = node->keyword->Search(obj.keywords[i]);
            if (cutIdx == INDEX_NOT_FOUND) 
                continue;
            auto cut = node->keyword->cuts[cutIdx];
            if (!visited[cutIdx]) {
                visited[cutIdx] = true;
                match(obj, i + 1, node->keyword->children[cutIdx].get(), out);
            }
        }
        // Search in dummy cut if it exists
        if (node->dummy.get()) 
            match(obj, offset, node->dummy.get(), out);

    } else if (node->type == Node::SPATIAL) {
        // Find the cell which covers the location of object
        auto cellIdx = node->spatial->GetCellIndex(obj.location);
        if (cellIdx.x != INDEX_NOT_FOUND && cellIdx.y != INDEX_NOT_FOUND) { // index is valid
            size_t vecIdx = cellIdx.y + cellIdx.x * node->spatial->nPartY;
            match(obj, offset, node->spatial->cells[vecIdx].get(), out);
        }
        // Search in dummy cell if it exists
        if (node->dummy.get())
            match(obj, offset, node->dummy.get(), out);
    }
}

/*
    Index Maintenance Strategy
    Since index maintenance is just briefly introduced in original paper, it should be elaborated on for practical implementation

    Algorithm:
    # Return new node if reconstructed, or the original one if not
    register(node, newQry): 
        if node is q-node:
            if node.Q.size() + newQry.size() < theta_q:
                node.Q = Merge node.Q and newQry
                return node
            else:
                subQry := Merge node.Q and newQry
                build(newNode, subQry)
                return newNode
        else: # the same for k-node and s-node
            Get statistics of queries
            D_KL := Compute KL-Divergence of new node
            if D_KL > theta_KL:
                subQry := Merge node.Q and newQry
                build(newNode, subQry)
                return newNode
            else:
                dmyQry := Find dummy queries in newQry
                    if dmyQry is not empty:
                        if node has dummy:
                            newDummy := register(dummy, dmyQry)
                            if newDummy != dummy:
                                Reset dummy node
                        else:
                            build(dummy, dmyQry)
                for each bucket:
                    bcktNewQry := Find matched queries in newQry
                    newBucket := register(bucket, bcktNewQry)
                    if newBucket != bucket:
                        Reset bucket node
                return node

    Modification:
    1. Change prototype of private build function to void APTree::build(Node *node, const std::vector<QueryNested *> &subQry)
    size_t offset, bool useKw, bool useSp become members in Node
    2. Add members recording weights(numbers of queries related) in Node for KL-Divergence computation

*/

void APTree::Register(const std::vector<Query> &newQry) 
{
    // Convert queries to nested type
    std::vector<QueryNested> nstdQry;
    for (const auto q : newQry) {
        QueryNested nq{ q.region, std::vector<size_t>() };
        for (const auto kw : q.keywords)
            nq.keywords.push_back(dictIndex[kw]);
        nstdQry.push_back(std::move(nq));
    }

    // Remove repeated queries
    std::sort(nstdQry.begin(), nstdQry.end());
    auto newEnd = std::unique(nstdQry.begin(), nstdQry.end());
    nstdQry.erase(newEnd, nstdQry.end());

    // Get nested query pointers
    std::vector<QueryNested *> nstdQryPtr;
    for (auto &q : nstdQry)
        nstdQryPtr.push_back(&q);

    // Begin registration recursion
    auto newRoot = regist(root, nstdQryPtr); // the whole tree may be reconstructed
    if (newRoot != root) {
        delete root;
        root = newRoot;
    }
    
    // COUT(collectAndMerge(root, std::vector<QueryNested *>()).size())
}

APTree::Node * APTree::regist(Node *node, const std::vector<QueryNested *> &newQry)
{
    if (node->type == Node::QUERY) {
        auto &query = node->query->queries;
        size_t newSize = query.size() + newQry.size();

        if (newSize < theta_Q || (!node->useKw && !node->useSp)) { // keep q-node structure
            // Append new queries to q-node
            for (const auto q : newQry)
                query.push_back(*q);
            // Remove repeated queries
            std::sort(query.begin(), query.end());
            auto newEnd = std::unique(query.begin(), query.end());
            query.erase(newEnd, query.end());
            return node; // not reconstructed
        }
        else { // should reconstruct
            auto recQry = collectAndMerge(node, newQry);
            auto newNode = new Node(*node);
            build(newNode, recQry);
            return newNode;
        }
    }
    else if (node->type == Node::KEYWORD) {
        // Get statistics of newly added queries
        std::vector<QueryNested *> dmyQry;
        auto nPart = node->keyword->nPart;
        auto cutQryStat = std::make_unique<std::vector<QueryNested *>[]>(nPart);
        for (size_t i = 0; i < nPart; i++) // initialize stat map
            cutQryStat[i] = std::vector<QueryNested *>();
        for (const auto q : newQry) {
            if (q->keywords.size() <= node->offset) // is dummy query
                dmyQry.push_back(q);
            else {
                auto cutIdx = node->keyword->Search(q->keywords[node->offset]);
                if (cutIdx != INDEX_NOT_FOUND)  // cut is valid
                    cutQryStat[cutIdx].push_back(q);
                else { // current node should be rebuilt because query falls in gap
                    auto recQry = collectAndMerge(node, newQry);
                    auto newNode = new Node(*node);
                    build(newNode, recQry);
                    return newNode;
                }
            }
        }

        // Compute KL-Divergence D_KL(w_old|w) of current node
        // See https://en.wikipedia.org/wiki/Kullback-Leibler_divergence for definition of KL-Divergence
        // D_KL(P|Q) = \sum_{x\in X} P(x)*log(P(x)/Q(x))
        // In terms of AP-Tree, D_KL(w_old|w) = \sum_{i=1}^f w_old(B_i)*log(w_old(B_i)/w(B_i))
        double D_KL = 0;
        for (size_t i = 0; i < node->keyword->nPart; i++) {
            node->keyword->nAdd[i] += cutQryStat[i].size(); // update number of newly added queries
            double wOld = node->keyword->nOld[i], wAdd = node->keyword->nAdd[i];
            double bKL = wOld * std::log10(1.0 + wAdd / wOld);
            D_KL += bKL;
        }

        if (D_KL > theta_KL) {
            // Reconstruct node if D_KL exceeds threshold theta_KL
            auto recQry = collectAndMerge(node, newQry);
            auto newNode = new Node(*node);
            build(newNode, recQry);
            return newNode;
        }
        else { 
            // Keep current node but modify it
            // Build dummy node if possible
            if (dmyQry.size() > 0) {
                if (node->dummy.get()) { // previous node has dummy cut
                    auto newDmyNode = regist(node->dummy.get(), dmyQry);
                    if (newDmyNode != node->dummy.get())
                        node->dummy.reset(newDmyNode);
;               }
                else {
                    node->dummy = std::make_unique<Node>(node->offset + 1, node->bound, false, node->useSp);
                    build(node->dummy.get(), dmyQry);
                }
            }

            // Distribute queries to corresponding cuts
            for (size_t i = 0; i < node->keyword->nPart; i++) {
                auto cutNewQry = cutQryStat[i];
                auto oldNode = node->keyword->children[i].get();
                auto newNode = regist(oldNode, cutNewQry);
                if (newNode != oldNode) 
                    node->keyword->children[i].reset(newNode);
            }

            return node;
        }
    }
    else {
        // Get statistics of newly added queries
        size_t nPartX = node->spatial->nPartX, nPartY = node->spatial->nPartY;
        std::vector<QueryNested *> dmyQry;
        auto cellQryStat = std::make_unique<std::vector<QueryNested *>[]>(nPartX * nPartY);
        for (size_t idx = 0; idx < nPartX * nPartY; idx++)
            cellQryStat[idx] = std::vector<QueryNested *>();
        for (const auto qry : newQry) {
            // Pick out dummy query
            if (qry->region.Contains(node->bound))
                dmyQry.push_back(qry);
            // Add query to overlapping cell
            for (size_t idx = 0; idx < nPartX * nPartY; idx++)
                if (qry->region.Overlaps(node->spatial->cells[idx]->bound))
                    cellQryStat[idx].push_back(qry);
        }

        // Compute KL-Divergence of current node
        double D_KL = 0;
        for (size_t idx = 0; idx < nPartX * nPartY; idx++) {
            node->spatial->nAdd[idx] += cellQryStat[idx].size();
            double wOld = node->spatial->nOld[idx], wAdd = node->spatial->nAdd[idx];
            if (wOld == 0) continue; // account for empty cells
            double bKL = wOld * std::log10(1.0 + wAdd / wOld);
            D_KL += bKL;
        }

        if (D_KL > theta_KL) {
            // Reconstruct current node
            auto recQry = collectAndMerge(node, newQry);
            auto newNode = new Node(*node);
            build(newNode, recQry);
            return newNode;
        }
        else {
            // Keep current node
            // Build dummy node if possible
            if (dmyQry.size() > 0) {
                if (node->dummy.get()) {
                    auto newDmyNode = regist(node->dummy.get(), dmyQry);
                    if (newDmyNode != node->dummy.get())
                        node->dummy.reset(newDmyNode);
                }
                else {
                    node->dummy = std::make_unique<Node>(node->offset, node->bound, node->useKw, false);
                    build(node->dummy.get(), dmyQry);
                }
            }

            // Distribute queries to corresponding cells
            for (size_t idx = 0; idx < nPartX * nPartY; idx++) {
                auto cellNewQry = cellQryStat[idx];
                auto newNode = regist(node->spatial->cells[idx].get(), cellNewQry);
                if (newNode != node->spatial->cells[idx].get())
                    node->spatial->cells[idx].reset(newNode);
            }

            return node;
        }
    }
}

// Only collect query pointers for selected node, and merge with newly added ones.
// The actual value is fetched only when q-node is constructed.
std::vector<APTree::QueryNested *> APTree::collectAndMerge(const Node *node, const std::vector<QueryNested *> &newQry) const
{
    // Begin collecting recursion
    std::vector<QueryNested *> out;
    collect(node, out);

    // Merge newly added queries
    out.insert(out.end(), newQry.begin(), newQry.end());
    
    // Remove repeated elements
    std::sort(out.begin(), out.end(), [](auto p1, auto p2) { return *p1 < *p2; });
    auto newEnd = std::unique(out.begin(), out.end(), [](auto p1, auto p2) { return *p1 == *p2; });
    out.erase(newEnd, out.end());

    return out;
}

void APTree::collect(const Node *node, std::vector<QueryNested *> &out) const
{
    if (node->type == Node::QUERY) {
        auto &queries = node->query->queries;
        for (auto &q : queries)
            out.push_back(&q);
    }
    else if (node->type == Node::KEYWORD) {
        if (node->dummy.get())
            collect(node->dummy.get(), out);
        for (size_t i = 0; i < node->keyword->nPart; i++)
            if (node->keyword->children[i].get())
                collect(node->keyword->children[i].get(), out);
    }
    else if (node->type == Node::SPATIAL) {
        if (node->dummy.get())
            collect(node->dummy.get(), out);
        for (size_t idx = 0; idx < node->spatial->nPartX * node->spatial->nPartY; idx++)
            collect(node->spatial->cells[idx].get(), out);
    }
}