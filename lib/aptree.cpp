#include "../include/aptree.hpp"

struct APTree::QueryNested {
    Pointf location;
    std::set<size_t> wordIndexes; // stores indexes of keywords in dictionart for algorithm efficiency
};

struct APTree::STObjectNested {
    Boundf region;
    std::set<size_t> wordIndexes; // the same as QueryNested
};

struct APTree::KeywordCut {
    size_t start, end; // both included
    KeywordCut(size_t start, size_t end) : start(start), end(end) {}
    bool operator < (const KeywordCut &other) const { return start < other.start; }
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
            std::unique_ptr<Node> quad[2][2]; // [x > divide.x][y > divide.y]
            std::unique_ptr<std::vector<QueryNested>> dummy = nullptr; // can be null
        };
    };
};