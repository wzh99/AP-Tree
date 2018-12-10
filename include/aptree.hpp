#pragma once

#include "primitive.hpp"
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

struct STObject {
    Pointf location;
    std::set<std::string> keywords; 
};

struct Query {
    Boundf region;
    std::set<std::string> keywords;
};

class APTree {
public:
    APTree(const std::vector<std::string> &vocab, const std::vector<Query> &qry);
    std::set<Query> Match(const STObject &obj);
    void Register(const Query &query);
    void Deregister(const Query &query);

private:
    struct QueryNested;
    struct STObjectNested;
    struct Node;
    struct KeywordCut;

    std::unique_ptr<Node> root; // root node of the whole AP-Tree
    std::vector<std::string> dict; // stores all the vocabulary
    std::unordered_map<std::string, size_t> dictIndex; // stores index of keywords in dictionart
};