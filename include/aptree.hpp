#pragma once

#include "primitive.hpp"
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

struct STObject {
    Pointf loc;
    std::set<std::string> kw; 
};

struct Query {
    Boundf rgn;
    std::set<std::string> kw;
};

class APTree {
public:
    APTree(const std::vector<std::string> &vol, const std::vector<Query> &qry);
    std::set<Query> Match(const STObject &obj);
    void Register(const Query &query);
    void Deregister(const Query &query);

private:
};