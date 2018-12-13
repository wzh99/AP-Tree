#pragma once

#include "primitive.hpp"
#include <set>
#include <string>
#include <vector>

struct STObject {
    Pointf location;
    std::set<std::string> keywords; 
};

struct Query {
    Boundf region;
    std::set<std::string> keywords;
};

class STTree {
public:
    virtual std::vector<Query> Match(const STObject &obj) const = 0;
};