#pragma once

#include "primitive.hpp"
#include <set>
#include <string>
#include <vector>

struct STObject {
    Pointf location;
    std::set<std::string> keywords; 
};

inline std::ostream & operator << (std::ostream &os, const STObject &obj) {
    os << "{ location: " << obj.location << ", keywords: [ ";
    for (const auto &kw : obj.keywords)
        os << kw << ' ';
    os << "]\n";
    return os;
}

struct Query {
    Boundf region;
    std::set<std::string> keywords;
};

inline std::ostream & operator << (std::ostream &os, const Query &qry) {
    os << "{ region: " << qry.region << ", keywords: [ ";
    for (const auto &kw : qry.keywords)
        os << kw << ' ';
    os << "]\n";
    return os;
}

class STTree {
public:
    virtual std::vector<Query> Match(const STObject &obj) const = 0;
};