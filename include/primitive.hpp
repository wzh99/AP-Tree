#pragma once

#include <algorithm>
#include <cinttypes>
#include <iostream>

using uint = uint32_t;
using ulong = uint64_t;

template <class Type>
struct Point {
    using PointType = Point<Type>;
    Type x = 0, y = 0;

    Point() {}
    template <class OtherType>
    Point(OtherType x, OtherType y) : x(x), y(y) {}
    template <class OtherType>
    Point(const Point<OtherType> &pt) : x(pt.x), y(pt.y) {}

    PointType & operator = (const PointType &pt) { x = pt.x; y = pt.y; return *this; }

    bool operator == (const PointType &pt) const noexcept { return x == pt.x && y == pt.y; }
    bool operator != (const PointType &pt) const noexcept { return !(pt == *this); }
};

template <class Type>
std::ostream & operator << (std::ostream &os, const Point<Type> &pt) {
    os << "[ " << pt.x << ' ' << pt.y << " ]";
    return os;
}

using Pointi = Point<int>;
using Pointu = Point<size_t>;
using Pointf = Point<double>;

template <class Type>
struct Bound {
    using PointType = Point<Type>;
    using BoundType = Bound<Type>;
    PointType min, max;

    Bound() : min(), max() {}
    Bound(const PointType &p1, const PointType &p2) : min(p1), max(p2) {
        if (min.x > max.x) std::swap(min.x, max.x);
        if (min.y > max.y) std::swap(min.y, max.y);
    }
    Bound(const BoundType &bd) : min(bd.min), max(bd.max) {}

    bool operator == (const BoundType &bd) const noexcept { return min == bd.min && max == bd.max; }
    bool operator != (const BoundType &bd) const noexcept { return !(bd == *this); }

    bool Contains(const PointType &pt) const noexcept {
        return pt.x > min.x && pt.y > min.y && pt.x < max.x && pt.y < max.y;
    }

    bool Contains(const BoundType &bd) const noexcept {
        return bd.min.x > min.x && bd.min.y > min.y && bd.max.x < max.x && bd.max.y < max.y;
    }

    bool Overlaps(const BoundType &bd) const noexcept {
        return !(bd.min.x >= max.x || bd.min.y >= max.y || bd.max.x <= min.x || bd.max.y <= min.y);
    }

    Type Area() const { return (max.x - min.x) * (max.y - min.y); }
};

template <class Type>
std::ostream & operator << (std::ostream &os, const Bound<Type> &bd) {
    os << "[ " << bd.min << ' ' << bd.max << " ]";
    return os;
}

using Boundi = Bound<int>;
using Boundu = Bound<size_t>;
using Boundf = Bound<double>;
