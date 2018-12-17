#pragma once

#include <algorithm>
#include <cinttypes>

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
    PointType operator + (const PointType &pt) const { return {x + pt.x, y + pt.y}; }
    PointType operator - (const PointType &pt) const { return (x - pt.x, y - pt.y); }
    PointType operator * (Type c) const { return {c * x, c * y}; }

    bool operator == (const PointType &pt) const { return x == pt.x && y == pt.y; }
    bool operator != (const PointType &pt) const { return !(pt == *this); }
};

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

    bool operator == (const BoundType &bd) const { return min == bd.min && max == bd.max; }
    bool operator != (const BoundType &bd) const { return !(bd == *this); }

    bool Contains(const PointType &pt) const {
        return pt.x > min.x && pt.y > min.y && pt.x < max.x && pt.y < max.y;
    }

    bool Contains(const BoundType &bd) const {
        return bd.min.x > min.x && bd.min.y > min.y && bd.max.x < max.x && bd.max.y < max.y;
    }

    bool Overlaps(const BoundType &bd) const {
        return !(bd.min.x >= max.x || bd.min.y >= max.y || bd.max.x <= min.x || bd.max.y <= min.y);
    }

    Type Area() const { return (max.x - min.x) * (max.y - min.y); }
};

using Boundi = Bound<int>;
using Boundu = Bound<size_t>;
using Boundf = Bound<double>;