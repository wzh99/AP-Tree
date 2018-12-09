#pragma once

#include "DynamicalSearchTable.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <vector>

using uint = uint32_t;
using ulong = uint64_t;

enum HashFunction {
    REMAINDER, // Divide key by size and get its remainder
    DIGIT_ANALYSIS, // Select digits that vary most greatly
    SQUARE_MEDIAN, // Take square of key and select digits in the middle
    FOLD, // Add digits by group
    NUM_HASH_FUNCTIONS
};

enum CollisionSolution {
    LINEAR_DETECTION,
    SQUARE_DETECTION,
    NUM_COLLISION_SOLUTION
};

template <class KeyType, class ValueType>
class ClosedHashTable : public DynamicSearchTable<KeyType, ValueType> {
public:
    using PairType = Pair<KeyType, ValueType>;
    
    ClosedHashTable(const PairType *data, uint size, uint nKeyDigits, HashFunction hf, CollisionSolution cs)
        : pSize(nextPrime(2 * size)),
          // doubled size ensures insertion success
          // for efficiency, table size should be next prime of data size
          cs(cs) 
    {
        // Allocate table
        table = std::unique_ptr<Node[]>(new Node[pSize]);

        // Construct hash function objects
        switch (hf) {
        case REMAINDER: 
            hash = [] (KeyType key) { return ulong(key); }; break;
        case DIGIT_ANALYSIS:
            hash = DigitAnalysis(data, size, nKeyDigits); break;
        case SQUARE_MEDIAN:
            hash = SquareMedian(data, size, nKeyDigits); break;
        case FOLD:
            hash = Fold(data, size, nKeyDigits); break;
        default: 
            break;
        }
    }

    PairType * find(const KeyType &key) override {
        int pos = 0;
        switch (cs) {
        case LINEAR_DETECTION:
            pos = linearDetect(key, FIND); break;
        case SQUARE_DETECTION:
            pos = squareDetect(key, FIND); break;
        default: 
            break;
        }

        if (pos >= 0) 
            return &table[pos].data;
        return nullptr;
    }

    void insert(const PairType &pair) override {
        int pos = -1;
        switch (cs) {
        case LINEAR_DETECTION:
            pos = linearDetect(pair.key, INSERT); break;
        case SQUARE_DETECTION:
            pos = squareDetect(pair.key, INSERT); break;
        default:
            break;
        }

        if (pos >= 0) {
            table[pos].data = pair;
            table[pos].state = Node::ACTIVE;
        } else 
            std::cout << "Insertion failed. \n";
    }

    void remove(const KeyType &key) override {
        int pos = -1;
        switch (cs) {
        case LINEAR_DETECTION:
            pos = linearDetect(key, REMOVE); break;
        case SQUARE_DETECTION:
            pos = squareDetect(key, REMOVE); break;
        default:
            break;
        }

        if (pos >= 0)
            table[pos].state = Node::REMOVED;
    }

    static inline uint getDigit(KeyType key, uint index) {
        while (index-- > 0) { key /= 10; } 
        return key % 10;
    }

    static inline uint countDigits(KeyType key) { return std::log10(double(key)) + 1; }

    static inline uint nextPrime(uint n) {
        if (n % 2 == 0) n++;
        for (; ; n += 2) {
            for (uint i = 3; i * i <= n; i += 2)
                if (n % i == 0)
                    goto OuterLoop;
                return n;
            OuterLoop: ;
        }
    }

private:
    struct Node {
        PairType data;
        enum State {
            EMPTY,
            ACTIVE,
            REMOVED
        } state;
        Node() : state(EMPTY) {}
    };
    
    std::unique_ptr<Node[]> table;
    uint pSize; // table size in prime
    std::function<ulong(KeyType)> hash;
    CollisionSolution cs;

    // Hash function class definitions
    class DigitAnalysis {
    public:
        DigitAnalysis(const PairType *data, uint size, uint nKeyDigits) {
            uint nSelectedDigits = countDigits(nextPrime(2 * size));

            // Analyze digit variation
            std::vector<std::set<uint>> digitVariation(nKeyDigits);
            for (uint i = 0; i < size; i++) {
                auto key = data[i].key;
                for (uint j = 0; j < nKeyDigits; j++) {
                    uint digit = key % 10;
                    digitVariation[j].insert(digit);
                    key /= 10;
                }
            }

            // Select digits that vary most greatly
            std::vector<std::pair<uint, uint>> digitStats;
            for (uint i = 0; i < nKeyDigits; i++) 
                digitStats.push_back({ i, uint(digitVariation[i].size()) });
            std::sort(digitStats.begin(), digitStats.end(), 
                      [] (const auto &p1, const auto &p2) { return p1.second > p2.second; });

            // Build digit weight map
            for (uint i = 0; i < nSelectedDigits; i++) {
                std::for_each(digitWeights.begin(), digitWeights.end(), [] (auto &p) { p.second *= 10; });
                auto digit = digitStats[i].first;
                digitWeights.insert({digit, 1});
            }
        }

        ulong operator () (KeyType key) const {
            ulong ret = 0;
            for (const auto &p : digitWeights)
                ret += getDigit(key, p.first) * p.second;
            return ret;
        }

    private:
        std::map<uint, uint> digitWeights; // first - index, second - weight
    };

    class SquareMedian {
    public:
        SquareMedian(const PairType *data, uint size, uint nKeyDigits) {
            uint nSelectedDigits = countDigits(nextPrime(2 * size));
            uint begin = (nKeyDigits - nSelectedDigits) / 2;
            for (uint digit = begin; digit < (begin + nSelectedDigits); digit++) {
                std::for_each(digitWeights.begin(), digitWeights.end(), 
                              [] (auto &p) { p.second *= 10; });
                digitWeights.insert({digit, 1});
            }
        }

        ulong operator () (KeyType key) const {
            key = key * key;
            ulong ret = 0;
            for (const auto &p : digitWeights)
                ret += getDigit(key, p.first) * p.second;
            return ret;
        }

    private:
        std::map<uint, uint> digitWeights;
    };

    class Fold {
    public:
        Fold(const PairType *data, uint size, uint nKeyDigits) 
            : step(countDigits(nextPrime(2 * size))), nKeyDigits(nKeyDigits)
        {
            for (uint i = 0; i < step; i++) {
                std::for_each(digitWeights.begin(), digitWeights.end(), 
                              [] (auto &p) { p.second *= 10; });
                digitWeights.insert({i, 1});
            }
        }

        ulong operator () (KeyType key) const {
            ulong ret = 0;
            for (uint i = 0; i < nKeyDigits; i++) 
                ret += getDigit(key, i) * digitWeights.at(i % step);
            return ret;
        }

    private:
        uint step;
        uint nKeyDigits;
        std::map<uint, uint> digitWeights;
    };

    enum Operation {
        FIND,
        INSERT,
        REMOVE
    };

    int linearDetect(const KeyType &key, Operation op) const { // -1 indicates "NOT FOUND"
        int initPos, pos;
        initPos = pos = hash(key) % pSize;
        do {
            switch (op) {
            case FIND: case REMOVE: 
                if (table[pos].state == Node::EMPTY) return -1; 
                if (table[pos].state == Node::ACTIVE && table[pos].data.key == key) 
                    return pos;
                break;

            case INSERT:
                if (table[pos].state != Node::ACTIVE)
                    return pos;
            }

            if (++pos == pSize) pos -= pSize;
        } while (pos != initPos);

        return -1;
    }

    int squareDetect(const KeyType &key, Operation op) const {
        int i = 1, pos = hash(key) % pSize;
        do {
            switch (op) {
            case FIND: case REMOVE: 
                if (table[pos].state == Node::EMPTY) return -1; 
                if (table[pos].state == Node::ACTIVE && table[pos].data.key == key) 
                    return pos;
                break;

            case INSERT:
                if (table[pos].state != Node::ACTIVE)
                    return pos;
            }

            pos += (2 * i - 1); 
            if (pos >= pSize) pos -= pSize;
            i++;
        } while ((2 * i - 1) < pSize);
        
        return -1;
    }
};