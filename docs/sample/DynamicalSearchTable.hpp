#pragma once

template <class KeyType, class ValueType>
struct Pair {
    KeyType key;
    ValueType value;
};

template <class KeyType, class ValueType>
class DynamicSearchTable {
public:
    virtual Pair<KeyType, ValueType> * find(const KeyType &key) = 0;
    virtual void insert(const Pair<KeyType, ValueType> &pair) = 0;
    virtual void remove(const KeyType &key) = 0;
    virtual ~DynamicSearchTable() {}
};