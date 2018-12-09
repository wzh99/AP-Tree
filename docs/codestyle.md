# Code Style

## 4 Spaces as Tab

## Don't wrap braces

## Naming conventions
* Upper camel case for class names
* Lower camel case for function and variable names
* All upper case for macros and enumeration entries

## No `using namespace`

## Use `#pragma once` instead of `#ifndef`

## Example
```c++
#pragma once

#include <memory>

#define PI 3.14159

enum Fruit {
    APPLE,
    BANANA,
    ORANGE
};

template <class Type>
class Stack {
public:
    Stack();
    void push(const Type &element);
    Type pop();
    Type top() const {
        return data[size - 1];
    }
private:
    std::unique_ptr<Type> data;
    size_t size;
};
```