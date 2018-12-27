#pragma once

#include <chrono>
#include <iostream>

#define TIME_COUNT(code, ite, duration)  { \
    Clock clock; \
    for (int i = 0; i < ite; i++) code; \
    std::cout << clock.count<std::chrono::duration>() / ite << std::endl; }

class Clock {
public:
    Clock() : begin(std::chrono::steady_clock::now()) {}

    template<class duration = std::chrono::microseconds>
    long long count() {
        auto end = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<duration>(end - begin).count();
    }

private:
    const std::chrono::steady_clock::time_point begin;
};
