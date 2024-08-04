#pragma once

#include <iostream>

typedef unsigned int uint;

#define EXTRA_BITS 5

template <uint width>
class CSD {
    public:
    int bits[width + EXTRA_BITS]; // the CSD number
    uint len; // the length of the CSD number
    int num; // the integer value of the CSD number
    CSD(); // all zero initialization
    CSD(int x); // initialization from integer
    CSD& operator=(const CSD& other); // assignment operator
    bool operator==(const CSD& other); // equal operator
    bool operator!=(const CSD& other); // not equal operator
    CSD operator -(); // negation operator
    CSD get_sub(int l, int r); // get the sub word from [l, r)

    friend std::ostream& operator<<(std::ostream& os, const CSD& num) {
        for (int i = num.len - 1; i >= 0; i--)
            os << num.bits[i];
        return os;
    }
};

template <uint width>
CSD<width>::CSD() {
    for (int i = 0; i < width + EXTRA_BITS; i++) {
        bits[i] = 0;
    }
    len = 0;
    num = 0;
}

template <uint width>
CSD<width>::CSD(int x) {
    memset(bits, 0, sizeof(bits));
    num = x;
    len = 0;
    bool neg = false;
    if (x < 0) {
        x = -x;
        neg = true;
    }
    while (x) {
        bits[len++] = x & 1;
        x >>= 1;
    }
    for (int i = 0; i < (int)len - 1; i++) {
        if (bits[i] == 1 && bits[i + 1] == 1) {
            bits[i] = -1;
            bits[i + 1] = 0;
            bits[i + 2] += 1;
            if (i == len - 2)
                len++;
        }
        uint k = i + 2;
        if (k >= len) continue;
        while (bits[k] == 2) {
            bits[k] = 0;
            bits[k + 1] += 1;
            k++;
        }
        len = std::max(len, k + 1);
    }
    if (neg) {
        for (int i = 0; i < len; i++) {
            bits[i] = -bits[i];
        }
    }
}

template <uint width>
CSD<width>& CSD<width>::operator=(const CSD& other) {
    if (this == &other) {
        return *this;
    }
    len = other.len;
    num = other.num;
    for (int i = 0; i < len; i++) {
        bits[i] = other.bits[i];
    }
    return *this;
}

template <uint width>
bool CSD<width>::operator==(const CSD& other) {
    if (len != other.len) return false;
    for (int i = 0; i < len; i++) {
        if (bits[i] != other.bits[i]) return false;
    }
    return true;
}

template <uint width>
bool CSD<width>::operator!=(const CSD& other) {
    return !(*this == other);
}

template <uint width>
CSD<width> CSD<width>::operator-() {
    CSD<width> res;
    for (int i = 0; i < len; i++) {
        res.bits[i] = -bits[i];
    }
    res.len = len;
    res.num = -num;
    return res;
}

template <uint width>
CSD<width> CSD<width>::get_sub(int l, int r) {
    CSD<width> res;
    for (int i = l; i < r; i++) {
        res.bits[i - l] = bits[i];
    }
    res.len = r - l;
    res.num = 0;
    for (int i = res.len - 1; i >= 0; i--) {
        res.num = res.num * 2 + res.bits[i];
    }
    return res;
}