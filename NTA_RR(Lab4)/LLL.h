#pragma once
#include <iostream>
#include <vector>

using namespace std;

const int s1 = 30;

class LatticeBasis
{
private:
    int n;
    std::vector<int64_t> b;

public:
    LatticeBasis(int size) : n(size), b(size* size, 0) {}

    int64_t& operator()(int i, int j);

    void swap_rows(int row1, int row2);

    void reduce_vector(int k, int j, int64_t q);

    int64_t dot_product(int i, int j) const;

    double comput_mu(int i, int j) const;
};

void GSh(const LatticeBasis& B, int n, vector<vector<double>>& mu, vector<double>& b_);

void LLL(LatticeBasis& B, int n, double delta = 0.75);
