#include "LLL.h"

int64_t& LatticeBasis::operator()(int i, int j)
{
    return b[i * n + j];
}

void LatticeBasis::swap_rows(int row1, int row2)
{
    for (int j = 0; j < n; ++j)
    {
        swap(b[row1 * n + j], b[row2 * n + j]);
    }
}

void LatticeBasis::reduce_vector(int k, int j, int64_t q)
{
    for (int k = 0; k < n; ++k)
    {
        b[k * n + k] -= q * b[j * n + k];
    }
}

int64_t LatticeBasis::dot_product(int i, int j) const
{
    int64_t sum = 0;
    for (int k = 0; k < n; ++k)
    {
        sum += b[i * n + k] * b[j * n + k];
    }
    return sum;
}

double LatticeBasis::comput_mu(int i, int j) const
{
    return double(dot_product(i, j)) / dot_product(i, i);
}

int64_t ocruglenie(double mu)
{
    if (mu >= 0.0)
    {
        return int64_t(mu + 0.5);
    }
    else
    {
        return int64_t(mu - 0.5);
    }
}

void GSh(const LatticeBasis& B, int n, vector<vector<double>>& mu, vector<double>& b_)
{
    for (int i = 0; i < n; ++i)
    {
        double norm = static_cast<double>(B.dot_product(i, i));

        for (int j = 0; j < i; ++j)
        {
            double temp = static_cast<double>(B.dot_product(i, j));

            for (int k = 0; k < j; ++k)
            {
                temp -= mu[i][k] * mu[j][k] * b_[k];
            }

            mu[i][j] = temp / b_[j];

            norm -= mu[i][j] * mu[i][j] * b_[j];
        }

        b_[i] = norm;

        mu[i][i] = 1.0;
    }
}

void LLL(LatticeBasis& B, int n, double delta = 0.75)
{
    vector<vector<double>> mu(n, std::vector<double>(n, 0.0));
    vector<double> b_(n, 0.0);
}
