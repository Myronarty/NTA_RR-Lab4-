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
    for (int t = 0; t < n; t++)
    {
        b[k * n + t] -= q * b[j * n + t];
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

bool chek_lovas(const vector<double>& b_, double mu, int i, double delta)
{
    return b_[i] >= (delta - (mu * mu)) * b_[i - 1];
}

bool chek_min(double mu)
{
    return abs(mu) > 0.5;
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

void LLL(LatticeBasis& B, int n, double delta)
{
    vector<vector<double>> mu(n, vector<double>(n, 0.0));
    vector<double> b_(n, 0.0);

    GSh(B, n, mu, b_);

    int k = 1;
    while (k < n)
    {
        for (int j = k - 1; j > -1; j--)
        {
            if (chek_min(mu[k][j]))
            {
                int64_t q = std::round(mu[k][j]);
                B.reduce_vector(k, j, q);

                GSh(B, n, mu, b_);
            }
        }

        if (chek_lovas(b_, mu[k][k - 1], k, delta))
        {
            k++;
        }
        else
        {
            B.swap_rows(k, k - 1);
            GSh(B, n, mu, b_);
            k = max(1, k - 1);
        }
    }
}

double trikutnik(const LatticeBasis& B, int n)
{
    vector<vector<double>> mu(n, vector<double>(n, 0.0));
    vector<double> b_(n, 0.0);

    GSh(B, n, mu, b_);

    double det = 1.0;
    for (int i = 0; i < n; ++i)
    {
        det *= b_[i];
    }

    return det;
}

double big_prod(const LatticeBasis& B, int n)
{
    double rez = 1;
    for (int i = 0; i < n; i++)
    {
        rez *= sqrt(B.dot_product(i, i));
    }
    return rez;
}

double Adamar(const LatticeBasis& B, int n)
{
    double det = trikutnik(B, n);

    double ob = std::sqrt(det);

    double bb = big_prod(B, n);

    double H = ob / bb;

    return pow(H, 1.0 / n);
}
