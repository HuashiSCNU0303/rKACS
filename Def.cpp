#include "Def.h"

int gcd(int num, int deno)
{
    int a = MAX(num, deno), b = MIN(num, deno), c = 1;
    while (b != 0)
    {
        c = a;
        a = b;
        b = c % b;
    }
    return a;
}

Edge GetEdge(int u, int v)
{
    if (u < v) return make_pair(u, v);
    return make_pair(v, u);
}

bool SortEdge(const Edge& e1, const Edge& e2)
{
    if (e1.first != e2.first)
    {
        return e1.first < e2.first;
    }
    return e1.second < e2.second;
}

bool SortEdge2(const Edge& e1, const Edge& e2)
{
    return e1.first < e2.first;
}

bool CompareClique(const Clique& c1, const Clique& c2)
{
    int c1Size = c1.size(), c2Size = c2.size();
    int smallSize = MIN(c1Size, c2Size);

    for (int i = 0; i < smallSize; i++)
    {
        if (c1[i] != c2[i])
        {
            return c1[i] < c2[i];
        }
    }
    return smallSize == c1Size;
}

bool CompareLongClique(const vector<long>& a, const vector<long>& b)
{
    return a.size() > b.size();
}

uint64_t MurmurHash2(const vector<int>& data)
{
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;

    int size = data.size(), len = size * 4;

    uint64_t h = 0 ^ (len * m);

    int i = 0;
    for (; i < size - 1; i += 2)
    {
        int_combine newInt(data[i], data[i + 1]);
        uint64_t k = *((uint64_t*)(&newInt));

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    if (i == size - 1)
    {
        int num = data[size - 1];
        unsigned char* c = (unsigned char*)&num;
        h ^= uint64_t(c[3]) << 24;
        h ^= uint64_t(c[2]) << 16;
        h ^= uint64_t(c[1]) << 8;
        h ^= uint64_t(c[0]);
        h *= m;
    }

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}