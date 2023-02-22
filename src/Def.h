#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iterator>
#include <string>
#include <queue>
#include <functional>
#include <fstream>
#include <tuple>
#include <random>
#include <chrono>
#include <phmap.h>

#if defined(_WIN64)
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")

#elif defined(__linux__)
#include <unistd.h>
#include <sys/resource.h>

#endif

#include "Snap.h"

using namespace std;
using namespace phmap;

constexpr auto GRAPH_FILE = "./DataSets/{0}/graph.txt";
constexpr auto TRUSSNESS_FILE = "./DataSets/{0}/edge_trussness.txt";
constexpr auto ATTRIBUTE_FILE = "./DataSets/{0}/user_attributes_int.txt";
constexpr auto SIMILARITY_FILE = "./DataSets/{0}/edge_similarity.bin";
constexpr auto QUERY_FILE = "./DataSets/{0}/{1}_query.txt";
constexpr auto AVG_RESULT_FILE = "./Results/{0}/{1}/{2}/avg_result.txt";
constexpr auto RESULT_FILE = "./Results/{0}/{1}/{2}/{2}_{3}_{4}.txt";
constexpr auto QUALITY_RESULT_FILE = "./Results/{0}/quality/{1}_quality_result.txt";

constexpr auto TIME_LIMIT = 3600;

constexpr auto BASIC = 1;
constexpr auto KRCORE = 2;
constexpr auto MCMEI = 3;
constexpr auto IMCE = 4;
constexpr auto NIEMCH = 5;
constexpr auto NIEMC = 6;
constexpr auto INCREMENTAL = 7;
const string METHOD_NAMES[7] = { "Basic", "KRCore", "MCMEI", "IMCE", "NIEMCH", "NIEMC", "Incremental" };

// compute the hash value for pair<int, int>
template <typename T>
inline void hash_combine(std::size_t& seed, const T& val) 
{
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
template <typename T> inline void hash_val(std::size_t& seed, const T& val) 
{
    hash_combine(seed, val);
}
template <typename T, typename... Types>
inline void hash_val(std::size_t& seed, const T& val, const Types&... args) 
{
    hash_combine(seed, val);
    hash_val(seed, args...);
}
template <typename... Types>
inline std::size_t hash_val(const Types&... args) 
{
    std::size_t seed = 0;
    hash_val(seed, args...);
    return seed;
}
struct pair_hash 
{
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const 
    {
        return hash_val(p.first, p.second);
    }
};

struct int_combine
{
    int num1, num2;
    int_combine(int n1, int n2)
    {
        num1 = n1;
        num2 = n2;
    }
};
uint64_t MurmurHash2(const vector<int>& data);
struct clique_hash
{
    uint64_t operator()(const vector<int>& c) const
    {
        return MurmurHash2(c); 
    }
};

struct clique_eq
{
    bool operator ()(const vector<int>& lhs, const vector<int>& rhs) const
    {
        return lhs == rhs;
    }
};

// here we use fractions to denote score(u, v, Q_W) and CSQ
int gcd(int num, int deno);
typedef struct Weight 
{
    int numerator;
    int denominator;
    Weight()
    {

    }
    Weight(int num, int deno)
    {
        int c = gcd(num, deno);
        numerator = num / c;
        denominator = deno / c;
    }
    bool operator < (const Weight& b) const 
    {
        return numerator * b.denominator - denominator * b.numerator < 0;
    }
    bool operator > (const Weight& b) const
    {
        return numerator * b.denominator - denominator * b.numerator > 0;
    }
    bool operator == (const Weight& b) const
    {
        return numerator == b.numerator && denominator == b.denominator;
    }
    Weight operator * (const Weight& b) const
    {
        return Weight(numerator * b.numerator, denominator * b.denominator);
    }
    Weight operator * (const int& b) const
    {
        return Weight(numerator * b, denominator);
    }
    Weight Dist() const // Çó1 - a/b
    {
        return Weight(denominator - numerator, denominator);
    }
};
const Weight WEIGHT_ZERO(0, 1);
const Weight WEIGHT_MAX(INT_MAX - 1, 1);

typedef pair<int, int> Edge;
Edge GetEdge(int u, int v);
bool SortEdge(const Edge& e1, const Edge& e2);
bool SortEdge2(const Edge& e1, const Edge& e2);

typedef vector<int> Clique;
bool CompareClique(const Clique& c1, const Clique& c2);
bool CompareLongClique(const vector<long>& a, const vector<long>& b);

typedef int Attribute;
typedef pair<Weight, int> Index;
typedef flat_hash_map<Edge, int> EdgeMap;
typedef flat_hash_set<Edge> EdgeSet;
typedef map<Weight, vector<Clique>, greater<Weight>> ResultMap;
typedef map<Weight, vector<Edge>, greater<Weight>> WeightEdgeMap;
typedef double Time;
typedef int Memory;

struct Counter
{
    struct value_type { template<typename T> value_type(const T&) { } };
    void push_back(const value_type&) { ++count; }
    int count = 0;
};