#pragma once

#include "Def.h"

int GetIndex(int i, int j, int nodeCnt);
void SplitInt(const string& s, vector<int>& v, const string& c);

// 往有序集合clique中插入顶点
void InsertNode(Clique& clique, int node);
void SeqInsertNode(Clique& clique, int node);
void BinInsertNode(Clique& clique, int node);

// 在有序集合clique中查找顶点
int SearchNode(const Clique& clique, int node);
int BinSearchNode(const Clique& clique, int node, int left, int right);

// 在有序集合clique中删除顶点
bool DeleteNode(Clique& clique, int node);

TIntV ToTIntV(Clique& clique);
Clique ToClique(TIntV& nodes);

bool MaxEval(const Clique& cand, vector<int>& E, PUNGraph G);
bool MaxEval(const Clique& cand, flat_hash_set<int>& E, PUNGraph G);
pair<int, int> GetMinDegNode(flat_hash_map<int, TUNGraph::TNodeI>& nodeInfo, const Clique& cand);

bool MaxEval(const Clique& cand, vector<int>& E, flat_hash_map<int, Clique>& nodeInfo);
bool MaxEval(const Clique& cand, flat_hash_set<int>& E, flat_hash_map<int, Clique>& nodeInfo);
pair<int, int> GetMinDegNode(flat_hash_map<int, Clique>& nodeInfo, const Clique& cand);

void GetNonAscDeg(PUNGraph G, flat_hash_map<int, int>& deg, set<pair<int, int>, greater<pair<int, int>>>& nonAscDeg);
void MaintainDeg_DelNode(flat_hash_map<int, int>& deg, set<pair<int, int>, greater<pair<int, int>>>& nonAscDeg, int node);

void GetNonAscDeg(PUNGraph G, flat_hash_map<int, pair<int, int>>& deg, vector<Clique>& nonAscDeg);
void MaintainDeg_DelNode(flat_hash_map<int, pair<int, int>>& deg, vector<Clique>& nonAscDeg, int node);
pair<int, int> GetMaxDegNode(vector<Clique>& nonAscDeg);

void GetNonAscDeg(PUNGraph G, vector<pair<int, int>>& deg);

PUNGraph GetSubGraph(PUNGraph graph, Clique& nodes);
flat_hash_map<int, Clique>& GetSubGraph1(PUNGraph graph, Clique& nodes, flat_hash_map<int, Clique>& subgraph);

void RemoveReplicas(const vector<Clique>& cliques, vector<int>& candIndexes, int candBegin, int candEnd, int level, int size, vector<int>& cliqueIndexes);

bool IsIn(const Clique& smallSet, const Clique& bigSet);