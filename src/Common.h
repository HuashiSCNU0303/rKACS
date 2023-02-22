#pragma once

#include "Def.h"

int GetIndex(int i, int j, int nodeCnt);
void SplitInt(const string& s, vector<int>& v, const string& c);

// insert node in the sorted node set
void InsertNode(Clique& clique, int node);
void SeqInsertNode(Clique& clique, int node);
void BinInsertNode(Clique& clique, int node);

// search node in the sorted node set
int SearchNode(const Clique& clique, int node);
int BinSearchNode(const Clique& clique, int node, int left, int right);

// delete node in the sorted node set
bool DeleteNode(Clique& clique, int node);

TIntV ToTIntV(Clique& clique);
Clique ToClique(TIntV& nodes);

// MaxEval, using PUNGraph to model graph
bool MaxEval(const Clique& cand, vector<int>& E, PUNGraph G);
pair<int, int> GetMinDegNode(flat_hash_map<int, TUNGraph::TNodeI>& nodeInfo, const Clique& cand);

// MaxEval, using hashtable to model graph
bool MaxEval(const Clique& cand, vector<int>& E, flat_hash_map<int, Clique>& nodeInfo);
bool MaxEval(const Clique& cand, flat_hash_set<int>& E, flat_hash_map<int, Clique>& nodeInfo);
pair<int, int> GetMinDegNode(flat_hash_map<int, Clique>& nodeInfo, const Clique& cand);

void GetNonAscDeg(PUNGraph G, flat_hash_map<int, pair<int, int>>& deg, vector<Clique>& nonAscDeg);
void MaintainDeg_DelNode(flat_hash_map<int, pair<int, int>>& deg, vector<Clique>& nonAscDeg, int node);

// get the degree of each node, and sorted in descending order
void GetNonAscDeg(PUNGraph G, vector<pair<int, int>>& deg);

// an modified implementation of TSnap::GetSubGraph
PUNGraph GetSubGraph(PUNGraph graph, Clique& nodes);
flat_hash_map<int, Clique>& GetSubGraph1(PUNGraph graph, Clique& nodes, flat_hash_map<int, Clique>& subgraph);

// remove the duplicates in cliques, deprecated
void RemoveReplicas(const vector<Clique>& cliques, vector<int>& candIndexes, int candBegin, int candEnd, int level, int size, vector<int>& cliqueIndexes);

// judge smallSet \subseteq bigSet?
bool IsIn(const Clique& smallSet, const Clique& bigSet);