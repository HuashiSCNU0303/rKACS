#pragma once

#include "Common.h"

// 计算图中某条边的支持度
int GetEdgeSupport(PUNGraph graph, EdgeMap& sup, Edge e);
void GetAllEdgeSupport(PUNGraph G, EdgeMap& sup, vector<Edge>& supUqEdges, int& supQEdgeCnt, int k);
void GetAllEdgeSupport1(PUNGraph G, EdgeMap& sup, vector<EdgeSet>& nonDecSup);

// 删边时维护支持度，用于truss的相关计算
void MaintainSupport_DelEdge(PUNGraph graph, EdgeMap& sup, vector<EdgeSet>& nonDecSup, Edge e);
void MaintainSupport_AddEdge(PUNGraph graph, EdgeMap& sup, vector<EdgeSet>& nonDecSup, Edge e, int k);

int TrussDecomposition(PUNGraph G, EdgeMap& edgeTrussness);

// 找maximal k-truss，原始方法
PUNGraph GetMaxKTruss(PUNGraph G, int k);
PUNGraph GetMaxKTruss1(PUNGraph G, int k, PUNGraph weightEqPairs);

// 找maximal k-truss，只算新增加的边的支持度
unordered_set<Edge, pair_hash> GetMaxKTrussInc(PUNGraph preKTruss, int k, EdgeMap& sup, vector<EdgeSet>& nonDecSup);

vector<Clique> GetMaxTCKTruss(PUNGraph maximalKTruss);
vector<Clique> GetMaxTCKTruss1(PUNGraph maximalKTruss);
vector<Clique> GetMaxTCKTruss2(PUNGraph maximalKTruss);