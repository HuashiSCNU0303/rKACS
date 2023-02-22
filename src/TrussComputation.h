#pragma once

#include "Common.h"

int GetEdgeSupport(PUNGraph graph, EdgeMap& sup, Edge e);
void GetAllEdgeSupport(PUNGraph G, EdgeMap& sup, vector<Edge>& supUqEdges, int& supQEdgeCnt, int k);
void GetAllEdgeSupport1(PUNGraph G, EdgeMap& sup, vector<EdgeSet>& nonDecSup);

void MaintainSupport_DelEdge(PUNGraph graph, EdgeMap& sup, vector<EdgeSet>& nonDecSup, Edge e);
void MaintainSupport_AddEdge(PUNGraph graph, EdgeMap& sup, vector<EdgeSet>& nonDecSup, Edge e, int k);

// compute the trussness value for each edge e¡ÊE(G)
int TrussDecomposition(PUNGraph G, EdgeMap& edgeTrussness);

// get the maximal k-truss of G
PUNGraph GetMaxKTruss(PUNGraph G, int k);
PUNGraph GetMaxKTruss1(PUNGraph G, int k, PUNGraph weightEqPairs); // use the new node pairs added in each iteration to skip unpromising truss computation 

// incrementally get the maximal k-truss, only compute the support of newly added edges
unordered_set<Edge, pair_hash> GetMaxKTrussInc(PUNGraph preKTruss, int k, EdgeMap& sup, vector<EdgeSet>& nonDecSup);

// get triangle-connected maximal k-trusses, based on maximal k-truss
vector<Clique> GetMaxTCKTruss(PUNGraph maximalKTruss); // use hashset to store edges
vector<Clique> GetMaxTCKTruss2(PUNGraph maximalKTruss); // use matrix to store edges