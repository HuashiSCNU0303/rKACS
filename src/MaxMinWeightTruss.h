#pragma once

#include "Common.h"
#include "KRCore.h"
#include "DynamicClique.h"
#include "DynamicTCTruss.h"
#include "TrussComputation.h"

// rKACS Query Algorithms
class MaxMinWeightTruss
{
public:
	PUNGraph G, maxKTruss; // social network G, and T_k[G]
	flat_hash_map<int, int> nodeComAttrCnt; // |A(u) ” Q_W| for each node u （ V(G)
	vector<vector<Edge>> trussnessToEdges; // the corresponding edge set for each trussness value
	vector<Clique> trussnessToNodes; 
	flat_hash_map<Edge, Weight> edgeSims; // attribute similarity of the edges in G
	flat_hash_map<int, vector<Attribute>> attrsOfNode; // attribute set for each node v （ V(G), should be sorted
	flat_hash_map<Attribute, set<int>> nodesInAttr; // vertex set V_w containing attribute w
	WeightEdgeMap edgeWeightInSN; // score(u, v, Q_W) for each edge (u, v) （ E(T_k[G]), sorted in descending order

	bool isTimeout = false;
	clock_t gStartTime = 0; // the start time of the query algorithm

	MaxMinWeightTruss();
	// get T_k[G] by k-truss index
	void GetMaxKTrussByIndex(int k); 
	// compute the attribute similarity between node u and v, 
	// here, parameter node represents node u, parameter nodeAttrs represents the attributes of node v 
	Weight ComputeEdgeSim(flat_hash_set<int>& nodeAttrs, int node);
	void ComputeComAttrCnt(vector<Attribute>& queryAttrs);

	// compute score(u, v, Q_W) for any two nodes u, v （ nodes, and store the result in "results"
	void ComputeWeight(Clique& nodes, WeightEdgeMap& results, bool addEdges = true); 
	// compute edgeWeightInSN
	void ComputeWeight(); 

	// test score(clique, Q_W) == weight?
	bool TestTrussWeightEqual(Clique& clique, Weight weight); 

	void LoadVertexAttribute(const char* fileName);
	int LoadTrussness(const char* fileName);
	int LoadSimilarity(const char* fileName);

	void PrepareQuery(int k, vector<int> queryAttrs);

	// for maximal clique enumeration
	void GetDegeneracyOrder(PUNGraph G, Clique& vertexOrder);
	void BK(flat_hash_map<int, Clique>& G, Clique& R, Clique& P, Clique& X, vector<Clique>& results);
	void GetMaximalCliques(PUNGraph G, vector<Clique>& results);

	// part of the query algorithms
	void QueryByPairs_Basic(int k, int r, const vector<Edge>& pairs, const vector<Index>& weightIndex, ResultMap& results);
	void QueryByPairs_IEMC(int k, int r, DynamicClique* dynamicClique, const vector<Edge>& pairs, const vector<Index>& weightIndex, ResultMap& results, int method);

	// find maximal triangle-connected k-trusses based on the returned maximal cliques
	void GetResults(int k, Weight weight, vector<Clique>& newCliques, PUNGraph updateEdgeGraph, ResultMap& results);
	void GetResults_KRCore(int k, Weight weight, vector<Clique>& newCliques, ResultMap& results);

	// query algorithms, queryAttr = Q_W
	ResultMap rKACS_Basic(int k, int r, vector<int> queryAttrs); // Basic
	ResultMap rKACS_KRCore(int k, int r, vector<int> queryAttrs); // Basic+KRCore
	ResultMap rKACS_IEMC(int k, int r, vector<int> queryAttrs, int method); // Basic+IEMC (the parameter "method" is for specifying the IEMC method(MCMEI/IMCE/NIEMC))
	ResultMap rKACS_Incremental(int k, int r, vector<int> queryAttrs); // Incremental
	ResultMap rKACS_Incremental1(int k, int r, vector<int> queryAttrs, int method); // Incremental, using different methods for IEMC (the parameter "method" is for specifying the IEMC method(MCMEI/IMCE/NIEMC))

	// only to test the running time, and ignore the results
	Time Query_TimeWrapper(int method, int k, int r, vector<int> queryAttrs);

	// VAC specified for |Q| = 0
	// our method can be directly applied with little modification
	ResultMap VACQuery(int k, vector<int>& queryAttrs);

	// Keyword-centric Community
	ResultMap KCQuery(int k, vector<int>& queryAttrs);
	bool KCQuery_IsAllAttrIn(PUNGraph G, vector<int>& queryAttrs);
	int KCQuery_GetKDistU(PUNGraph G, int u, vector<int>& queryAttrs);
	double ComputeKwdC(Clique& queryAttrs, Clique& candCom);
};

