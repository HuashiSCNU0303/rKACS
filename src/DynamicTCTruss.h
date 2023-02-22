#pragma once

#include "Common.h"

// dynamic computation of triangle-connected k-truss after edge insertion
class DynamicTCTruss
{
private:
	PUNGraph T; // maximal k-truss
	vector<PUNGraph> maxTCTrusses; // the maximal triangle-connected k-trusses in T
	EdgeMap TCTrussOfEdge; // the index of maximal triangle-connected k-truss of each edge
	queue<int> emptyIndex; // the empty indices in maxTCTrusses
	flat_hash_map<int, vector<int>> newPairs;

public:
	DynamicTCTruss();
	~DynamicTCTruss();
	void AddTCTruss(PUNGraph newTCTruss);
	void AddEdges(unordered_set<Edge, pair_hash>& edges); // insert new edges, and maintain the triangle-connected k-trusses
	flat_hash_map<int, vector<int>>& GetNewPairs();
	void ResetNewPairs();

	PUNGraph GetTruss();
};

