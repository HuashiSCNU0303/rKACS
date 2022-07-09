#pragma once

#include "Common.h"

class DynamicTCTruss
{
private:
	PUNGraph T; // maximal k-truss
	vector<PUNGraph> maxTCTrusses; // 当前的maximal triangle-connected k-truss
	EdgeMap TCTrussOfEdge; // 每条边所在的maxTCTruss，唯一（用下标代替指针）
	queue<int> emptyIndex; // 对应位置没有maxTCTruss的数组下标
	flat_hash_map<int, vector<int>> newPairs;

public:
	DynamicTCTruss();
	~DynamicTCTruss();
	void AddTCTruss(PUNGraph newTCTruss);
	void AddEdges(unordered_set<Edge, pair_hash>& edges);
	flat_hash_map<int, vector<int>>& GetNewPairs();
	void ResetNewPairs();

	PUNGraph GetTruss();
};

