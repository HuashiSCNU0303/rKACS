#pragma once

#include "Common.h"
#include "Trie.h"

class DynamicClique
{
public:
	PUNGraph G;
	vector<Clique> maxCliques; // 当前的maximal clique
	flat_hash_map<int, flat_hash_set<int>> cliquesOfNode; // 每个顶点所在的maximal clique的集合（用下标代替指针）
	flat_hash_set<int> notEmptyIndex, emptyIndex, resultCliques; // 对应位置有clique的数组下标，没有clique的数组下标，以及添加/删除边后要返回的clique集合
	
	clock_t startTime = 0, maxTime = 0;
	bool isTimeout = false;
	unordered_map<Clique, int, clique_hash> cliqueToIndex; // 每个maximal clique对应的index，使用MurmurHash2，只有IMCE使用

	void SetMaxTime(clock_t maxTime);

	DynamicClique();
	DynamicClique(PUNGraph G, vector<Clique> maxCliques);
	~DynamicClique();
	
	void AddClique(Clique& clique);
	void AddClique2(Clique& clique);
	void DelClique(int index);

	void AddEdge(int u, int v);
	vector<Clique> AddEdges_SInsert(PUNGraph edgesToInsert); // 单条边插入
	vector<Clique> AddEdges_BInsert(PUNGraph edgesToInsert); // 我们的方法
	vector<Clique> AddEdges_BInsertH(PUNGraph edgesToInsert); // 使用哈希去重
	vector<Clique> AddEdges_BInsertM(PUNGraph edgesToInsert); // 找新clique的时候使用MCE，小的时候还可以，大的时候不行了，就不用这个方法了

	vector<Clique> AddEdges_IMCE(PUNGraph edgesToInsert); // VLDBJ 2019，感觉有点吹牛逼
	void TTTExcludeEdges(flat_hash_map<int, Clique>& G, Clique& K, Clique& cand, Clique& fini, flat_hash_set<Edge>& e, vector<Clique> &results);
	vector<Clique> IMCE_NewCliques(PUNGraph edgesToInsert);
	void IMCE_SubCliques(PUNGraph edgesToInsert, vector<Clique>& newCliques);
	

	vector<Clique> DelEdges(PUNGraph edgesToDelete);
};

