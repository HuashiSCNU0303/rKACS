#pragma once

#include "Common.h"

class DynamicClique
{
public:
	PUNGraph G;
	vector<Clique> maxCliques; // maximal cliques of G
	flat_hash_map<int, flat_hash_set<int>> cliquesOfNode; // the indices of maximal cliques containing each node in maxCliques
	flat_hash_set<int> notEmptyIndex, emptyIndex, resultCliques; // the non-empty indices in maxCliques; the empty indices in maxCliques
	
	clock_t startTime = 0, maxTime = 0;
	bool isTimeout = false;
	unordered_map<Clique, int, clique_hash> cliqueToIndex; // the index of each maximal cliques in maxCliques(only used in IMCE)

	void SetMaxTime(clock_t maxTime);

	DynamicClique();
	DynamicClique(PUNGraph G, vector<Clique> maxCliques);
	~DynamicClique();
	
	void AddClique(Clique& clique);
	void AddClique2(Clique& clique);
	void DelClique(int index);

	void AddEdge(int u, int v); // MCMEI
	vector<Clique> AddEdges_SInsert(PUNGraph edgesToInsert); 
	vector<Clique> AddEdges_NIEMC(PUNGraph edgesToInsert); // NIEMC, use RemoveReplicas() to remove the duplicates
	vector<Clique> AddEdges_NIEMCH(PUNGraph edgesToInsert); // NIEMC, use MurmurHash2 to remove the duplicates
	vector<Clique> AddEdges_NIEMCM(PUNGraph edgesToInsert); // NIEMC, use MCE to generate new maximal cliques

	// IMCE
	vector<Clique> AddEdges_IMCE(PUNGraph edgesToInsert); 
	void TTTExcludeEdges(flat_hash_map<int, Clique>& G, Clique& K, Clique& cand, Clique& fini, flat_hash_set<Edge>& e, vector<Clique> &results);
	vector<Clique> IMCE_NewCliques(PUNGraph edgesToInsert);
	void IMCE_SubCliques(PUNGraph edgesToInsert, vector<Clique>& newCliques);
};

