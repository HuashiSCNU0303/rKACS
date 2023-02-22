#include "DynamicClique.h"

DynamicClique::DynamicClique()
{
	G = TUNGraph::New();
}

DynamicClique::DynamicClique(PUNGraph graph, vector<Clique> cliques)
{
	this->G = new TUNGraph(*graph);
	for (auto& clique : cliques)
	{
		AddClique2(clique);
	}
	resultCliques.clear();
}

DynamicClique::~DynamicClique()
{
	G.Clr();
	vector<Clique>().swap(maxCliques);
	flat_hash_map<int, flat_hash_set<int>>().swap(cliquesOfNode);
	flat_hash_set<int>().swap(resultCliques);
	flat_hash_set<int>().swap(notEmptyIndex);
	flat_hash_set<int>().swap(emptyIndex);
}

void DynamicClique::AddClique(Clique& clique)
{
	int index;
	if (emptyIndex.empty())
	{
		index = maxCliques.size();
		maxCliques.push_back(clique);
	}
	else
	{
		index = *(emptyIndex.begin());
		maxCliques[index] = clique;
		emptyIndex.erase(index);
	}
	for (auto node : clique)
	{
		cliquesOfNode[node].insert(index);
	}
	notEmptyIndex.insert(index);
	resultCliques.insert(index);
}

void DynamicClique::AddClique2(Clique& clique)
{
	int index;
	if (emptyIndex.empty())
	{
		index = maxCliques.size();
		maxCliques.push_back(clique);
	}
	else
	{
		index = *(emptyIndex.begin());
		maxCliques[index] = clique;
		emptyIndex.erase(index);
	}
	for (auto node : clique)
	{
		cliquesOfNode[node].insert(index);
	}
	notEmptyIndex.insert(index);
	resultCliques.insert(index);
	cliqueToIndex[clique] = index;
}

void DynamicClique::DelClique(int index)
{
	for (auto node : maxCliques[index])
	{
		cliquesOfNode[node].erase(index);
	}
	notEmptyIndex.erase(index);
	resultCliques.erase(index); 
	emptyIndex.insert(index);
	maxCliques[index].clear();
}

void DynamicClique::AddEdge(int u, int v)
{
	// directly add a new clique {u, v}
	if (!G->IsNode(u) || !G->IsNode(v))
	{
		G->AddEdge2(u, v);
		vector<int> clique;
		if (u < v)
		{
			clique.push_back(u);
			clique.push_back(v);
		}
		else
		{
			clique.push_back(v);
			clique.push_back(u);
		}
		AddClique(clique);
		return;
	}
	G->AddEdge(u, v);
	// assume that d(u, G) < d(v, G)
	if (G->GetNI(u).GetDeg() > G->GetNI(v).GetDeg())
	{
		swap(u, v);
	}

	auto vI = G->GetNI(v);
	auto vDeg = vI.GetDeg();
	flat_hash_set<int> vNbrs;
	for (int i = 0; i < vDeg; i++)
	{
		vNbrs.insert(vI.GetNbrNId(i));
	}

	flat_hash_set<int> cliquesU(cliquesOfNode[u]), cliquesV(cliquesOfNode[v]);
	flat_hash_set<Clique, clique_hash> candCliques;
	candCliques.reserve(cliquesU.size());
	for (auto cliqueIndex : cliquesU)
	{
		auto& clique = maxCliques[cliqueIndex];
		vector<int> cand;
		cand.reserve(MIN(clique.size(), vNbrs.size()));
		for (auto node : clique)
		{
			if (vNbrs.find(node) != vNbrs.end())
			{
				cand.push_back(node);
			}
		}
		InsertNode(cand, v);
		Clique C(clique);
		InsertNode(C, v);
		if (cand == C)
		{
			InsertNode(clique, v);
			cliquesOfNode[v].insert(cliqueIndex);
			resultCliques.insert(cliqueIndex);
		}
		else
		{
			if (MaxEval(cand, cand, G))
			{
				candCliques.emplace(std::move(cand));
			}
		}
	}

	for (auto& c : candCliques)
	{
		Clique clique = c;
		AddClique(clique);
	}

	auto uI = G->GetNI(u);
	auto uDeg = uI.GetDeg();
	flat_hash_set<int> uNbrs;
	for (int i = 0; i < uDeg; i++)
	{
		uNbrs.insert(uI.GetNbrNId(i));
	}

	for (auto cliqueIndex : cliquesV)
	{
		auto& clique = maxCliques[cliqueIndex];
		vector<int> cand;
		cand.reserve(MIN(clique.size(), uNbrs.size()));
		for (auto node : clique)
		{
			if (uNbrs.find(node) != uNbrs.end())
			{
				cand.push_back(node);
			}
		}

		InsertNode(cand, u);
		Clique C(clique);
		InsertNode(C, u);
		if (cand == C)
		{
			DelClique(cliqueIndex);
		}
	}
}

vector<Clique> DynamicClique::AddEdges_SInsert(PUNGraph edgesToInsert)
{
	for (TUNGraph::TEdgeI EI = edgesToInsert->BegEI(); EI != edgesToInsert->EndEI(); EI++)
	{
		int u = EI.GetSrcNId(), v = EI.GetDstNId();
		AddEdge(u, v);
	}

	vector<Clique> newCliques;
	newCliques.reserve(resultCliques.size());
	for (auto cliqueIndex : resultCliques)
	{
		newCliques.push_back(maxCliques[cliqueIndex]);
	}
	resultCliques.clear();
	return newCliques;
}

vector<Clique> DynamicClique::AddEdges_NIEMC(PUNGraph edgesToInsert)
{
	vector<pair<int, int>> deg;
	GetNonAscDeg(edgesToInsert, deg);
	while (edgesToInsert->GetEdges() != 0)
	{
		auto& deg_u = deg.back();
		deg.pop_back();
		int u = deg_u.second;
		auto uI = edgesToInsert->GetNI(u);
		int uDeg = uI.GetDeg();
		if (uDeg <= 0)
		{
			for (int i = 0; i < uDeg; i++)
			{
				int nbr = uI.GetNbrNId(i);
				AddEdge(u, nbr);
			}
		}
		else
		{
			flat_hash_set<int> N; // N(u, G)¡ÈN(u, G1)
			vector<int> newNbrs, N_2; // N_2 is the actual N(u, G1), since |N(u, G) ¡É N(u, G1)| may ¡Ù 0 here
			N_2.reserve(uDeg);
			newNbrs.reserve(uDeg);
			for (int i = 0; i < uDeg; i++)
			{
				int nbr = uI.GetNbrNId(i);
				if (!G->IsNode(nbr))
				{
					newNbrs.push_back(nbr);
				}
				else
				{
					N.insert(nbr);
				}
			}

			flat_hash_map<int, Clique> nbrGraph;
			Clique nbrNodes;

			if (G->IsNode(u))
			{
				flat_hash_set<int> N_1; // N(u, G), not the N_1 in Algorithm 3
				auto uGI = G->GetNI(u);
				int uGDeg = uGI.GetDeg();
				N_1.reserve(uGDeg);
				for (int i = 0; i < uGDeg; i++)
				{
					int nbr = uGI.GetNbrNId(i);
					N_1.insert(nbr);
					N.insert(nbr);
				}

				nbrNodes.reserve(N.size());
				for (auto nbr : N)
				{
					nbrNodes.push_back(nbr);
					if (N_1.find(nbr) == N_1.end()) // N_2 = N \ N_1, to get the actual N(u, G1), so that to eliminate the duplicate node pairs generated by GetNewPairs
					{
						N_2.push_back(nbr);
					}
				}
				sort(nbrNodes.begin(), nbrNodes.end());
				sort(N_2.begin(), N_2.end());
				nbrGraph = GetSubGraph1(G, nbrNodes, nbrGraph); // G[N]

				auto& uCliques = cliquesOfNode.at(u); // the cliques containing u
				for (auto it = uCliques.begin(); it != uCliques.end(); )
				{
					int cliqueIndex = *it;
					auto& clique = maxCliques[cliqueIndex];
					DeleteNode(clique, u);
					if (!MaxEval(clique, N_1, nbrGraph))
					{
						DelClique(cliqueIndex);
						it = uCliques.erase(it);
					}
					else
					{
						InsertNode(clique, u);
						it++;
					}
				}

				for (int node : N_2)
				{
					G->AddEdge(u, node);
				}
			}
			else
			{
				TIntV nbrs;
				nbrNodes.reserve(N.size());
				nbrs.Reserve(N.size());
				for (auto nbr : N)
				{
					nbrNodes.push_back(nbr);
					nbrs.Add(nbr);
				}
				sort(nbrNodes.begin(), nbrNodes.end());
				nbrGraph = GetSubGraph1(G, nbrNodes, nbrGraph);
				G->AddNode(u, nbrs);
				N_2 = std::move(nbrNodes);
			}

			vector<Clique> cliqueDivision;
			cliqueDivision.reserve(N_2.size());
			flat_hash_set<int> cliquesToIntersect;
			for (auto n : N_2)
			{
				auto& cliques = cliquesOfNode.at(n);
				Clique nCliques;
				nCliques.reserve(cliques.size());
				for (auto it = cliques.begin(); it != cliques.end(); it++)
				{
					int index = *it;
					if (cliquesToIntersect.insert(index).second) 
					{
						nCliques.push_back(index);
					}
				}
				cliqueDivision.push_back(std::move(nCliques));
			}

			vector<Clique> candCliques;
			vector<int> candIndexes;
			candCliques.reserve(notEmptyIndex.size());
			candIndexes.reserve(notEmptyIndex.size());
			for (auto& cliques : cliqueDivision)
			{
				unordered_map<int, Clique> cliquesOfSize;
				for (auto& cliqueIndex : cliques)
				{
					auto& clique = maxCliques[cliqueIndex];

					vector<int> intersections;
					intersections.reserve(MIN(clique.size(), N.size()) + 1);
					for (auto node : clique)
					{
						if (N.find(node) != N.end())
						{
							intersections.push_back(node);
						}
					}

					if (intersections.size() == clique.size())
					{
						InsertNode(clique, u);
						cliquesOfNode[u].insert(cliqueIndex);
						resultCliques.insert(cliqueIndex);
					}
					else
					{
						int newIndex = candCliques.size();
						candIndexes.push_back(newIndex);
						cliquesOfSize[intersections.size()].push_back(newIndex);
						candCliques.push_back(std::move(intersections));
					}
				}

				for (auto& size_indexes : cliquesOfSize)
				{
					int size = size_indexes.first;
					auto& indexes = size_indexes.second;
					if (indexes.size() > 1)
					{
						RemoveReplicas(candCliques, indexes, 0, indexes.size() - 1, 0, size, candIndexes);
					}
				}
			}

			for (auto& index : candIndexes)
			{
				if (index == -1)
				{
					continue;
				}
				auto& clique = candCliques[index];
				if (MaxEval(clique, clique, nbrGraph))
				{
					InsertNode(clique, u);
					AddClique(clique);
				}
			}

			for (auto nbr : newNbrs)
			{
				AddEdge(u, nbr);
			}
		}
		edgesToInsert->DelNode(u);
	}

	vector<Clique> newCliques;
	newCliques.reserve(resultCliques.size());
	for (auto cliqueIndex : resultCliques)
	{
		newCliques.push_back(maxCliques[cliqueIndex]);
	}
	resultCliques.clear();
	return newCliques;
}

vector<Clique> DynamicClique::AddEdges_NIEMCH(PUNGraph edgesToInsert)
{
	vector<pair<int, int>> deg;
	GetNonAscDeg(edgesToInsert, deg);
	while (edgesToInsert->GetEdges() != 0)
	{
		auto& deg_u = deg.back();
		deg.pop_back();
		int u = deg_u.second;
		auto uI = edgesToInsert->GetNI(u);
		int uDeg = uI.GetDeg();
		if (uDeg <= 0)
		{
			for (int i = 0; i < uDeg; i++)
			{
				int nbr = uI.GetNbrNId(i);
				AddEdge(u, nbr);
			}
		}
		else
		{
			flat_hash_set<int> N; // N(u, G)¡ÈN(u, G1)
			vector<int> newNbrs, N_2; // N_2 is the actual N(u, G1), since |N(u, G) ¡É N(u, G1)| may ¡Ù 0 here
			N_2.reserve(uDeg);
			newNbrs.reserve(uDeg);
			for (int i = 0; i < uDeg; i++)
			{
				int nbr = uI.GetNbrNId(i);
				if (!G->IsNode(nbr))
				{
					newNbrs.push_back(nbr);
				}
				else
				{
					N.insert(nbr);
				}
			}

			flat_hash_map<int, Clique> nbrGraph;
			Clique nbrNodes;

			// lines 5-9 of Algorithm 3
			if (G->IsNode(u))
			{
				flat_hash_set<int> N_1; // N(u, G), not the N_1 in Algorithm 3
				auto uGI = G->GetNI(u);
				int uGDeg = uGI.GetDeg();
				N_1.reserve(uGDeg);
				for (int i = 0; i < uGDeg; i++)
				{
					int nbr = uGI.GetNbrNId(i);
					N_1.insert(nbr);
					N.insert(nbr);
				}

				nbrNodes.reserve(N.size());
				for (auto nbr : N)
				{
					nbrNodes.push_back(nbr);
					if (N_1.find(nbr) == N_1.end()) // N_2 = N \ N_1, to get the actual N(u, G1), so that to eliminate the duplicate node pairs generated by GetNewPairs
					{
						N_2.push_back(nbr);
					}
				}
				sort(nbrNodes.begin(), nbrNodes.end());
				sort(N_2.begin(), N_2.end());
				nbrGraph = GetSubGraph1(G, nbrNodes, nbrGraph); // G[N]

				auto& uCliques = cliquesOfNode.at(u); // the cliques containing u
				for (auto it = uCliques.begin(); it != uCliques.end(); )
				{
					int cliqueIndex = *it;
					auto& clique = maxCliques[cliqueIndex];
					DeleteNode(clique, u);
					if (!MaxEval(clique, N_1, nbrGraph))
					{
						DelClique(cliqueIndex);
						it = uCliques.erase(it);
					}
					else
					{
						InsertNode(clique, u);
						it++;
					}
				}

				for (int node : N_2)
				{
					G->AddEdge(u, node);
				}
			}
			else
			{
				TIntV nbrs;
				nbrNodes.reserve(N.size());
				nbrs.Reserve(N.size());
				for (auto nbr : N)
				{
					nbrNodes.push_back(nbr);
					nbrs.Add(nbr);
				}
				sort(nbrNodes.begin(), nbrNodes.end());
				nbrGraph = GetSubGraph1(G, nbrNodes, nbrGraph);
				G->AddNode(u, nbrs);
				N_2 = std::move(nbrNodes);
			}

			// line 13 of Algorithm 3
			flat_hash_set<int> cliquesToIntersect;
			for (auto n : N_2)
			{
				auto& cliques = cliquesOfNode.at(n);
				Clique nCliques;
				nCliques.reserve(cliques.size());
				for (auto it = cliques.begin(); it != cliques.end(); it++)
				{
					int index = *it;
					cliquesToIntersect.insert(index);
				}
			}

			flat_hash_set<Clique, clique_hash> candCliques; // employ Murmurhash when inserting new clique
			candCliques.reserve(cliquesToIntersect.size());
			for (auto& cliqueIndex : cliquesToIntersect)
			{
				auto& clique = maxCliques[cliqueIndex];

				vector<int> intersections;
				intersections.reserve(MIN(clique.size(), N.size()) + 1);
				for (auto node : clique)
				{
					if (N.find(node) != N.end())
					{
						intersections.push_back(node);
					}
				}

				if (intersections.size() == clique.size())
				{
					// lines 16-17
					InsertNode(clique, u);
					cliquesOfNode[u].insert(cliqueIndex);
					resultCliques.insert(cliqueIndex);
				}
				else
				{
					// line 18
					candCliques.emplace(std::move(intersections));
				}
			}

			for (auto& c: candCliques)
			{
				Clique clique = c;
				if (MaxEval(clique, clique, nbrGraph))
				{
					InsertNode(clique, u);
					AddClique(clique);
				}
			}

			for (auto nbr : newNbrs)
			{
				AddEdge(u, nbr);
			}
		}
		edgesToInsert->DelNode(u);
	}

	vector<Clique> newCliques;
	newCliques.reserve(resultCliques.size());
	for (auto cliqueIndex : resultCliques)
	{
		newCliques.push_back(maxCliques[cliqueIndex]);
	}
	resultCliques.clear();
	return newCliques;
}

// BK for maximal clique enumeration
void BK(flat_hash_map<int, Clique>& G, Clique& R, Clique& P, Clique& X, const flat_hash_set<int>& fini, vector<Clique>& results)
{
	if (P.empty())
	{
		if (X.empty())
		{
			results.push_back(std::move(R));
		}
		return;
	}

	int maxCnt = -1, pivot = -1;
	for (auto u : P)
	{
		auto& uNbrs = G[u];
		Counter c;
		set_intersection(uNbrs.begin(), uNbrs.end(), P.begin(), P.end(), back_inserter(c));
		if (c.count > maxCnt)
		{
			maxCnt = c.count;
			pivot = u;
		}
	}

	auto& pNbrs = G[pivot];
	Clique cand;
	cand.reserve(P.size());
	set_difference(P.begin(), P.end(), pNbrs.begin(), pNbrs.end(), back_inserter(cand)); // P / N(pivot)

	for (auto u : cand)
	{
		if (fini.find(u) == fini.end())
		{
			auto& uNbrs = G[u];
			Clique R_ = R;
			InsertNode(R_, u);
			Clique P_;
			P_.reserve(P.size());
			set_intersection(P.begin(), P.end(), uNbrs.begin(), uNbrs.end(), back_inserter(P_));
			Clique X_;
			X_.reserve(X.size());
			set_intersection(X.begin(), X.end(), uNbrs.begin(), uNbrs.end(), back_inserter(X_));
			BK(G, R_, P_, X_, fini, results);
		}
		DeleteNode(P, u);
		InsertNode(X, u);
		
	}
}

vector<Clique> DynamicClique::AddEdges_NIEMCM(PUNGraph edgesToInsert)
{
	vector<pair<int, int>> deg;
	GetNonAscDeg(edgesToInsert, deg);
	while (edgesToInsert->GetEdges() != 0)
	{
		auto& deg_u = deg.back();
		deg.pop_back();
		int u = deg_u.second;
		auto uI = edgesToInsert->GetNI(u);
		int uDeg = uI.GetDeg();
		if (uDeg <= 0)
		{
			for (int i = 0; i < uDeg; i++)
			{
				int nbr = uI.GetNbrNId(i);
				AddEdge(u, nbr);
			}
		}
		else
		{
			flat_hash_set<int> N; // N(u, G)¡ÈN(u, G1)
			vector<int> newNbrs, N_2;
			N_2.reserve(uDeg);
			newNbrs.reserve(uDeg);
			for (int i = 0; i < uDeg; i++)
			{
				int nbr = uI.GetNbrNId(i);
				if (!G->IsNode(nbr))
				{
					newNbrs.push_back(nbr);
				}
				else
				{
					N.insert(nbr);
				}
			}

			flat_hash_map<int, Clique> nbrGraph;
			Clique nbrNodes;

			if (G->IsNode(u))
			{
				flat_hash_set<int> N_1; // N(u, G)
				auto uGI = G->GetNI(u);
				int uGDeg = uGI.GetDeg();
				N_1.reserve(uGDeg);
				for (int i = 0; i < uGDeg; i++)
				{
					int nbr = uGI.GetNbrNId(i);
					N_1.insert(nbr);
					N.insert(nbr);
				}

				nbrNodes.reserve(N.size());
				for (auto nbr : N)
				{
					nbrNodes.push_back(nbr);
					if (N_1.find(nbr) == N_1.end())
					{
						N_2.push_back(nbr);
					}
				}
				sort(nbrNodes.begin(), nbrNodes.end());
				sort(N_2.begin(), N_2.end());
				nbrGraph = GetSubGraph1(G, nbrNodes, nbrGraph);

				auto& uCliques = cliquesOfNode.at(u);
				for (auto it = uCliques.begin(); it != uCliques.end(); )
				{
					int cliqueIndex = *it;
					auto& clique = maxCliques[cliqueIndex];
					cliqueToIndex.erase(clique);
					DeleteNode(clique, u);
					if (!MaxEval(clique, N_1, nbrGraph))
					{
						DelClique(cliqueIndex);
						it = uCliques.erase(it);
					}
					else
					{
						InsertNode(clique, u);
						cliqueToIndex[clique] = cliqueIndex;
						it++;
					}
				}

				for (int node : N_2)
				{
					G->AddEdge(u, node);
				}
			}
			else
			{
				TIntV nbrs;
				nbrNodes.reserve(N.size());
				nbrs.Reserve(N.size());
				for (auto nbr : N)
				{
					nbrNodes.push_back(nbr);
					nbrs.Add(nbr);
					N_2.push_back(nbr);
				}
				sort(nbrNodes.begin(), nbrNodes.end());
				sort(N_2.begin(), N_2.end());
				nbrGraph = GetSubGraph1(G, nbrNodes, nbrGraph);
				G->AddNode(u, nbrs);
			}

			vector<Clique> cliques;
			flat_hash_set<int> X_;
			X_.reserve(N_2.size());
			for (auto v : N_2)
			{
				Clique R, P, X;
				R.push_back(v);
				P = nbrGraph.at(v);
				BK(nbrGraph, R, P, X, X_, cliques);
				X_.emplace(v);
			}

			for (auto& clique : cliques)
			{
				if (clique.empty())
				{
					continue;
				}
				auto it = cliqueToIndex.find(clique);
				if (it != cliqueToIndex.end())
				{
					DelClique(it->second);
					cliqueToIndex.erase(it);
				}

				InsertNode(clique, u);
				AddClique2(clique);
			}

			for (auto nbr : newNbrs)
			{
				G->AddEdge2(u, nbr);
				vector<int> clique;
				if (u < nbr)
				{
					clique.push_back(u);
					clique.push_back(nbr);
				}
				else
				{
					clique.push_back(nbr);
					clique.push_back(u);
				}
				AddClique2(clique);
			}
		}
		edgesToInsert->DelNode(u);
	}

	vector<Clique> newCliques;
	newCliques.reserve(resultCliques.size());
	for (auto cliqueIndex : resultCliques)
	{
		newCliques.push_back(maxCliques[cliqueIndex]);
	}
	resultCliques.clear();
	return newCliques;
}

void DynamicClique::SetMaxTime(clock_t maxTime)
{
	this->maxTime = maxTime;
}

vector<Clique> DynamicClique::AddEdges_IMCE(PUNGraph edgesToInsert)
{
	isTimeout = false;
	startTime = clock();

	vector<Clique> newCliques = IMCE_NewCliques(edgesToInsert);
	if (isTimeout)
	{
		return newCliques;
	}

	IMCE_SubCliques(edgesToInsert, newCliques);
	for (auto& newClique : newCliques)
	{
		AddClique2(newClique);
	}
	edgesToInsert->Clr();
	return newCliques;
}

void DynamicClique::TTTExcludeEdges(flat_hash_map<int, Clique>& G, Clique& K, Clique& cand, Clique& fini, flat_hash_set<Edge>& edges, vector<Clique>& results)
{
	if (clock() - startTime > maxTime)
	{
		isTimeout = true;
		return;
	}
	if (cand.empty() && fini.empty())
	{
		results.push_back(K);
		return;
	}

	int maxCnt = -1, pivot = -1;
	for (auto u: cand)
	{
		auto& uNbrs = G[u];
		Counter c;
		set_intersection(uNbrs.begin(), uNbrs.end(), cand.begin(), cand.end(), back_inserter(c));
		if (c.count > maxCnt)
		{
			maxCnt = c.count;
			pivot = u;
		}
	}
	for (auto u : fini)
	{
		auto& uNbrs = G[u];
		Counter c;
		set_intersection(uNbrs.begin(), uNbrs.end(), cand.begin(), cand.end(), back_inserter(c));
		if (c.count > maxCnt)
		{
			maxCnt = c.count;
			pivot = u;
		}
	}

	vector<int> ext;
	auto& pNbrs = G[pivot];
	ext.reserve(MAX(cand.size(), pNbrs.size()));
	set_difference(cand.begin(), cand.end(), pNbrs.begin(), pNbrs.end(), back_inserter(ext));

	for (auto q : ext)
	{
		vector<int> K_q = K;
		InsertNode(K_q, q);
		bool skipFlag = false;
		int KqCnt = K_q.size();
		for (int i = 0; i < KqCnt; i++)
		{
			int u = K_q[i];
			for (int j = i + 1; j < KqCnt; j++)
			{
				int v = K_q[j];
				if (edges.find(make_pair(u, v)) != edges.end())
				{
					skipFlag = true;
					break;
				}
			}
			if (skipFlag)
			{
				break;
			}
		}
		if (skipFlag)
		{
			DeleteNode(cand, q);
			InsertNode(fini, q);
			continue;
		}
		auto& qNbrs = G[q];
		vector<int> cand_q, fini_q;
		cand_q.reserve(cand.size());
		fini_q.reserve(fini.size());
		set_intersection(cand.begin(), cand.end(), qNbrs.begin(), qNbrs.end(), back_inserter(cand_q));
		set_intersection(fini.begin(), fini.end(), qNbrs.begin(), qNbrs.end(), back_inserter(fini_q));

		TTTExcludeEdges(G, K_q, cand_q, fini_q, edges, results);

		DeleteNode(cand, q);
		InsertNode(fini, q);

		if (isTimeout)
		{
			break;
		}
	}

}

vector<Clique> DynamicClique::IMCE_NewCliques(PUNGraph edgesToInsert)
{
	vector<Edge> newEdges;
	newEdges.reserve(edgesToInsert->GetEdges());
	flat_hash_set<Edge> finishedEdges;
	finishedEdges.reserve(edgesToInsert->GetEdges());

	vector<Clique> results;
	for (TUNGraph::TEdgeI EI = edgesToInsert->BegEI(); EI != edgesToInsert->EndEI(); EI++)
	{
		int src = EI.GetSrcNId(), dst = EI.GetDstNId();
		G->AddEdge2(src, dst);
		newEdges.push_back(GetEdge(src, dst));
	}
	
	flat_hash_map<int, Clique> nbrMap;
	for (auto& e : newEdges)
	{
		int u = e.first, v = e.second;
		if (nbrMap.find(u) == nbrMap.end())
		{
			vector<int> uNbrs;
			auto uI = G->GetNI(u);
			int uDeg = uI.GetDeg();
			uNbrs.reserve(uDeg);
			for (int i = 0; i < uDeg; i++)
			{
				uNbrs.push_back(uI.GetNbrNId(i));
			}
			nbrMap.emplace(u, std::move(uNbrs));
		}

		if (nbrMap.find(v) == nbrMap.end())
		{
			vector<int> vNbrs;
			auto vI = G->GetNI(v);
			int vDeg = vI.GetDeg();
			vNbrs.reserve(vDeg);
			for (int i = 0; i < vDeg; i++)
			{
				vNbrs.push_back(vI.GetNbrNId(i));
			}
			nbrMap.emplace(v, std::move(vNbrs));
		}
	}

	for (auto& e : newEdges)
	{
		int u = e.first, v = e.second;
		auto& uNbrs = nbrMap[u];
		auto& vNbrs = nbrMap[v];

		vector<int> comNbrs, cand;
		comNbrs.reserve(MIN(uNbrs.size(), vNbrs.size()) + 2);
		set_intersection(uNbrs.begin(), uNbrs.end(), vNbrs.begin(), vNbrs.end(), back_inserter(comNbrs));
		cand = comNbrs;
		InsertNode(comNbrs, u);
		InsertNode(comNbrs, v);

		flat_hash_map<int, Clique> nbrGraph;
		nbrGraph = GetSubGraph1(G, comNbrs, nbrGraph);
		vector<int> K;
		K.push_back(u < v ? u : v);
		K.push_back(u < v ? v : u);

		vector<int> fini;
		TTTExcludeEdges(nbrGraph, K, cand, fini, finishedEdges, results);
		if (isTimeout)
		{
			break;
		}

		finishedEdges.emplace(e);
	}

	return results;
}

void DynamicClique::IMCE_SubCliques(PUNGraph edgesToInsert, vector<Clique>& newCliques)
{
	for (auto& c : newCliques)
	{
		if (clock() - startTime > maxTime)
		{
			isTimeout = true;
			return;
		}
		flat_hash_set<Clique, clique_hash> S;
		S.emplace(c);
		int nodeCnt = c.size();
		vector<Edge> intersections;
		intersections.reserve(MIN(edgesToInsert->GetEdges(), 0.5 * nodeCnt * nodeCnt));
		for (int i = 0; i < nodeCnt; i++)
		{
			int u = c[i];
			for (int j = i + 1; j < nodeCnt; j++)
			{
				int v = c[j];
				if (edgesToInsert->IsEdge(u, v))
				{
					intersections.push_back(GetEdge(u, v));
				}
			}
		}
		for (auto& e : intersections)
		{
			int u = e.first, v = e.second;
			vector<Clique> candInsertion;
			candInsertion.reserve(S.size());
			for (flat_hash_set<Clique, clique_hash>::iterator it = S.begin(); it != S.end();)
			{
				int uIndex = SearchNode(*it, u);
				if (uIndex != -1)
				{
					int vIndex = SearchNode(*it, v);
					if (vIndex != -1)
					{
						Clique c1 = *it, c2 = *it;
						c1.erase(c1.begin() + uIndex);
						c2.erase(c2.begin() + vIndex);
						candInsertion.push_back(std::move(c1));
						candInsertion.push_back(std::move(c2));
						it = S.erase(it);
					}
					else
					{
						it++;
					}
				}
				else
				{
					it++;
				}
			}
			S.reserve(S.size() + candInsertion.size());
			for (auto& clique : candInsertion)
			{
				S.emplace(std::move(clique));
			}
		}
		for (auto& c_ : S)
		{
			auto it = cliqueToIndex.find(c_);
			if (it != cliqueToIndex.end())
			{
				DelClique(it->second);
				cliqueToIndex.erase(it);
			}
		}
	}
}