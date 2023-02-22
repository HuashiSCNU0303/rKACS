#include "MaxMinWeightTruss.h"

clock_t maintainTime = 0;

MaxMinWeightTruss::MaxMinWeightTruss()
{
	G = TUNGraph::New();
	maxKTruss = TUNGraph::New();
}

void MaxMinWeightTruss::LoadVertexAttribute(const char* fileName)
{
	FILE* attributeF = fopen(fileName, "r");
	while (!feof(attributeF)) {
		int node;
		char attr[20005];
		fscanf(attributeF, "%d\t%s\n", &node, &attr);
		if (!G->IsNode(node)) continue;

		vector<Attribute> vec;
		string str = string(attr);
		if (str.size() == 0) continue;
	
		SplitInt(str, vec, ",");
		for (int i = 0; i < vec.size(); i++) {
			attrsOfNode[node].push_back(vec[i]);
			nodesInAttr[vec[i]].insert(node);
		}
	}
	fclose(attributeF);

	return;
}

int MaxMinWeightTruss::LoadTrussness(const char* fileName)
{
	if (fileName != NULL)
	{
		ifstream file;
		file.open(fileName);
		if (file)
		{
			FILE* FedgeTruss = fopen(fileName, "r");
			bool flag = false;
			while (!feof(FedgeTruss))
			{
				if (!flag)
				{
					int kmax;
					fscanf(FedgeTruss, "%d\n", &kmax);
					trussnessToEdges.resize(kmax + 1);
					trussnessToNodes.resize(kmax + 1);
					flag = true;
				}
				else
				{
					int u, v, truss;
					fscanf(FedgeTruss, "%d\t%d\t%d\n", &u, &v, &truss);
					trussnessToEdges[truss].push_back(GetEdge(u, v));
				}
			}
			flat_hash_set<int> nodes;
			for (int i = trussnessToEdges.size() - 1; i >= 0; i--)
			{
				auto& edges = trussnessToEdges[i];
				sort(edges.begin(), edges.end(), SortEdge);
				for (auto& e : edges)
				{
					if (nodes.insert(e.first).second)
					{
						trussnessToNodes[i].push_back(e.first);
					}
					if (nodes.insert(e.second).second)
					{
						trussnessToNodes[i].push_back(e.second);
					}
				}
			}
			fclose(FedgeTruss);
			return 1;
		}
		else
		{
			EdgeMap trussness;
			int kmax = TrussDecomposition(G, trussness);
			// save the edge trussness to file
			FILE* FedgeTruss = fopen(fileName, "w");
			fprintf(FedgeTruss, "%d\n", kmax);
			for (auto edge_truss : trussness)
			{
				int u = edge_truss.first.first;
				int v = edge_truss.first.second;
				int truss = edge_truss.second;
				fprintf(FedgeTruss, "%d\t%d\t%d\n", u, v, truss);
			}
			fclose(FedgeTruss);
			return 0;
		}
	}
	return -1;
}

int MaxMinWeightTruss::LoadSimilarity(const char* fileName)
{
	if (fileName != NULL)
	{
		ifstream file(fileName, ios::in | ios::binary);
		if (file)
		{
			while (!file.eof())
			{
				int u, v;
				Weight sim;
				file.read((char*)&u, sizeof(int));
				file.read((char*)&v, sizeof(int));
				file.read((char*)&sim, sizeof(Weight));
				edgeSims.emplace(GetEdge(u, v), sim);
			}
			file.close();
			return 1;
		}
		else
		{
			flat_hash_map<Edge, Weight> similarity;

			for (TUNGraph::TEdgeI EI = G->BegEI(); EI != G->EndEI(); EI++)
			{
				int u = EI.GetSrcNId(), v = EI.GetDstNId();
				auto& uAttrs = attrsOfNode.at(u);
				auto& vAttrs = attrsOfNode.at(v);
				Counter c;
				set_intersection(uAttrs.begin(), uAttrs.end(), vAttrs.begin(), vAttrs.end(), back_inserter(c));
				int comAttrCnt = c.count;
				Weight sim = comAttrCnt != 0 ? Weight(comAttrCnt, uAttrs.size() + vAttrs.size() - comAttrCnt) : WEIGHT_ZERO;
				similarity.emplace(GetEdge(u, v), sim);
			}
			
			// save the similarity to file
			ofstream simFile(fileName, ios::out | ios::binary);
			for (auto edge_sim : similarity)
			{
				int u = edge_sim.first.first;
				int v = edge_sim.first.second;
				Weight sim = edge_sim.second;
				simFile.write((char*)&u, sizeof(int));
				simFile.write((char*)&v, sizeof(int));
				simFile.write((char*)&sim, sizeof(Weight));
			}
			simFile.close();
			return 0;
		}
	}
	return -1;
}

Weight MaxMinWeightTruss::ComputeEdgeSim(flat_hash_set<int>& nodeAttrs, int node)
{
	auto& v2Attrs = attrsOfNode.at(node);
	int v2Size = v2Attrs.size();
	int comAttrCnt = 0;
	for (Attribute attr : v2Attrs)
	{
		if (nodeAttrs.find(attr) != nodeAttrs.end())
		{
			comAttrCnt++;
		}
	}
	return comAttrCnt != 0 ? Weight(comAttrCnt, nodeAttrs.size() + v2Size - comAttrCnt) : WEIGHT_ZERO;
}

void MaxMinWeightTruss::ComputeComAttrCnt(vector<Attribute>& queryAttrs)
{
	if (queryAttrs[0] != -1)
	{
		for (Attribute attr : queryAttrs)
		{
			for (int node : nodesInAttr[attr])
			{
				if (maxKTruss->IsNode(node))
				{
					nodeComAttrCnt[node]++;
				}
			}
		}
	}
	else
	{
		for (TUNGraph::TNodeI NI = maxKTruss->BegNI(); NI != maxKTruss->EndNI(); NI++)
		{
			nodeComAttrCnt[NI.GetId()] = 1;
		}
	}
}

void MaxMinWeightTruss::ComputeWeight(Clique& nodes, WeightEdgeMap& results, bool addEdges)
{
	flat_hash_set<Attribute> nodeAttrs;
	int nodeCnt = nodes.size();
	for (int i = 0; i < nodeCnt; i++)
	{
		if ((double)(clock() - gStartTime) > TIME_LIMIT * CLOCKS_PER_SEC)
		{
			isTimeout = true;
			break;
		}
		int u = nodes[i];
		int uComAttrCnt = nodeComAttrCnt[u];
		nodeAttrs.insert(attrsOfNode[u].begin(), attrsOfNode[u].end());
		for (int j = i + 1; j < nodeCnt; j++)
		{
			int v = nodes[j];
			int comAttrCnt = uComAttrCnt + nodeComAttrCnt[v];
			Weight sim = ComputeEdgeSim(nodeAttrs, v);
			if (sim == WEIGHT_ZERO)
			{
				continue;
			}
			Weight weight = sim * comAttrCnt;
			if (addEdges)
			{
				results[weight].push_back(GetEdge(u, v));
			}
			else
			{
				results.emplace(weight, vector<Edge>());
			}
		}
		nodeAttrs.clear();
	}
}

void MaxMinWeightTruss::ComputeWeight()
{
	TIntV nodes;
	maxKTruss->GetNIdV(nodes);
	nodes.Sort();
	int nodeCnt = nodes.Len();

	for (int m = 0; m < nodeCnt; m++)
	{
		int u = nodes[m];
		auto uI = maxKTruss->GetNI(u);
		auto uDeg = uI.GetDeg();
		int uComAttrCnt = nodeComAttrCnt.at(u);
		for (int i = 0; i < uDeg; i++)
		{
			int nbr = uI.GetNbrNId(i);
			if (nbr < u)
			{
				continue;
			}
			int comAttrCnt = uComAttrCnt + nodeComAttrCnt.at(nbr);
			Edge e = make_pair(u, nbr);
			Weight weight = edgeSims[e] * comAttrCnt;
			edgeWeightInSN[weight].push_back(e);
		}
	}
}

bool MaxMinWeightTruss::TestTrussWeightEqual(Clique& clique, Weight weight)
{
	int nodeCnt = clique.size();
	flat_hash_set<Attribute> nodeAttrs;
	for (int i = 0; i < nodeCnt; i++)
	{
		int u = clique[i], uComAttrCnt = nodeComAttrCnt[u];
		nodeAttrs.insert(attrsOfNode[u].begin(), attrsOfNode[u].end());
		for (int j = i + 1; j < nodeCnt; j++)
		{
			int v = clique[j];
			int comAttrCnt = uComAttrCnt + nodeComAttrCnt[v];
			Weight sim = ComputeEdgeSim(nodeAttrs, v);
			Weight w = sim * comAttrCnt;
			if (w == weight)
			{
				return true;
			}
		}
		nodeAttrs.clear();
	}
	return false;
}

void MaxMinWeightTruss::GetMaxKTrussByIndex(int k)
{
	int trussValueCnt = trussnessToEdges.size();
	int nodeCnt = 0;
	for (int i = k; i < trussValueCnt; i++)
	{
		nodeCnt += trussnessToNodes[i].size();
	}

	maxKTruss->Reserve(nodeCnt, -1);
	for (int i = k; i < trussValueCnt; i++)
	{
		auto& nodes = trussnessToNodes[i];
		for (auto node : nodes)
		{
			maxKTruss->AddNodeUnchecked(node);
		}
	}

	for (int i = k; i < trussValueCnt; i++)
	{
		for (auto& edge : trussnessToEdges[i])
		{
			int u = edge.first, v = edge.second;
			maxKTruss->AddEdgeUnchecked(u, v);
		}
	}
	maxKTruss->SortNodeAdjV();
}

void MaxMinWeightTruss::PrepareQuery(int k, vector<int> queryAttrs)
{
	maxKTruss->Clr();
	GetMaxKTrussByIndex(k);

	// compute |A(u) ¡É Q_W| for each node u ¡Ê V(G)
	nodeComAttrCnt.clear();
	ComputeComAttrCnt(queryAttrs);

	// delete the irrelevant nodes
	vector<int> nodes;
	nodes.reserve(maxKTruss->GetNodes());
	for (auto& node_cnt : nodeComAttrCnt)
	{
		if (node_cnt.second > 0)
		{
			nodes.push_back(node_cnt.first);
		}
	}
	sort(nodes.begin(), nodes.end());
	maxKTruss = GetSubGraph(maxKTruss, nodes);

	// delete the edges with attribute similarity = 0
	vector<Edge> edges;
	edges.reserve(2000);
	for (TUNGraph::TEdgeI EI = maxKTruss->BegEI(); EI != maxKTruss->EndEI(); EI++)
	{
		int u = EI.GetSrcNId(), v = EI.GetDstNId();
		Edge e = make_pair(u, v);
		if (edgeSims[e] == WEIGHT_ZERO)
		{
			edges.push_back(e);
		}
	}
	for (auto& e : edges)
	{
		maxKTruss->DelEdge(e.first, e.second);
	}

	maxKTruss = GetMaxKTruss(maxKTruss, k);
	edgeWeightInSN.clear();
}

// Get the degeneracy order of the vertices in G
void MaxMinWeightTruss::GetDegeneracyOrder(PUNGraph G, Clique& vertexOrder)
{
	PUNGraph graph = new TUNGraph(*G);
	flat_hash_map<int, pair<int, int>> deg;
	vector<Clique> nonAscDeg;
	GetNonAscDeg(graph, deg, nonAscDeg);
	vertexOrder.reserve(G->GetNodes());
	while (graph->GetNodes() != 0)
	{
		int minDegNode = -1;
		for (int i = 0; i < nonAscDeg.size(); i++)
		{
			auto& nodes = nonAscDeg[i];
			while (!nodes.empty())
			{
				int u = nodes.back();
				nodes.pop_back();
				if (u != -1)
				{
					minDegNode = u;
					break;
				}
			}
			if (minDegNode != -1)
			{
				break;
			}
		}
		auto uI = graph->GetNI(minDegNode);
		int uDeg = uI.GetDeg();
		for (int i = 0; i < uDeg; i++)
		{
			int nbr = uI.GetNbrNId(i);
			MaintainDeg_DelNode(deg, nonAscDeg, nbr);
		}
		vertexOrder.push_back(minDegNode);
		graph->DelNode(minDegNode);
	}
}

// BK
void MaxMinWeightTruss::BK(flat_hash_map<int, Clique>& G, Clique& R, Clique& P, Clique& X, vector<Clique>& results)
{
	if ((double)(clock() - gStartTime) > TIME_LIMIT * CLOCKS_PER_SEC)
	{
		isTimeout = true;
		return;
	}
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
		auto& uNbrs = G[u];
		Clique R_ = R;
		InsertNode(R_, u);
		Clique P_;
		P_.reserve(P.size());
		set_intersection(P.begin(), P.end(), uNbrs.begin(), uNbrs.end(), back_inserter(P_));
		Clique X_;
		X_.reserve(X.size());
		set_intersection(X.begin(), X.end(), uNbrs.begin(), uNbrs.end(), back_inserter(X_));
		BK(G, R_, P_, X_, results);
		
		DeleteNode(P, u);
		InsertNode(X, u);
		
		if (isTimeout)
		{
			break;
		}
	}
}

void MaxMinWeightTruss::GetMaximalCliques(PUNGraph G, vector<Clique>& results)
{
	flat_hash_map<int, Clique> graph;
	for (TUNGraph::TNodeI NI = G->BegNI(); NI != G->EndNI(); NI++)
	{
		int u = NI.GetId(), uDeg = NI.GetDeg();
		auto& uNbrs = graph[u];
		uNbrs.reserve(uDeg);
		for (int i = 0; i < uDeg; i++)
		{
			uNbrs.push_back(NI.GetNbrNId(i));
		}
	}

	Clique vertexOrder;
	GetDegeneracyOrder(G, vertexOrder);
	int nodeCnt = vertexOrder.size();
	for (int i = 0; i < nodeCnt; i++)
	{
		int v = vertexOrder[i];
		auto& vNbrs = graph[v];
		Clique R, preV, postV;
		R.push_back(v);

		postV.reserve(nodeCnt - i - 1);
		preV.reserve(i);
		for (int j = 0; j < nodeCnt; j++)
		{
			int vertex = vertexOrder[j];
			if (j < i)
			{
				preV.push_back(vertex);
			}
			else if (j > i)
			{
				postV.push_back(vertex);
			}
		}
		sort(postV.begin(), postV.end());
		sort(preV.begin(), preV.end());

		Clique P, X;
		P.reserve(vNbrs.size());
		X.reserve(vNbrs.size());
		set_intersection(postV.begin(), postV.end(), vNbrs.begin(), vNbrs.end(), back_inserter(P));
		set_intersection(preV.begin(), preV.end(), vNbrs.begin(), vNbrs.end(), back_inserter(X));

		BK(graph, R, P, X, results);
		if (isTimeout)
		{
			break;
		}
	}
}

void MaxMinWeightTruss::QueryByPairs_Basic(int k, int r, const vector<Edge>& pairs, const vector<Index>& weightIndex, ResultMap& results)
{
	int pairLIndex = 0;
	PUNGraph G = TUNGraph::New();
	for (auto& weight_cnt : weightIndex)
	{
		// lines 5-6 of Algorithm 1
		Weight weight = weight_cnt.first;
		int cnt = weight_cnt.second;
		PUNGraph updateEdgeGraph = TUNGraph::New();

		for (int i = pairLIndex; i < pairLIndex + cnt; i++)
		{
			int u = pairs[i].first, v = pairs[i].second;
			G->AddEdge2(u, v);
			updateEdgeGraph->AddEdge2(u, v);
		}

		vector<Clique> resultCliques;
		GetMaximalCliques(G, resultCliques);
		if (isTimeout)
		{
			break;
		}

		// lines 7 of Algorithm 1
		GetResults(k, weight, resultCliques, updateEdgeGraph, results);
		pairLIndex = pairLIndex + cnt;

		if (results.size() == r)
		{
			return;
		}
	}
}

void MaxMinWeightTruss::QueryByPairs_IEMC(int k, int r, DynamicClique* dynamicClique, 
	                                          const vector<Edge>& pairs, const vector<Index>& weightIndex, 
	                                          ResultMap& results, int method)
{
	int pairLIndex = 0;
	PUNGraph updateEdgeGraph = TUNGraph::New();
	for (auto& weight_cnt : weightIndex)
	{
		Weight weight = weight_cnt.first;
		int cnt = weight_cnt.second;

		for (int i = pairLIndex; i < pairLIndex + cnt; i++)
		{
			updateEdgeGraph->AddEdge2(pairs[i].first, pairs[i].second);
		}

		PUNGraph updateEdgeGraph1 = new TUNGraph(*updateEdgeGraph);
		vector<Clique> newCliques;
		clock_t startTime = clock();
		if (method == MCMEI)
		{
			newCliques = dynamicClique->AddEdges_SInsert(updateEdgeGraph);
		}
		else if (method == IMCE)
		{
			auto LIMIT = TIME_LIMIT * CLOCKS_PER_SEC - maintainTime;
			dynamicClique->SetMaxTime(LIMIT);
			newCliques = dynamicClique->AddEdges_IMCE(updateEdgeGraph);
		}
		else if (method == NIEMCH)
		{
			auto LIMIT = TIME_LIMIT * CLOCKS_PER_SEC - maintainTime;
			dynamicClique->SetMaxTime(LIMIT);
			newCliques = dynamicClique->AddEdges_NIEMCH(updateEdgeGraph);
		}
		else if (method == NIEMC)
		{
			newCliques = dynamicClique->AddEdges_NIEMC(updateEdgeGraph);
		}
		maintainTime += (clock() - startTime);
		if (maintainTime >= TIME_LIMIT * CLOCKS_PER_SEC)
		{
			updateEdgeGraph.Clr();
			isTimeout = true;
			return;
		}

		GetResults(k, weight, newCliques, updateEdgeGraph1, results);
		pairLIndex = pairLIndex + cnt;

		updateEdgeGraph->Clr();
		updateEdgeGraph1.Clr();
		if (results.size() == r)
		{
			updateEdgeGraph.Clr();
			return;
		}
	}
}

// GetKAC
void MaxMinWeightTruss::GetResults(int k, Weight weight, vector<Clique>& newCliques, PUNGraph updateEdgeGraph, ResultMap& results)
{
	vector<Clique> cliques;
	vector<int> candIndexes, indexes;
	int newCliqueCnt = newCliques.size();
	candIndexes.reserve(newCliqueCnt);
	cliques.reserve(newCliqueCnt);
	indexes.reserve(newCliqueCnt);

	flat_hash_map<int, Clique> cliquesOfSize;
	for (auto& newClique : newCliques)
	{
		PUNGraph cliqueGraph = GetSubGraph(maxKTruss, newClique);
		PUNGraph weightEqPairs = GetSubGraph(updateEdgeGraph, newClique);
		cliqueGraph = GetMaxKTruss1(cliqueGraph, k, weightEqPairs);
		if (cliqueGraph->GetEdges() == 0)
		{
			cliqueGraph.Clr();
			weightEqPairs.Clr();
			continue;
		}
		vector<Clique> TCT = GetMaxTCKTruss2(cliqueGraph);
		for (auto& clique : TCT)
		{
			if (TestTrussWeightEqual(clique, weight))
			{
				if (clique.size() == newClique.size())
				{
					// already maximal
					results[weight].push_back(clique);
					cliques.push_back(clique);
				}
				else
				{
					// line 6
					int newIndex = cliques.size();
					cliques.push_back(clique);
					candIndexes.push_back(newIndex);
					cliquesOfSize[clique.size()].push_back(newIndex);
				}
				indexes.push_back(indexes.size());
			}
		}
		cliqueGraph.Clr();
		weightEqPairs.Clr();
	}

	// remove identical triangle-connected k-truss, you can also use Murmurhash2
	for (auto& size_indexes : cliquesOfSize)
	{
		int size = size_indexes.first;
		auto& indexes_ = size_indexes.second;
		if (indexes_.size() > 1)
		{
			RemoveReplicas(cliques, indexes_, 0, indexes_.size() - 1, 0, size, indexes);
		}
	}
	
	vector<pair<int, int>> sortedIndexes;
	sortedIndexes.reserve(indexes.size());
	for (auto index : indexes)
	{
		if (index == -1)
		{
			continue;
		}
		auto& clique = cliques[index];
		sortedIndexes.push_back(make_pair(clique.size(), index));
	}
	sort(sortedIndexes.begin(), sortedIndexes.end(), SortEdge);

	for (auto& index : candIndexes)
	{
		if (indexes[index] == -1)
		{
			continue;
		}
		auto& clique = cliques[index];
		int comNbrCnt = clique.size();
		bool isMaximal = true;
		for (int m = sortedIndexes.size() - 1; m >= 0; m--)
		{
			int size = sortedIndexes[m].first;
			if (size <= comNbrCnt)
			{
				break;
			}
			auto& bigIndex = sortedIndexes[m].second;
			auto& bigClique = cliques[bigIndex];

			if (IsIn(clique, bigClique)) // such clique T'(bigClique) exists, then T \subseteq T', so T(clique) is not maximal
			{
				isMaximal = false; 
				break;
			}
		}
		if (isMaximal)
		{
			results[weight].push_back(clique);
		}
	}
} 

// GetKAC, slightly modified for KRCore
void MaxMinWeightTruss::GetResults_KRCore(int k, Weight weight, vector<Clique>& newCliques, ResultMap& results)
{
	vector<Clique> cliques;
	vector<int> candIndexes, indexes;
	int newCliqueCnt = newCliques.size();
	candIndexes.reserve(newCliqueCnt);
	cliques.reserve(newCliqueCnt);
	indexes.reserve(newCliqueCnt);

	flat_hash_map<int, Clique> cliquesOfSize;
	for (auto& newClique : newCliques)
	{
		sort(newClique.begin(), newClique.end());
		PUNGraph cliqueGraph = GetSubGraph(maxKTruss, newClique);
		cliqueGraph = GetMaxKTruss(cliqueGraph, k);
		if (cliqueGraph->GetEdges() == 0)
		{
			cliqueGraph.Clr();
			continue;
		}
		vector<Clique> TCT = GetMaxTCKTruss2(cliqueGraph);
		for (auto& clique : TCT)
		{
			if (TestTrussWeightEqual(clique, weight))
			{
				if (clique.size() == newClique.size())
				{
					results[weight].push_back(clique);
					cliques.push_back(clique);
				}
				else
				{
					int newIndex = cliques.size();
					cliques.push_back(clique);
					candIndexes.push_back(newIndex);
					cliquesOfSize[clique.size()].push_back(newIndex);
				}
				indexes.push_back(indexes.size());
			}
		}
		cliqueGraph.Clr();
	}

	for (auto& size_indexes : cliquesOfSize)
	{
		int size = size_indexes.first;
		auto& indexes_ = size_indexes.second;
		if (indexes_.size() > 1)
		{
			RemoveReplicas(cliques, indexes_, 0, indexes_.size() - 1, 0, size, indexes);
		}
	}

	vector<pair<int, int>> sortedIndexes;
	sortedIndexes.reserve(indexes.size());
	for (auto index : indexes)
	{
		if (index == -1)
		{
			continue;
		}
		auto& clique = cliques[index];
		sortedIndexes.push_back(make_pair(clique.size(), index));
	}
	sort(sortedIndexes.begin(), sortedIndexes.end(), SortEdge);

	for (auto& index : candIndexes)
	{
		if (indexes[index] == -1)
		{
			continue;
		}
		auto& clique = cliques[index];
		int comNbrCnt = clique.size();
		bool isMaximal = true;
		for (int m = sortedIndexes.size() - 1; m >= 0; m--)
		{
			int size = sortedIndexes[m].first;
			if (size <= comNbrCnt)
			{
				break;
			}
			auto& bigIndex = sortedIndexes[m].second;
			auto& bigClique = cliques[bigIndex];

			if (IsIn(clique, bigClique))
			{
				isMaximal = false;
				break;
			}
		}
		if (isMaximal)
		{
			results[weight].push_back(clique);
		}
	}
}

ResultMap MaxMinWeightTruss::rKACS_Basic(int k, int r, vector<int> queryAttrs)
{
	PrepareQuery(k, queryAttrs);
	ResultMap results;

	auto maxTCKTruss = GetMaxTCKTruss(maxKTruss);
	WeightEdgeMap edgeWeightMap;
	vector<Edge> pairs;
	vector<Index> weightIndex;
	for (auto& t : maxTCKTruss)
	{
		ComputeWeight(t, edgeWeightMap);
		if (isTimeout)
		{
			return results;
		}
	}

	// organize the node pairs in edgeWeightMap
	// pairs are sorted in descending order of its pairwise weight
	// the weight value are stored in weightIndex, since many pairs share the same weight value
	// each pair in weightIndex is <weight value s, count of node pair with weight = s>
	for (auto weight_edges = edgeWeightMap.begin(); weight_edges != edgeWeightMap.end();)
	{
		Weight weight = weight_edges->first;
		pairs.insert(pairs.end(), weight_edges->second.begin(), weight_edges->second.end());
		weightIndex.push_back(make_pair(weight, weight_edges->second.size()));
		edgeWeightMap.erase(weight_edges++);
	}
	WeightEdgeMap().swap(edgeWeightMap);

	QueryByPairs_Basic(k, r, pairs, weightIndex, results);
	return results;
}

ResultMap MaxMinWeightTruss::rKACS_KRCore(int k, int r, vector<int> queryAttrs)
{
	PrepareQuery(k, queryAttrs);
	ComputeWeight();
	ResultMap results;

	KRCore krCore;
	krCore.initAttrs(maxKTruss, attrsOfNode, queryAttrs);

	auto maxTCKTruss = GetMaxTCKTruss(maxKTruss);
	WeightEdgeMap weights;
	for (auto& t : maxTCKTruss)
	{
		ComputeWeight(t, weights, false);
		if (isTimeout)
		{
			return results;
		}
	}

	for (auto& weight_edges: weights)
	{
		if ((double)(clock() - gStartTime) > TIME_LIMIT * CLOCKS_PER_SEC)
		{
			isTimeout = true;
			break;
		}
		Weight weight = weight_edges.first;
		double value = (double)weight.numerator / (weight.denominator * 2 * queryAttrs.size());

		vector<Edge> edges;
		edges.reserve(maxKTruss->GetEdges() * 2);
		for (auto& weight_edges : edgeWeightInSN)
		{
			if (weight_edges.first < weight)
			{
				break;
			}
			auto& edges_ = weight_edges.second;
			for (auto& edge : edges_)
			{
				int u = edge.first, v = edge.second;
				edges.push_back(make_pair(u, v));
				edges.push_back(make_pair(v, u));
			}
		}
		sort(edges.begin(), edges.end(), SortEdge);

		krCore.init(k - 1, value);
		krCore.dataInput(edges);
		krCore.enumerateKRCore();
		auto newCliques = krCore.getKRCore();
		GetResults_KRCore(k, weight, newCliques, results);
		if (results.size() == r)
		{
			break;
		}
	}
	
	return results;
}

ResultMap MaxMinWeightTruss::rKACS_IEMC(int k, int r, vector<int> queryAttrs, int method)
{
	PrepareQuery(k, queryAttrs);
	ResultMap results;

	auto maxTCKTruss = GetMaxTCKTruss(maxKTruss);
	WeightEdgeMap edgeWeightMap;
	vector<Edge> pairs;
	vector<Index> weightIndex;
	for (auto& t : maxTCKTruss)
	{
		ComputeWeight(t, edgeWeightMap);
		if (isTimeout)
		{
			return results;
		}
	}
	for (auto weight_edges = edgeWeightMap.begin(); weight_edges != edgeWeightMap.end();)
	{
		Weight weight = weight_edges->first;
		pairs.insert(pairs.end(), weight_edges->second.begin(), weight_edges->second.end());
		weightIndex.push_back(make_pair(weight, weight_edges->second.size()));
		edgeWeightMap.erase(weight_edges++);
	}
	WeightEdgeMap().swap(edgeWeightMap);

	DynamicClique* dynamicClique = new DynamicClique();
	QueryByPairs_IEMC(k, r, dynamicClique, pairs, weightIndex, results, method);
	delete dynamicClique;
	return results;
}

ResultMap MaxMinWeightTruss::rKACS_Incremental(int k, int r, vector<int> queryAttrs)
{
	PrepareQuery(k, queryAttrs);
	ComputeWeight(); 

	ResultMap results;

	// nonDecSup here consists of two EdgeSet, nonDecSup[1] stores the edges with support >= k - 2, 
    // and nonDecSup[0] stores the edges with support < k - 2
	EdgeMap sup;
	vector<EdgeSet> nonDecSup;
	nonDecSup.resize(2);
	PUNGraph S = TUNGraph::New(); // we incrementally add edges to build G'_{>=s}
	WeightEdgeMap pairsLeqUb; // P_0
	auto it = edgeWeightInSN.begin();
	Weight ub(3 * queryAttrs.size(), 1); 

	DynamicTCTruss* dynamicTCTruss = new DynamicTCTruss();
	DynamicClique* dynamicClique = new DynamicClique();
	PUNGraph updateEdgeGraph = TUNGraph::New();

	while (it != edgeWeightInSN.end() || !pairsLeqUb.empty())
	{
		Weight w = it != edgeWeightInSN.end() ? it->first : WEIGHT_ZERO;
		vector<Edge> P_1;
		vector<Edge> P_2; 
		vector<Index> P_1_WeightIndex;
		if (it != edgeWeightInSN.end())
		{
			for (auto& e : it->second)
			{
				MaintainSupport_AddEdge(S, sup, nonDecSup, e, k);
			}

			unordered_set<Edge, pair_hash> newEdges = GetMaxKTrussInc(dynamicTCTruss->GetTruss(), k, sup, nonDecSup);
			if (!newEdges.empty())
			{
				dynamicTCTruss->AddEdges(newEdges);
				flat_hash_map<int, vector<int>> newPairs = dynamicTCTruss->GetNewPairs(); // P
				flat_hash_set<int> nodeAttrs;
				// compute the pairwise weight for each pair in P
				for (auto& pair : newPairs)
				{
					int u = pair.first, uComAttrCnt = nodeComAttrCnt[u];
					nodeAttrs.insert(attrsOfNode[u].begin(), attrsOfNode[u].end());
					auto& nodes = pair.second;
					for (int v : nodes)
					{
						int comAttrCnt = uComAttrCnt + nodeComAttrCnt[v];
						Weight sim = ComputeEdgeSim(nodeAttrs, v);
						if (sim == WEIGHT_ZERO)
						{
							continue;
						}
						Weight weight = sim * comAttrCnt;

						if (weight < ub)
						{
							pairsLeqUb[weight].push_back(make_pair(u, v));
						}
						else
						{
							P_2.push_back(make_pair(u, v));
						}
					}
					nodeAttrs.clear();
				}
				dynamicTCTruss->ResetNewPairs();
			}
			it++;
		}

		// obtain P_1 and P_0(= P_0 \ P_1)
		for (auto weight_edges = pairsLeqUb.begin(); weight_edges != pairsLeqUb.end();)
		{
			Weight weight = weight_edges->first;
			if (weight >= w && weight < ub)
			{
				P_1.insert(P_1.end(), weight_edges->second.begin(), weight_edges->second.end());
				P_1_WeightIndex.push_back(make_pair(weight, weight_edges->second.size()));
				pairsLeqUb.erase(weight_edges++);
			}
			else
			{
				weight_edges++;
			}
		}

		for (auto& pair : P_2)
		{
			updateEdgeGraph->AddEdge2(pair.first, pair.second);
		}
		dynamicClique->AddEdges_NIEMCH(updateEdgeGraph);
		updateEdgeGraph->Clr();
		vector<Edge>().swap(P_2);

		// lines 14-18 of Algorithm 6
		QueryByPairs_IEMC(k, r, dynamicClique, P_1, P_1_WeightIndex, results, NIEMCH);

		if (results.size() == r)
		{
			break;
		}

		ub = w;
	}

	S.Clr();
	delete dynamicClique;
	delete dynamicTCTruss;
	return results;
}

ResultMap MaxMinWeightTruss::rKACS_Incremental1(int k, int r, vector<int> queryAttrs, int method)
{
	PrepareQuery(k, queryAttrs);
	ComputeWeight(); 

	ResultMap results;

	EdgeMap sup;
	vector<EdgeSet> nonDecSup;
	nonDecSup.resize(2); 
	PUNGraph S = TUNGraph::New(); 
	WeightEdgeMap pairsLeqUb; 
	auto it = edgeWeightInSN.begin();
	Weight ub(3 * queryAttrs.size(), 1); 

	DynamicTCTruss* dynamicTCTruss = new DynamicTCTruss();
	DynamicClique* dynamicClique = new DynamicClique();
	PUNGraph updateEdgeGraph = TUNGraph::New();

	while (it != edgeWeightInSN.end() || !pairsLeqUb.empty())
	{
		Weight w = it != edgeWeightInSN.end() ? it->first : WEIGHT_ZERO;
		vector<Edge> P_1;
		vector<Edge> P_2; 
		vector<Index> P_1_WeightIndex;
		if (it != edgeWeightInSN.end())
		{
			for (auto& e : it->second)
			{
				MaintainSupport_AddEdge(S, sup, nonDecSup, e, k);
			}

			unordered_set<Edge, pair_hash> newEdges = GetMaxKTrussInc(dynamicTCTruss->GetTruss(), k, sup, nonDecSup);
			if (!newEdges.empty())
			{
				dynamicTCTruss->AddEdges(newEdges);
				flat_hash_map<int, vector<int>> newPairs = dynamicTCTruss->GetNewPairs();
				flat_hash_set<int> nodeAttrs;
				for (auto& pair : newPairs)
				{
					int u = pair.first, uComAttrCnt = nodeComAttrCnt[u];
					nodeAttrs.insert(attrsOfNode[u].begin(), attrsOfNode[u].end());
					auto& nodes = pair.second;
					for (int v : nodes)
					{
						int comAttrCnt = uComAttrCnt + nodeComAttrCnt[v];
						Weight sim = ComputeEdgeSim(nodeAttrs, v);
						if (sim == WEIGHT_ZERO)
						{
							continue;
						}
						Weight weight = sim * comAttrCnt;

						if (weight < ub)
						{
							pairsLeqUb[weight].push_back(make_pair(u, v));
						}
						else
						{
							P_2.push_back(make_pair(u, v));
						}
					}
					nodeAttrs.clear();
				}
				dynamicTCTruss->ResetNewPairs();
			}
			it++;
		}

		for (auto weight_edges = pairsLeqUb.begin(); weight_edges != pairsLeqUb.end();)
		{
			Weight weight = weight_edges->first;
			if (weight >= w && weight < ub)
			{
				P_1.insert(P_1.end(), weight_edges->second.begin(), weight_edges->second.end());
				P_1_WeightIndex.push_back(make_pair(weight, weight_edges->second.size()));
				pairsLeqUb.erase(weight_edges++);
			}
			else
			{
				weight_edges++;
			}
		}

		for (auto& pair : P_2)
		{
			updateEdgeGraph->AddEdge2(pair.first, pair.second);
		}		
		clock_t startTime = clock();
		if (method == MCMEI)
		{
			dynamicClique->AddEdges_SInsert(updateEdgeGraph);
		}
		else if (method == IMCE)
		{
			auto LIMIT = TIME_LIMIT * CLOCKS_PER_SEC - maintainTime;
			dynamicClique->SetMaxTime(LIMIT);
			dynamicClique->AddEdges_IMCE(updateEdgeGraph);
		}
		else if (method == NIEMCH)
		{
			dynamicClique->AddEdges_NIEMCH(updateEdgeGraph);
		}
		maintainTime += (clock() - startTime);
		if (maintainTime >= TIME_LIMIT * CLOCKS_PER_SEC)
		{
			updateEdgeGraph.Clr();
			isTimeout = true;
			return results;
		}

		updateEdgeGraph->Clr();
		vector<Edge>().swap(P_2);

		QueryByPairs_IEMC(k, r, dynamicClique, P_1, P_1_WeightIndex, results, method);

		if (results.size() == r)
		{
			break;
		}

		ub = w;
	}

	S.Clr();
	delete dynamicClique;
	delete dynamicTCTruss;
	return results;
}

Time MaxMinWeightTruss::Query_TimeWrapper(int method, int k, int r, vector<int> queryAttrs)
{
	Time time = 0;
	gStartTime = clock();
	isTimeout = false;
	if (method == BASIC)
	{
		// Basic
		rKACS_Basic(k, r, queryAttrs);
		time = (double)(clock() - gStartTime) / CLOCKS_PER_SEC;
	}
	else if (method == KRCORE)
	{
		// Basic+KRCore
		rKACS_KRCore(k, r, queryAttrs);
		time = (double)(clock() - gStartTime) / CLOCKS_PER_SEC;
	}
	else if (method == NIEMC)
	{
		// Basic+NIEMC
		rKACS_IEMC(k, r, queryAttrs, NIEMCH);
		time = (double)(clock() - gStartTime) / CLOCKS_PER_SEC;
	}
	else if (method == INCREMENTAL)
	{
		// Incremental
		rKACS_Incremental(k, r, queryAttrs);
		time = (double)(clock() - gStartTime) / CLOCKS_PER_SEC;
	}
	else
	{
		// test the time for different IEMC methods(MCMEI/IMCE/NIEMC)
		rKACS_Incremental1(k, r, queryAttrs, method);
		time = (double)maintainTime / CLOCKS_PER_SEC;
	}
	maintainTime = 0;
	return isTimeout? -1: time;
}

ResultMap MaxMinWeightTruss::VACQuery(int k, vector<int>& queryAttrs)
{
	ResultMap results;

	maxKTruss->Clr();
	GetMaxKTrussByIndex(k);

	if (!queryAttrs.empty())
	{
		nodeComAttrCnt.clear();
		ComputeComAttrCnt(queryAttrs);

		vector<int> nodes;
		nodes.reserve(maxKTruss->GetNodes());
		for (auto& node_cnt : nodeComAttrCnt)
		{
			if (node_cnt.second > 0)
			{
				nodes.push_back(node_cnt.first);
			}
		}
		sort(nodes.begin(), nodes.end());
		maxKTruss = GetSubGraph(maxKTruss, nodes);
		maxKTruss = GetMaxKTruss(maxKTruss, k);

		if (maxKTruss->Empty())
		{
			return results;
		}
	}

	nodeComAttrCnt.clear();
	for (auto NI = maxKTruss->BegNI(); NI != maxKTruss->EndNI(); NI++)
	{
		nodeComAttrCnt[NI.GetId()] = 1;
	}

	auto maxTCKTruss = GetMaxTCKTruss(maxKTruss);
	WeightEdgeMap edgeWeightMap;
	vector<Edge> pairs;
	vector<Index> weightIndex;
	for (auto& t : maxTCKTruss)
	{
		ComputeWeight(t, edgeWeightMap);
		if (isTimeout)
		{
			return results;
		}
	}
	for (auto weight_edges = edgeWeightMap.begin(); weight_edges != edgeWeightMap.end();)
	{
		Weight weight = weight_edges->first;
		pairs.insert(pairs.end(), weight_edges->second.begin(), weight_edges->second.end());
		weightIndex.push_back(make_pair(weight, weight_edges->second.size()));
		edgeWeightMap.erase(weight_edges++);
	}
	WeightEdgeMap().swap(edgeWeightMap);

	DynamicClique* dynamicClique = new DynamicClique();
	int pairLIndex = 0;
	PUNGraph updateEdgeGraph = TUNGraph::New();
	for (auto& weight_cnt : weightIndex)
	{
		Weight weight = weight_cnt.first;
		int cnt = weight_cnt.second;

		for (int i = pairLIndex; i < pairLIndex + cnt; i++)
		{
			updateEdgeGraph->AddEdge2(pairs[i].first, pairs[i].second);
		}

		PUNGraph updateEdgeGraph1 = new TUNGraph(*updateEdgeGraph);

		clock_t startTime = clock();
		auto newCliques = dynamicClique->AddEdges_NIEMC(updateEdgeGraph);

		vector<Clique> cliques;
		vector<int> candIndexes, indexes;
		int newCliqueCnt = newCliques.size();
		candIndexes.reserve(newCliqueCnt);
		cliques.reserve(newCliqueCnt);
		indexes.reserve(newCliqueCnt);

		flat_hash_map<int, Clique> cliquesOfSize;
		for (auto& newClique : newCliques)
		{
			PUNGraph cliqueGraph = GetSubGraph(maxKTruss, newClique);
			PUNGraph weightEqPairs = GetSubGraph(updateEdgeGraph1, newClique);
			cliqueGraph = GetMaxKTruss1(cliqueGraph, k, weightEqPairs);
			if (cliqueGraph->GetEdges() == 0)
			{
				cliqueGraph.Clr();
				weightEqPairs.Clr();
				continue;
			}
			vector<Clique> TCT = GetMaxTCKTruss2(cliqueGraph);
			for (auto& clique : TCT)
			{
				if (TestTrussWeightEqual(clique, weight))
				{
					if (clique.size() == newClique.size())
					{
						results[weight].push_back(clique);
						cliques.push_back(clique);
					}
					else
					{
						int newIndex = cliques.size();
						cliques.push_back(clique);
						candIndexes.push_back(newIndex);
						cliquesOfSize[clique.size()].push_back(newIndex);
					}
					indexes.push_back(indexes.size());
				}
			}
			cliqueGraph.Clr();
			weightEqPairs.Clr();
		}

		for (auto& size_indexes : cliquesOfSize)
		{
			int size = size_indexes.first;
			auto& indexes_ = size_indexes.second;
			if (indexes_.size() > 1)
			{
				RemoveReplicas(cliques, indexes_, 0, indexes_.size() - 1, 0, size, indexes);
			}
		}

		vector<pair<int, int>> sortedIndexes;
		sortedIndexes.reserve(indexes.size());
		for (auto index : indexes)
		{
			if (index == -1)
			{
				continue;
			}
			auto& clique = cliques[index];
			sortedIndexes.push_back(make_pair(clique.size(), index));
		}
		sort(sortedIndexes.begin(), sortedIndexes.end(), SortEdge);

		for (auto& index : candIndexes)
		{
			if (indexes[index] == -1)
			{
				continue;
			}
			auto& clique = cliques[index];
			int comNbrCnt = clique.size();
			bool isMaximal = true;
			for (int m = sortedIndexes.size() - 1; m >= 0; m--)
			{
				int size = sortedIndexes[m].first;
				if (size <= comNbrCnt)
				{
					break;
				}
				auto& bigIndex = sortedIndexes[m].second;
				auto& bigClique = cliques[bigIndex];

				if (IsIn(clique, bigClique))
				{
					isMaximal = false; 
					break;
				}
			}
			if (isMaximal)
			{
				results[weight].push_back(clique);
			}
		}

		pairLIndex = pairLIndex + cnt;

		updateEdgeGraph1.Clr();
		if (!results.empty())
		{
			break;
		}
	}
	delete dynamicClique;
	return results;
}

ResultMap MaxMinWeightTruss::KCQuery(int k, vector<int>& queryAttrs)
{
	ResultMap results;

	maxKTruss->Clr();
	GetMaxKTrussByIndex(k);

	nodeComAttrCnt.clear();
	ComputeComAttrCnt(queryAttrs);

	vector<int> nodes;
	nodes.reserve(maxKTruss->GetNodes());
	for (auto& node_cnt : nodeComAttrCnt)
	{
		if (node_cnt.second > 0)
		{
			nodes.push_back(node_cnt.first);
		}
	}
	sort(nodes.begin(), nodes.end());
	maxKTruss = GetSubGraph(maxKTruss, nodes);
	maxKTruss = GetMaxKTruss(maxKTruss, k);

	if (maxKTruss->Empty() || !KCQuery_IsAllAttrIn(maxKTruss, queryAttrs))
	{
		return results;
	}

	int d = -1;
	unordered_map<int, int> kDistMap;
	for (auto NI = maxKTruss->BegNI(); NI != maxKTruss->EndNI(); NI++)
	{
		int dist = KCQuery_GetKDistU(maxKTruss, NI.GetId(), queryAttrs);
		kDistMap.emplace(NI.GetId(), dist);
		if (dist > d)
		{
			d = dist;
		}
	}
	PUNGraph KCResult = new TUNGraph(*maxKTruss);

	PUNGraph G1 = new TUNGraph(*maxKTruss);
	while (KCQuery_IsAllAttrIn(G1, queryAttrs))
	{
		vector<int> delNodes;
		delNodes.reserve(G1->GetNodes());
		for (auto NI = G1->BegNI(); NI != G1->EndNI(); NI++)
		{
			if (kDistMap.at(NI.GetId()) >= d)
			{
				delNodes.push_back(NI.GetId());
			}
		}
		for (auto node : delNodes)
		{
			G1->DelNode(node);
		}
		G1 = GetMaxKTruss(G1, k);
		if (G1->Empty())
		{
			break;
		}

		int d1 = -1;
		for (auto NI = G1->BegNI(); NI != G1->EndNI(); NI++)
		{
			int dist = KCQuery_GetKDistU(G1, NI.GetId(), queryAttrs);
			kDistMap[NI.GetId()] = dist; 
			if (dist > d1)
			{
				d1 = dist;
			}
		}

		if (d1 < d)
		{
			d = d1;
			KCResult = new TUNGraph(*G1);
		}
	}

	vector<int> KCNodes;
	KCNodes.reserve(KCResult->GetNodes());
	for (auto NI = KCResult->BegNI(); NI != KCResult->EndNI(); NI++)
	{
		KCNodes.push_back(NI.GetId());
	}
	sort(KCNodes.begin(), KCNodes.end());

	vector<Clique> temp;
	temp.push_back(KCNodes);
	results.emplace(Weight(d, 1), temp);
	return results;
}

int MaxMinWeightTruss::KCQuery_GetKDistU(PUNGraph G, int u, vector<int>& queryAttrs)
{
	flat_hash_set<int> queryAttrsH;
	for (auto attr : queryAttrs)
	{
		queryAttrsH.emplace(attr);
	}

	queue<pair<int, int>> Q; // <node, distance>
	flat_hash_set<int> visited;
	Q.push(make_pair(u, 0));
	visited.emplace(u);
	while (!Q.empty())
	{
		int node = Q.front().first, dist = Q.front().second;
		Q.pop();

		auto& nodeAttrs = attrsOfNode.at(node);
		for (auto attr : nodeAttrs)
		{
			if (queryAttrsH.find(attr) != queryAttrsH.end())
			{
				queryAttrsH.erase(attr);
			}
		}

		if (queryAttrsH.empty())
		{
			return dist;
		}

		auto nI = G->GetNI(node);
		int nDeg = nI.GetDeg();
		for (int i = 0; i < nDeg; i++)
		{
			int nbr = nI.GetNbrNId(i);
			if (visited.find(nbr) == visited.end())
			{
				Q.push(make_pair(nbr, dist + 1));
				visited.emplace(nbr);
			}
		}
	}

	return INT_MAX;
}

bool MaxMinWeightTruss::KCQuery_IsAllAttrIn(PUNGraph G, vector<int>& queryAttrs)
{
	for (auto attr : queryAttrs)
	{
		bool isAttrIn = false;
		auto& aNodes = nodesInAttr.at(attr);
		for (auto node : aNodes)
		{
			if (G->IsNode(node))
			{
				isAttrIn = true;
				break;
			}
		}
		if (!isAttrIn)
		{
			return false;
		}
	}
	return true;
}

double MaxMinWeightTruss::ComputeKwdC(Clique& queryAttrs, Clique& candCom)
{
	PUNGraph graph = TSnap::GetSubGraph(G, ToTIntV(candCom));
	int d = -1;
	for (auto NI = graph->BegNI(); NI != graph->EndNI(); NI++)
	{
		int dist = KCQuery_GetKDistU(graph, NI.GetId(), queryAttrs);
		if (dist > d)
		{
			d = dist;
		}
	}
	return d;
}