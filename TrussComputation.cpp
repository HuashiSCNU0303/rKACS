#include "TrussComputation.h"

int GetEdgeSupport(PUNGraph graph, EdgeMap& sup, Edge e)
{
	int src = e.first, dst = e.second;
	if (graph->GetNI(src).GetDeg() > graph->GetNI(dst).GetDeg())
	{
		swap(src, dst); 
	}

	int support = 0;
	auto srcI = graph->GetNI(src), dstI = graph->GetNI(dst);
	auto srcDeg = srcI.GetDeg();
	for (int i = 0; i < srcDeg; i++)
	{
		int nbr = srcI.GetNbrNId(i);
		if (dstI.IsNbrNId(nbr)) support++;
	}

	Edge edge = GetEdge(src, dst);
	sup[edge] = support;
	return support;
}

void GetAllEdgeSupport(PUNGraph G, EdgeMap& sup, vector<Edge>& supUqEdges, int& supQEdgeCnt, int k)
{
	sup.reserve(G->GetEdges());
	supUqEdges.reserve(G->GetEdges());

	int nodeCnt = G->GetNodes();
	flat_hash_map<int, Clique> nbrMap;
	Clique nodes;
	nbrMap.reserve(nodeCnt);
	nodes.resize(nodeCnt);
	int j = 0;
	for (TUNGraph::TNodeI NI = G->BegNI(); NI != G->EndNI(); NI++, j++)
	{
		int u = NI.GetId(), uDeg = NI.GetDeg();
		auto& uNbrs = nbrMap[u];
		uNbrs.resize(uDeg);
		for (int i = 0; i < uDeg; i++)
		{
			uNbrs[i] = NI.GetNbrNId(i);
		}
		nodes[j] = u;
	}
	sort(nodes.begin(), nodes.end());

	for (int m = 0; m < nodeCnt; m++)
	{
		int u = nodes[m];
		auto& uNbrs = nbrMap.at(u);
		int uDeg = uNbrs.size();
		for (int i = 0; i < uDeg; i++)
		{
			int v = uNbrs[i];
			if (v <= u)
			{
				continue;
			}
			auto& vNbrs = nbrMap.at(v);
			Counter c;
			set_intersection(uNbrs.begin(), uNbrs.end(), vNbrs.begin(), vNbrs.end(), back_inserter(c));

			Edge edge = make_pair(u, v);
			int support = c.count;
			sup.emplace(edge, support);
			if (support < k - 2)
			{
				supUqEdges.push_back(edge);
			}
			else
			{
				supQEdgeCnt++;
			}
		}
	}
}

void GetAllEdgeSupport(PUNGraph G, vector<Edge>& edgesToCompute, EdgeMap& sup, vector<Edge>& supUqEdges, int& supQEdgeCnt, int k)
{
	sup.reserve(edgesToCompute.size());
	supUqEdges.reserve(edgesToCompute.size());
	flat_hash_map<int, Clique> nbrMap;
	nbrMap.reserve(G->GetNodes());

	for (auto& e : edgesToCompute)
	{
		int u = e.first, v = e.second;
		if (nbrMap.find(u) == nbrMap.end())
		{
			auto uI = G->GetNI(u);
			int uDeg = uI.GetDeg();
			Clique uNbrs;
			uNbrs.resize(uDeg);
			for (int i = 0; i < uDeg; i++)
			{
				uNbrs[i] = uI.GetNbrNId(i);
			}
			nbrMap.emplace(u, uNbrs);
		}
		if (nbrMap.find(v) == nbrMap.end())
		{
			auto vI = G->GetNI(v);
			int vDeg = vI.GetDeg();
			Clique vNbrs;
			vNbrs.resize(vDeg);
			for (int i = 0; i < vDeg; i++)
			{
				vNbrs[i] = vI.GetNbrNId(i);
			}
			nbrMap.emplace(v, vNbrs);
		}
		auto& uNbrs = nbrMap.at(u);
		auto& vNbrs = nbrMap.at(v);

		Counter c;
		set_intersection(uNbrs.begin(), uNbrs.end(), vNbrs.begin(), vNbrs.end(), back_inserter(c));
		
		Edge edge = make_pair(u, v);
		int support = c.count;
		sup.emplace(edge, support);
		if (support < k - 2)
		{
			supUqEdges.push_back(edge);
		}
		else
		{
			supQEdgeCnt++;
		}
	}
}

// for truss decomposition
void GetAllEdgeSupport1(PUNGraph G, EdgeMap& sup, vector<EdgeSet>& nonDecSup)
{
	int maxSup = -1;
	sup.reserve(G->GetEdges());

	int nodeCnt = G->GetNodes();
	unordered_map<int, Clique> nbrMap;
	Clique nodes;
	nbrMap.reserve(nodeCnt);
	nodes.resize(nodeCnt);
	int j = 0;
	for (TUNGraph::TNodeI NI = G->BegNI(); NI != G->EndNI(); NI++, j++)
	{
		int u = NI.GetId(), uDeg = NI.GetDeg();
		Clique uNbrs;
		uNbrs.resize(uDeg);
		for (int i = 0; i < uDeg; i++)
		{
			uNbrs[i] = NI.GetNbrNId(i);
		}
		nbrMap[u] = uNbrs;
		nodes[j] = u;
	}
	sort(nodes.begin(), nodes.end());

	for (int m = 0; m < nodeCnt; m++)
	{
		int u = nodes[m];
		auto uI = G->GetNI(u);
		int uDeg = uI.GetDeg();
		auto& uNbrs = nbrMap[u];
		for (int i = 0; i < uDeg; i++)
		{
			int v = uI.GetNbrNId(i);
			if (v <= u)
			{
				continue;
			}
			auto& vNbrs = nbrMap[v];
			Counter c;
			set_intersection(uNbrs.begin(), uNbrs.end(), vNbrs.begin(), vNbrs.end(), back_inserter(c));

			Edge edge = make_pair(u, v);
			int support = c.count;
			sup[edge] = support;
			if (support > maxSup)
			{
				maxSup = support;
			}
		}
	}

	nonDecSup.resize(maxSup + 1);
	for (auto& e_sup : sup)
	{
		int s = e_sup.second;
		Edge e = e_sup.first;
		nonDecSup[s].insert(e);
	}
}

void MaintainSupport_DelEdge(PUNGraph graph, EdgeMap& sup, vector<EdgeSet>& nonDecSup, Edge e)
{
	int s = sup[e];
	int u = e.first, v = e.second;
	if (graph->GetNI(u).GetDeg() > graph->GetNI(v).GetDeg()) swap(u, v);

	auto uI = graph->GetNI(u), vI = graph->GetNI(v);
	auto uDeg = uI.GetDeg();
	for (int i = 0; i < uDeg; i++)
	{
		int nbr = uI.GetNbrNId(i);
		if (!vI.IsNbrNId(nbr))
		{
			continue;
		}
		Edge nbr_v = GetEdge(nbr, v), nbr_u = GetEdge(nbr, u);

		int& s1 = sup.at(nbr_v);
		nonDecSup[s1].erase(nbr_v);
		s1 -= 1;
		nonDecSup[s1].insert(nbr_v);

		int& s2 = sup.at(nbr_u);
		nonDecSup[s2].erase(nbr_u);
		s2 -= 1;
		nonDecSup[s2].insert(nbr_u);

	}
	sup.erase(e);
	nonDecSup[s].erase(e);
	graph->DelEdge(e.first, e.second);
}

void MaintainSupport_AddEdge(PUNGraph graph, EdgeMap& sup, vector<EdgeSet>& nonDecSup, Edge e, int k)
{
	int u = e.first, v = e.second;
	if (graph->IsNode(u) && graph->IsNode(v))
	{
		if (graph->GetNI(u).GetDeg() > graph->GetNI(v).GetDeg()) swap(u, v);

		auto uI = graph->GetNI(u), vI = graph->GetNI(v);
		auto uDeg = uI.GetDeg();
		for (int i = 0; i < uDeg; i++)
		{
			int nbr = uI.GetNbrNId(i);
			if (!vI.IsNbrNId(nbr))
			{
				continue;
			}
			Edge nbr_v = GetEdge(nbr, v), nbr_u = GetEdge(nbr, u);
			if (sup.find(nbr_v) != sup.end())
			{
				int s1 = sup[nbr_v];
				sup[nbr_v] += 1;
				if (s1 + 1 == k - 2)
				{
					nonDecSup[0].erase(nbr_v);
					nonDecSup[1].insert(nbr_v);
				}
			}

			if (sup.find(nbr_u) != sup.end())
			{
				int s1 = sup[nbr_u];
				sup[nbr_u] += 1;
				if (s1 + 1 == k - 2)
				{
					nonDecSup[0].erase(nbr_u);
					nonDecSup[1].insert(nbr_u);
				}
			}
		}

		// nonDecSup here consists of two EdgeSet, nonDecSup[1] stores the edges with support >= k - 2, 
		// and nonDecSup[0] stores the edges with support < k - 2
		int support = GetEdgeSupport(graph, sup, e);
		if (support >= k - 2)
		{
			nonDecSup[1].insert(e);
		}
		else
		{
			nonDecSup[0].insert(e);
		}
		graph->AddEdge(u, v);
	}
	else
	{
		graph->AddEdge2(u, v);
		sup[e] = 0;
		nonDecSup[0].insert(e);
	}
}

int TrussDecomposition(PUNGraph G, EdgeMap& edgeTrussness)
{
	PUNGraph graph = new TUNGraph(*G);
	EdgeMap sup;
	vector<EdgeSet> nonDecSup;
	GetAllEdgeSupport1(graph, sup, nonDecSup);

	int k = 2, nonDecSupSize = nonDecSup.size();
	while (true)
	{
		int minSup = 0;
		for (; minSup < nonDecSupSize; minSup++)
		{
			if (nonDecSup[minSup].size() > 0)
			{
				break;
			}
		}
		if (minSup == nonDecSupSize)
		{
			break;
		}
		
		if (minSup <= k - 2)
		{
			Edge e = *(nonDecSup[minSup].begin());
			MaintainSupport_DelEdge(graph, sup, nonDecSup, e);
			edgeTrussness[e] = k;
		}
		else
		{
			k += 1;
		}
	}
	graph.Clr();
	return k;
}

PUNGraph GetMaxKTruss(PUNGraph G, int k)
{
	if (G->GetEdges() < k * (k - 1) / 2)
	{
		return TUNGraph::New();
	}
	PUNGraph graph = new TUNGraph(*G);
	EdgeMap sup;
	vector<Edge> supUqEdges;
	int supQEdgeCnt = 0;
	GetAllEdgeSupport(graph, sup, supUqEdges, supQEdgeCnt, k);

	for (int i = 0; i < supUqEdges.size(); i++)
	{
		if (supQEdgeCnt == 0)
		{
			return TUNGraph::New();
		}
		Edge e = supUqEdges[i];
		int u = e.first, v = e.second;
		auto uI = graph->GetNI(u), vI = graph->GetNI(v);
		auto uDeg = uI.GetDeg();
		for (int i = 0; i < uDeg; i++)
		{
			int nbr = uI.GetNbrNId(i);
			if (!vI.IsNbrNId(nbr))
			{
				continue;
			}
			Edge nbr_v = GetEdge(nbr, v), nbr_u = GetEdge(nbr, u);

			int s1 = sup[nbr_v];
			sup[nbr_v] -= 1;
			if (s1 == k - 2)
			{
				supQEdgeCnt--;
				supUqEdges.push_back(nbr_v);
			}

			s1 = sup[nbr_u];
			sup[nbr_u] -= 1;
			if (s1 == k - 2)
			{
				supQEdgeCnt--;
				supUqEdges.push_back(nbr_u);
			}
		}
		sup.erase(e);
		graph->DelEdge(e.first, e.second);
	}

	TSnap::DelZeroDegNodes(graph);
	return graph;
}

PUNGraph GetMaxKTruss1(PUNGraph G, int k, PUNGraph weightEqPairs)
{
	if (G->GetEdges() < k * (k - 1) / 2 || weightEqPairs->GetEdges() == 0)
	{
		return TUNGraph::New();
	}
	PUNGraph graph = new TUNGraph(*G);
	EdgeMap sup;
	vector<Edge> supUqEdges;
	int supQEdgeCnt = 0;
	GetAllEdgeSupport(graph, sup, supUqEdges, supQEdgeCnt, k);

	for (int i = 0; i < supUqEdges.size(); i++)
	{
		if (supQEdgeCnt == 0 || weightEqPairs->GetEdges() == 0)
		{
			return TUNGraph::New();
		}
		Edge e = supUqEdges[i];
		int u = e.first, v = e.second;
		auto uI = graph->GetNI(u), vI = graph->GetNI(v);
		auto uDeg = uI.GetDeg(), vDeg = vI.GetDeg();

		for (int i = 0; i < uDeg; i++)
		{
			int nbr = uI.GetNbrNId(i);

			Edge nbr_v = GetEdge(nbr, v);
			auto s1 = sup.find(nbr_v);
			if (s1 == sup.end())
			{
				continue;
			}

			if (s1->second == k - 2)
			{
				supQEdgeCnt--;
				supUqEdges.push_back(nbr_v);
			}
			s1->second -= 1;

			Edge nbr_u = GetEdge(nbr, u);
			int& s2 = sup.at(nbr_u);
			if (s2 == k - 2)
			{
				supQEdgeCnt--;
				supUqEdges.push_back(nbr_u);
			}
			s2 -= 1;
		}
		sup.erase(e);
		graph->DelEdge(e.first, e.second);
		if (uDeg == 1 && weightEqPairs->IsNode(e.first))
		{
			weightEqPairs->DelNode(e.first);
		}
		if (vDeg == 1 && weightEqPairs->IsNode(e.second))
		{
			weightEqPairs->DelNode(e.second);
		}
	}

	if (weightEqPairs->GetEdges() == 0)
	{
		return TUNGraph::New();
	}

	Clique delNodes;
	delNodes.reserve(graph->GetNodes());
	for (TUNGraph::TNodeI NI = graph->BegNI(); NI != graph->EndNI(); NI++)
	{
		if (NI.GetDeg() == 0)
		{
			delNodes.push_back(NI.GetId());
		}
	}

	for (auto node : delNodes)
	{
		graph->DelNode(node);
		if (weightEqPairs->IsNode(node))
		{
			weightEqPairs->DelNode(node);
		}
	}
	return graph;
}

unordered_set<Edge, pair_hash> GetMaxKTrussInc(PUNGraph graph, int k, EdgeMap& sup, vector<EdgeSet>& nonDecSup)
{
	unordered_set<Edge, pair_hash> result;
	if (nonDecSup[1].size() == 0)
	{
		return result;
	}

	vector<Edge> newEdges, supUqEdges;
	newEdges.reserve(nonDecSup[1].size());
	supUqEdges.reserve(nonDecSup[1].size());
	for (auto& e : nonDecSup[1])
	{
		newEdges.push_back(e);
		graph->AddEdge2(e.first, e.second);
	}

	EdgeMap sup_new;
	int supQEdgeCnt = 0;
	GetAllEdgeSupport(graph, newEdges, sup_new, supUqEdges, supQEdgeCnt, k);

	for (int i = 0; i < supUqEdges.size(); i++)
	{
		if (supQEdgeCnt == 0)
		{
			return result;
		}
		Edge e = supUqEdges[i];
		int u = e.first, v = e.second;
		auto uI = graph->GetNI(u), vI = graph->GetNI(v);
		auto uDeg = uI.GetDeg();
		for (int i = 0; i < uDeg; i++)
		{
			int nbr = uI.GetNbrNId(i);
			if (!vI.IsNbrNId(nbr))
			{
				continue;
			}
			Edge nbr_v = GetEdge(nbr, v), nbr_u = GetEdge(nbr, u);

			auto it_v = sup_new.find(nbr_v);
			if (it_v != sup_new.end())
			{
				int s1 = it_v->second;
				it_v->second -= 1;
				if (s1 == k - 2)
				{
					supQEdgeCnt--;
					supUqEdges.push_back(nbr_v);
				}
			}

			auto it_u = sup_new.find(nbr_u);
			if (it_u != sup_new.end())
			{
				int s1 = it_u->second;
				it_u->second -= 1;
				if (s1 == k - 2)
				{
					supQEdgeCnt--;
					supUqEdges.push_back(nbr_u);
				}
			}
		}
		sup_new.erase(e);
		graph->DelEdge(e.first, e.second);
	}

	for (auto& e_sup : sup_new)
	{
		Edge e = e_sup.first;
		result.insert(e);
		nonDecSup[1].erase(e);
		sup.erase(e);
	}

	return result;
}

vector<Clique> GetMaxTCKTruss(PUNGraph maximalKTruss)
{
	vector<Clique> results;
	int nodeCnt = maximalKTruss->GetNodes();

	unordered_map<int, Clique> nbrMap;
	nbrMap.reserve(maximalKTruss->GetNodes());
	for (TUNGraph::TNodeI NI = maximalKTruss->BegNI(); NI != maximalKTruss->EndNI(); NI++)
	{
		int u = NI.GetId(), uDeg = NI.GetDeg();
		Clique uNbrs;
		uNbrs.resize(uDeg);
		for (int i = 0; i < uDeg; i++)
		{
			uNbrs[i] = NI.GetNbrNId(i);
		}
		nbrMap.emplace(u, uNbrs);
	}

	flat_hash_set<Edge> edges;
	for (TUNGraph::TEdgeI EI = maximalKTruss->BegEI(); EI != maximalKTruss->EndEI(); EI++)
	{
		edges.emplace(make_pair(EI.GetSrcNId(), EI.GetDstNId()));
	}

	while (!edges.empty())
	{
		unordered_set<int> result;
		queue<Edge> Q;
		Edge beginE = *edges.begin();
		Q.push(beginE);
		edges.erase(beginE);
		while (!Q.empty())
		{
			Edge e = Q.front();
			Q.pop();

			int u = e.first, v = e.second;
			result.insert(u);
			result.insert(v);

			auto& uNbrs = nbrMap.at(u);
			auto& vNbrs = nbrMap.at(v);

			vector<int> comNbrs;
			comNbrs.reserve(uNbrs.size());
			set_intersection(uNbrs.begin(), uNbrs.end(), vNbrs.begin(), vNbrs.end(), back_inserter(comNbrs));

			int first1 = 0, last1 = uNbrs.size(), first2 = 0, last2 = vNbrs.size();

			for (auto w : comNbrs)
			{
				Edge w_v = GetEdge(w, v);
				if (edges.erase(w_v) == 1)
				{
					Q.push(w_v);
				}

				Edge w_u = GetEdge(w, u);
				if (edges.erase(w_u) == 1)
				{
					Q.push(w_u);
				}
			}
		}

		Clique c;
		c.reserve(result.size());
		for (auto& n : result)
		{
			c.push_back(n);
		}
		sort(c.begin(), c.end());
		results.push_back(c);
	}
	return results;
}

vector<Clique> GetMaxTCKTruss2(PUNGraph maximalKTruss)
{
	vector<Clique> results;
	int nodeCnt = maximalKTruss->GetNodes();
	TIntV nodes;
	maximalKTruss->GetNIdV(nodes);
	nodes.Sort();
	flat_hash_map<int, int> nodeMap;
	for (int i = 0; i < nodeCnt; i++)
	{
		nodeMap.emplace(nodes[i], i);
	}

	flat_hash_map<int, Clique> nbrMap;
	nbrMap.reserve(nodeCnt);
	for (TUNGraph::TNodeI NI = maximalKTruss->BegNI(); NI != maximalKTruss->EndNI(); NI++)
	{
		int u = nodeMap.at(NI.GetId()), uDeg = NI.GetDeg();
		auto& uNbrs = nbrMap[u];
		uNbrs.resize(uDeg);
		for (int i = 0; i < uDeg; i++)
		{
			uNbrs[i] = nodeMap.at(NI.GetNbrNId(i));
		}
	}

	vector<char> visited; // build matrix to store the visited edges
	visited.resize(nodeCnt * nodeCnt);

	for (auto& src_nbrs : nbrMap)
	{
		int src = src_nbrs.first;
		auto& srcNbrs = src_nbrs.second;
		while (!srcNbrs.empty())
		{
			flat_hash_set<int> result;
			int nbr = srcNbrs[0];
			queue<Edge> Q;
			Q.push(GetEdge(src, nbr));
			while (!Q.empty())
			{
				Edge e = Q.front();
				Q.pop();

				int u = e.first, v = e.second;
				result.insert(u);
				result.insert(v);

				auto& uNbrs = nbrMap.at(u);
				auto& vNbrs = nbrMap.at(v);

				DeleteNode(uNbrs, v);
				DeleteNode(vNbrs, u);

				int first1 = 0, last1 = uNbrs.size(), first2 = 0, last2 = vNbrs.size();
				while (first1 != last1 && first2 != last2)
				{
					int uNbr = uNbrs[first1], vNbr = vNbrs[first2];
					if (uNbr < vNbr)
					{
						first1++;
					}
					else if (uNbr > vNbr)
					{
						first2++;
					}
					else
					{
						int w = uNbr;
						first1++;
						first2++;

						Edge w_v = GetEdge(w, v);
						int index1 = GetIndex(w_v.first, w_v.second, nodeCnt);
						if (visited[index1] == 0)
						{
							Q.push(w_v);
							visited[index1] = 1;
						}

						Edge w_u = GetEdge(w, u);
						int index2 = GetIndex(w_u.first, w_u.second, nodeCnt);
						if (visited[index2] == 0)
						{
							Q.push(w_u);
							visited[index2] = 1;
						}
					}
				}
			}

			Clique c;
			c.reserve(result.size());
			for (auto& n : result)
			{
				c.push_back(nodes[n]);
			}
			sort(c.begin(), c.end());
			results.push_back(c);
		}
	}
	return results;
}