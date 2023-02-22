#include "Common.h"

int GetIndex(int i, int j, int nodeCnt)
{
	return i * nodeCnt + j;
}

void SplitInt(const string& s, vector<int>& v, const string& c)
{
	string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (string::npos != pos2) 
	{
		int num = stoi(s.substr(pos1, pos2 - pos1), nullptr, 10);
		v.push_back(num);
		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length()) 
	{
		int num = stoi(s.substr(pos1), nullptr, 10);
		v.push_back(num);
	}
	return;
}

void InsertNode(Clique& clique, int node)
{
	int size = clique.size();
	if (size == 0)
	{
		clique.push_back(node);
		return;
	}
	BinInsertNode(clique, node);
}

void SeqInsertNode(Clique& clique, int node)
{
	int size = clique.size();
	for (int i = 0; i < size; i++)
	{
		if (clique[i] >= node)
		{
			clique.insert(clique.begin() + i, node);
			return;
		}
	}
	clique.push_back(node);
}

void BinInsertNode(Clique& clique, int node)
{
	int a = 0, b = clique.size() - 1;
	if (b >= 0 && clique[b] < node)
	{
		clique.push_back(node);
		return;
	}

	while (a < b)
	{
		int mid = (a + b) / 2;
		if (clique[mid] >= node)
		{
			b = mid;
		}
		else
		{
			a = mid + 1;
		}
	}
	clique.insert(clique.begin() + a, node);
}

int SearchNode(const Clique& clique, int node)
{
	int size = clique.size();
	if (size == 0)
	{
		return -1;
	}
	return BinSearchNode(clique, node, 0, size - 1);
}

int BinSearchNode(const Clique& clique, int node, int left, int right)
{
	if (left > right)
	{
		return -1;
	}
	int mid = (left + right) / 2;
	if (clique[mid] == node)
	{
		return mid;
	}
	else if (clique[mid] > node)
	{
		return BinSearchNode(clique, node, left, mid - 1);
	}
	else
	{
		return BinSearchNode(clique, node, mid + 1, right);
	}
}

bool DeleteNode(Clique& clique, int node)
{
	int index = SearchNode(clique, node);
	if (index == -1)
	{
		return false;
	}
	clique.erase(clique.begin() + index);
	return true;
}

TIntV ToTIntV(Clique& clique)
{
	TIntV results;
	for (auto node : clique)
	{
		results.Add(node);
	}
	return results;
}

Clique ToClique(TIntV& nodes)
{
	Clique clique;
	int nodeCnt = nodes.Len();
	clique.reserve(nodeCnt);
	for (int i = 0; i < nodeCnt; i++)
	{
		clique.push_back(nodes[i]);
	}
	return clique;
}

bool MaxEval(const Clique& cand, vector<int>& E, PUNGraph G)
{
	int candSize = cand.size();

	flat_hash_map<int, TUNGraph::TNodeI> nodeInfo;
	for (auto node : cand)
	{
		nodeInfo.emplace(node, G->GetNI(node));
	}

	auto node_deg = GetMinDegNode(nodeInfo, cand);
	int minDeg = node_deg.second, minDegNode = node_deg.first;
	if (minDeg == candSize - 1)
	{
		return true;
	}

	vector<int> nbrs, diff;
	auto& uI = nodeInfo[minDegNode];
	nbrs.resize(minDeg);
	diff.reserve(minDeg);
	for (int i = 0; i < minDeg; i++)
	{
		nbrs[i] = uI.GetNbrNId(i);
	}
	set_difference(nbrs.begin(), nbrs.end(), E.begin(), E.end(), back_inserter(diff));

	for (auto node : diff)
	{
		auto nodeI = G->GetNI(node);
		int nodeDeg = nodeI.GetDeg();
		if (nodeDeg >= candSize)
		{
			Clique uNbrs;
			uNbrs.resize(nodeDeg);
			for (int i = 0; i < nodeDeg; i++)
			{
				uNbrs[i] = nodeI.GetNbrNId(i);
			}
			if (IsIn(cand, uNbrs))
			{
				return false; 
			}
		}
	}
	return true;
}

pair<int, int> GetMinDegNode(flat_hash_map<int, TUNGraph::TNodeI>& nodeInfo, const Clique& cand)
{
	int minDeg = INT_MAX, minDegNode = -1;
	for (auto& node : cand)
	{
		int deg = nodeInfo[node].GetDeg();
		if (deg < minDeg)
		{
			minDeg = deg;
			minDegNode = node;
		}
	}
	return make_pair(minDegNode, minDeg);
}

bool MaxEval(const Clique& cand, vector<int>& E, flat_hash_map<int, Clique>& nodeInfo)
{
	int candSize = cand.size();

	auto node_deg = GetMinDegNode(nodeInfo, cand);
	int minDeg = node_deg.second, minDegNode = node_deg.first;
	if (minDeg == candSize - 1)
	{
		return true;
	}

	vector<int> diff;
	auto& nbrs = nodeInfo.at(minDegNode);
	diff.reserve(minDeg);
	set_difference(nbrs.begin(), nbrs.end(), E.begin(), E.end(), back_inserter(diff));

	for (auto u : diff)
	{
		auto& uNbrs = nodeInfo.at(u);
		if (uNbrs.size() >= candSize)
		{
			if (IsIn(cand, uNbrs))
			{
				return false; 
			}
		}
	}
	return true;
}

bool MaxEval(const Clique& cand, flat_hash_set<int>& E, flat_hash_map<int, Clique>& nodeInfo)
{
	int candSize = cand.size();
	auto node_deg = GetMinDegNode(nodeInfo, cand);
	int minDeg = node_deg.second, minDegNode = node_deg.first;
	if (minDeg == candSize - 1)
	{
		return true;
	}

	auto& nbrs = nodeInfo.at(minDegNode);
	for (auto u : nbrs)
	{
		if (E.find(u) == E.end()) // u ¡Ê N(minDegNode) / E
		{
			auto& uNbrs = nodeInfo.at(u);
			if (uNbrs.size() >= candSize)
			{
				if (IsIn(cand, uNbrs))
				{
					return false; 
				}
			}
		}
	}
	return true;
}

pair<int, int> GetMinDegNode(flat_hash_map<int, Clique>& nodeInfo, const Clique& cand)
{
	int minDeg = INT_MAX, minDegNode = -1;
	for (auto& node : cand)
	{
		auto it = nodeInfo.find(node);
		if (it == nodeInfo.end())
		{
			continue;
		}
		int deg = it->second.size();
		if (deg < minDeg)
		{
			minDeg = deg;
			minDegNode = node;
		}
		if (deg == cand.size() - 1)
		{
			return make_pair(minDegNode, minDeg);
		}
	}
	return make_pair(minDegNode, minDeg);
}

void GetNonAscDeg(PUNGraph G, flat_hash_map<int, pair<int, int>>& deg, vector<Clique>& nonAscDeg)
{
	int nodeCnt = G->GetNodes(), maxDeg = -1;
	deg.reserve(nodeCnt);
	for (TUNGraph::TNodeI NI = G->BegNI(); NI != G->EndNI(); NI++)
	{
		int degOfNode = NI.GetDeg(), node = NI.GetId();
		deg.emplace(node, make_pair(degOfNode, -1));
		if (degOfNode > maxDeg)
		{
			maxDeg = degOfNode;
		}
	}
	nonAscDeg.resize(maxDeg + 1);
	for (auto& u_uDeg : deg)
	{
		int u = u_uDeg.first;
		auto& uDeg_index = u_uDeg.second;
		auto& nodesOfDeg = nonAscDeg[uDeg_index.first];
		if (nodesOfDeg.empty())
		{
			nodesOfDeg.reserve(nodeCnt);
		}
		nodesOfDeg.push_back(u);
		u_uDeg.second.second = nodesOfDeg.size() - 1;
	}
}

void MaintainDeg_DelNode(flat_hash_map<int, pair<int, int>>& deg, vector<Clique>& nonAscDeg, int node)
{
	auto& deg_index = deg.at(node);
	int d = deg_index.first, index = deg_index.second;
	nonAscDeg[d][index] = -1;
	d--;
	if (nonAscDeg[d].empty())
	{
		nonAscDeg[d].reserve(nonAscDeg[d + 1].size());
	}
	nonAscDeg[d].push_back(node);
	deg_index.first--;
	deg_index.second = nonAscDeg[d].size() - 1;
}

void GetNonAscDeg(PUNGraph G, vector<pair<int, int>>& deg)
{
	int nodeCnt = G->GetNodes();
	deg.reserve(nodeCnt);
	for (TUNGraph::TNodeI NI = G->BegNI(); NI != G->EndNI(); NI++)
	{
		int degOfNode = NI.GetDeg(), node = NI.GetId();
		deg.push_back(make_pair(degOfNode, node));
	}
	sort(deg.begin(), deg.end(), SortEdge2);
}

PUNGraph GetSubGraph(PUNGraph graph, Clique& nodes)
{
	if (graph->GetNodes() == 0 || nodes.empty())
	{
		return TUNGraph::New();
	}
	PUNGraph NewGraphPt = TUNGraph::New();
	TUNGraph& NewGraph = *NewGraphPt;
	int nodeCnt = nodes.size(), maxNId = nodes.back();
	NewGraph.Reserve(nodeCnt, -1);
	flat_hash_set<int> nodeSet;
	nodeSet.reserve(nodeCnt);

	Clique tempC;
	tempC.reserve(nodeCnt);
	for (auto& n : nodes)
	{
		if (graph->IsNode(n))
		{
			nodeSet.emplace(n);
			NewGraph.AddNodeUnchecked(n);
			tempC.push_back(n);
		}
	}
	for (auto& n : tempC)
	{
		const TUNGraph::TNodeI NI = graph->GetNI(n);
		int nDeg = NI.GetDeg();
		for (int edge = 0; edge < nDeg; edge++) 
		{
			const int nbr = NI.GetNbrNId(edge);
			if (nbr > maxNId)
			{
				break;
			}
			if (nbr > n && nodeSet.find(nbr) != nodeSet.end())
			{
				NewGraph.AddEdgeUnchecked(n, nbr);
			}
		}
	}

	return NewGraphPt;
}

flat_hash_map<int, Clique>& GetSubGraph1(PUNGraph graph, Clique& nodes, flat_hash_map<int, Clique>& subgraph)
{
	if (graph->GetNodes() == 0 || nodes.empty())
	{
		return subgraph;
	}
	int nodeCnt = nodes.size(), maxNId = nodes.back();
	flat_hash_set<int> nodeSet;
	nodeSet.reserve(nodeCnt);
	subgraph.reserve(nodeCnt);

	for (auto& n : nodes)
	{
		nodeSet.emplace(n);
	}

	for (auto& n : nodes)
	{
		const TUNGraph::TNodeI NI = graph->GetNI(n);
		int nDeg = NI.GetDeg();
		vector<int> nbrs;
		nbrs.reserve(nDeg);
		for (int edge = 0; edge < nDeg; edge++)
		{
			const int nbr = NI.GetNbrNId(edge);
			if (nbr > maxNId)
			{
				break;
			}
			if (nodeSet.find(nbr) != nodeSet.end())
			{
				nbrs.push_back(nbr);
			}
		}
		subgraph.emplace(n, std::move(nbrs));
	}

	return subgraph;
}

void RemoveReplicas(const vector<Clique>& cliques, vector<int>& candIndexes, int candBegin, int candEnd,
					 int level, int size, vector<int>& cliqueIndexes)
{
	if (level >= size)
	{
		for (int i = candBegin + 1; i < candEnd + 1; i++)
		{
			int index = candIndexes[i];
			cliqueIndexes[index] = -1; 
		}
		return;
	}

	Clique indexes;
	indexes.reserve(candEnd - candBegin + 1);
	for (int i = candBegin; i < candEnd + 1; i++)
	{
		indexes.push_back(candIndexes[i]);
	}

	sort(indexes.begin(), indexes.end(), [&cliques, &level](const int& i1, const int& i2) {
		return cliques[i1][level] < cliques[i2][level];
	});

	int left = 0, indexCnt = indexes.size();
	for (int i = 0; i < indexCnt; i++)
	{
		if (cliques[indexes[left]][level] == cliques[indexes[i]][level])
		{
			continue;
		}
		else
		{
			if (i - left >= 2)
			{
				RemoveReplicas(cliques, indexes, left, i - 1, level + 1, size, cliqueIndexes);
			}
			left = i;
		}
	}

	if (left < indexCnt - 1)
	{
		RemoveReplicas(cliques, indexes, left, indexCnt - 1, level + 1, size, cliqueIndexes);
	}
}

bool IsIn(const Clique& smallSet, const Clique& bigSet)
{
	int first1 = 0, last1 = bigSet.size(), first2 = 0, last2 = smallSet.size();
	while (first1 != last1 && first2 != last2)
	{
		int uNbr = bigSet[first1], vNbr = smallSet[first2];
		if (uNbr < vNbr)
		{
			first1++;
		}
		else if (uNbr > vNbr)
		{
			return false;
		}
		else
		{
			first1++;
			first2++;
		}
	}
	return first2 == last2;
}