#include "DynamicTCTruss.h"

DynamicTCTruss::DynamicTCTruss()
{
	T = TUNGraph::New();
}

DynamicTCTruss::~DynamicTCTruss()
{
	T.Clr();
	for (auto &graph : maxTCTrusses)
	{
		graph.Clr();
	}
	queue<int>().swap(emptyIndex);
	EdgeMap().swap(TCTrussOfEdge);
	vector<PUNGraph>().swap(maxTCTrusses);
}

void DynamicTCTruss::AddTCTruss(PUNGraph newTCTruss)
{
	int index = emptyIndex.empty() ? maxTCTrusses.size() : emptyIndex.front();
	if (!emptyIndex.empty())
	{
		emptyIndex.pop();
		maxTCTrusses[index] = newTCTruss;
	}
	else
	{
		maxTCTrusses.push_back(newTCTruss);
	}

	for (TUNGraph::TEdgeI EI = newTCTruss->BegEI(); EI != newTCTruss->EndEI(); EI++)
	{
		TCTrussOfEdge[make_pair(EI.GetSrcNId(), EI.GetDstNId())] = index;
	}
}

void DynamicTCTruss::AddEdges(unordered_set<Edge, pair_hash>& edges)
{
	int index = maxTCTrusses.size();
	PUNGraph TCTrussGraph = TUNGraph::New();
	vector<vector<Edge>> newTCs;

	// do BFS, when traversing to existing edges, stop further traverse
	unordered_set<int> mergeTCTrusses;
	flat_hash_set<Edge> visited;
	while (!edges.empty())
	{
		vector<Edge> result;
		queue<Edge> Q;
		Edge beginE = *edges.begin();
		Q.push(beginE);
		visited.insert(beginE);
		while (!Q.empty())
		{
			Edge e = Q.front();
			Q.pop();

			auto it = TCTrussOfEdge.find(e);
			if (it != TCTrussOfEdge.end())
			{
				mergeTCTrusses.insert(it->second);
				continue;
			}

			int u = e.first, v = e.second;
			if (T->GetNI(u).GetDeg() > T->GetNI(v).GetDeg())
			{
				swap(u, v);
			}

			result.push_back(e);
			edges.erase(e);

			vector<int> comNeighbors;
			auto uI = T->GetNI(u), vI = T->GetNI(v);
			auto uDeg = uI.GetDeg();
			comNeighbors.reserve(uDeg);
			for (int i = 0; i < uDeg; i++)
			{
				int nbr = uI.GetNbrNId(i);
				if (vI.IsNbrNId(nbr))
				{
					comNeighbors.push_back(nbr);
				}
			}

			for (auto w : comNeighbors)
			{
				Edge w_v = GetEdge(w, v);
				if (visited.insert(w_v).second) 
				{
					Q.push(w_v);
				}
				Edge w_u = GetEdge(w, u);
				if (visited.insert(w_u).second)
				{
					Q.push(w_u);
				}
			}
		}

		newTCs.push_back(result);
		TCTrussGraph->AddNode(index);
		for (int TCIndex : mergeTCTrusses)
		{
			if (!TCTrussGraph->IsNode(TCIndex))
			{
				TCTrussGraph->AddNode(TCIndex);
			}
			TCTrussGraph->AddEdge(index, TCIndex);
		}
		index++;

		mergeTCTrusses.clear();
		visited.clear();
	}

	vector<vector<int>> newTCIndexes;
	typename PUNGraph::TObj::TNodeI NI;
	int graphSize = TCTrussGraph->GetNodes();
	THashSet<TInt> VisitedNId(graphSize + 1);
	TSnapQueue<int> NIdQ(graphSize + 1);
	vector<int> CcNIdV(1, 0);
	for (NI = TCTrussGraph->BegNI(); NI < TCTrussGraph->EndNI(); NI++)
	{
		if (NI.GetDeg() == 0)
		{
			const int NId = NI.GetId();
			VisitedNId.AddKey(NId);
			CcNIdV[0] = NId;
			newTCIndexes.push_back(CcNIdV);
		}
	}

	for (NI = TCTrussGraph->BegNI(); NI < TCTrussGraph->EndNI(); NI++)
	{
		const int NId = NI.GetId();
		if (!VisitedNId.IsKey(NId))
		{
			VisitedNId.AddKey(NId);
			NIdQ.Clr(false);
			NIdQ.Push(NId);
			CcNIdV.clear();
			CcNIdV.push_back(NId);
			while (!NIdQ.Empty())
			{
				const typename PUNGraph::TObj::TNodeI Node = TCTrussGraph->GetNI(NIdQ.Top());
				NIdQ.Pop();
				for (int e = 0; e < Node.GetOutDeg(); e++)
				{
					const int OutNId = Node.GetOutNId(e);
					if (!VisitedNId.IsKey(OutNId))
					{
						NIdQ.Push(OutNId);
						VisitedNId.AddKey(OutNId);
						CcNIdV.push_back(OutNId);
					}
				}
			}
			sort(CcNIdV.begin(), CcNIdV.end());
			newTCIndexes.push_back(CcNIdV);
		}
	}

	// Algorithm 5

	int TCTrussSize = maxTCTrusses.size();
	for (auto& newTC : newTCIndexes)
	{
		PUNGraph result = TUNGraph::New(); // T_N
		vector<int> newNodes; // N
		vector<vector<int>> preTCTrussNodes;
		int newTCSize = newTC.size();
		preTCTrussNodes.reserve(newTCSize);
		// we do the mapping of new and old triangle-connected k-trusses when dynamically computing the new triangle-connected k-trusses
		// so we directly obtain the corresponding \mathcal{T}_O
		for (auto& index : newTC)
		{
			if (index < TCTrussSize)
			{
				PUNGraph graph = maxTCTrusses[index];
				for (TUNGraph::TEdgeI EI = graph->BegEI(); EI != graph->EndEI(); EI++)
				{
					int u = EI.GetSrcNId(), v = EI.GetDstNId();
					result->AddNodeUnchecked(u);
					result->AddNodeUnchecked(v);
					result->AddEdgeUnchecked(u, v);
				}
				emptyIndex.push(index);
				vector<int> nodes;
				nodes.reserve(graph->GetNodes());
				for (TUNGraph::TNodeI NI = graph->BegNI(); NI != graph->EndNI(); NI++)
				{
					nodes.push_back(NI.GetId());
				}
				preTCTrussNodes.push_back(nodes);
			}
			else
			{
				for (auto& edge : newTCs[index - TCTrussSize])
				{
					int u = edge.first, v = edge.second;
					if (!result->IsNode(u))
					{
						result->AddNode(u);
						newNodes.push_back(u);
					}
					if (!result->IsNode(v))
					{
						result->AddNode(v);
						newNodes.push_back(v);
					}
					result->AddEdge(u, v);
				}
			}
		}
		AddTCTruss(result);

		unordered_set<int> oldNodeSet_1; // O_1
		unordered_map<int, int> nodeCntMap;
		for (auto& nodes : preTCTrussNodes)
		{
			for (auto node : nodes)
			{
				oldNodeSet_1.insert(node);
				nodeCntMap[node]++;
			}
		}

		// obtain P(N)
		int newNodeCnt = newNodes.size();
		for (int i = 0; i < newNodeCnt; i++)
		{
			for (int j = i + 1; j < newNodeCnt; j++)
			{
				int u = newNodes[i], v = newNodes[j];
				if (u < v)
				{
					newPairs[u].push_back(v);
				}
				else
				{
					newPairs[v].push_back(u);
				}
			}
		}

		// obtain {(u, v)|u ¡Ê N, v ¡Ê O_1}
		for (int i = 0; i < newNodeCnt; i++)
		{
			int u = newNodes[i];
			for (auto& v : oldNodeSet_1)
			{
				if (u < v)
				{
					newPairs[u].push_back(v);
				}
				else
				{
					newPairs[v].push_back(u);
				}
			}
		}

		// obtain P(O_1) \ \mathcal{P}(\mathcal{T}_O)
		for (auto& oldNodes : preTCTrussNodes)
		{
			vector<int> diffNodes; // D
			diffNodes.reserve(oldNodes.size());
			unordered_set<int> oldNodeSet_2(oldNodeSet_1); // O_2
			for (auto& oldNode : oldNodes)
			{
				nodeCntMap[oldNode]--;
				if (nodeCntMap[oldNode] == 0)
				{
					diffNodes.push_back(oldNode);
					oldNodeSet_1.erase(oldNode);
				}
				oldNodeSet_2.erase(oldNode);
			}

			for (auto& u : diffNodes)
			{
				for (auto& v : oldNodeSet_2)
				{
					if (u < v)
					{
						newPairs[u].push_back(v);
					}
					else
					{
						newPairs[v].push_back(u);
					}
				}
			}
		}
	}
}

flat_hash_map<int, vector<int>>& DynamicTCTruss::GetNewPairs()
{
	return newPairs;
}

void DynamicTCTruss::ResetNewPairs()
{
	flat_hash_map<int, vector<int>>().swap(newPairs);
	newPairs = flat_hash_map<int, vector<int>>();
}

PUNGraph DynamicTCTruss::GetTruss()
{
	return T;
}
