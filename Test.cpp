#include "Test.h"

void PruningTestOnFile()
{
    // 读取数据
    MaxMinWeightTruss maxWeightTruss;
    maxWeightTruss.G = TSnap::LoadEdgeList<PUNGraph>("./DataSets/graph.txt");
    maxWeightTruss.LoadTrussness("./DataSets/edge_trussness.txt");
    maxWeightTruss.LoadVertexAttribute("./DataSets/author_attributes_int.txt");

    ifstream queryFile("./DataSets/pruning_test.txt");
    int queryCnt;
    queryFile >> queryCnt;
    for (int i = 0; i < queryCnt; i++)
    {
        // 输入k, r
        int k, type, r;
        string queryAttrStr;
        vector<int> queryAttrs;
        queryFile >> k >> r >> queryAttrStr >> type;
        SplitInt(queryAttrStr, queryAttrs, ",");

        ResultMap results;
        if (type == 1)
        {
            clock_t startTime = clock();
            results = maxWeightTruss.BottomUpQuery_Maintain(k, r, queryAttrs, BINSERT);
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }
        /*else if (type == 2)
        {
            clock_t startTime = clock();
            results = maxWeightTruss.OnlyPruning(k, r, queryAttrs);
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }
        else if (type == 3)
        {
            clock_t startTime = clock();
            results = maxWeightTruss.OnlyPruning_Cand1(k, r, queryAttrs);
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }*/
    }
}

void AddEdgesTest()
{
    // 三个文件，一个原图的文件；一个clique的文件，然后基于clique的文件构建cliquesOfNode；一个是要添加的边的图的文件
    ifstream file("./AddEdgesTest/cliques.txt");
    vector<Clique> cliques;
    string line;
    while (getline(file, line)) // line中不包括每行的换行符
    {
        Clique clique;
        SplitInt(line, clique, ",");
        cliques.push_back(clique);
    }

    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("./AddEdgesTest/graph.txt");
    PUNGraph edgesToInsert = TSnap::LoadEdgeList<PUNGraph>("./AddEdgesTest/newEdges.txt");

    DynamicClique* dynamicClique = new DynamicClique(G, cliques);
    dynamicClique->AddEdges_BInsert(edgesToInsert);

    cout << "new cliques: " << endl;
    for (auto& clique : dynamicClique->maxCliques)
    {
        for (auto node : clique)
        {
            cout << node << ", ";
        }
        cout << endl;
    }
}

void DelEdgesTest()
{
    // 三个文件，一个原图的文件；一个clique的文件，然后基于clique的文件构建cliquesOfNode；一个是要添加的边的图的文件
    ifstream file("./DelEdgesTest/cliques.txt");
    vector<Clique> cliques;
    string line;
    while (getline(file, line)) // line中不包括每行的换行符
    {
        Clique clique;
        SplitInt(line, clique, ",");
        cliques.push_back(clique);
    }

    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("./DelEdgesTest/graph.txt");
    PUNGraph edgesToDelete = TSnap::LoadEdgeList<PUNGraph>("./DelEdgesTest/newEdges.txt");

    DynamicClique* dynamicClique = new DynamicClique(G, cliques);
    dynamicClique->DelEdges(edgesToDelete);

    cout << "new cliques: " << endl;
    for (auto& clique : dynamicClique->maxCliques)
    {
        for (auto node : clique)
        {
            cout << node << ", ";
        }
        cout << endl;
    }
}

void AddEdgeTest()
{
    // 三个文件，一个原图的文件；一个clique的文件，然后基于clique的文件构建cliquesOfNode；一个是要添加的边的图的文件
    ifstream file("./Test/AddEdgesTest/cliques.txt");
    vector<Clique> cliques;
    string line;
    while (getline(file, line)) // line中不包括每行的换行符
    {
        Clique clique;
        SplitInt(line, clique, ",");
        cliques.push_back(clique);
    }

    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("./Test/AddEdgesTest/graph.txt");
    PUNGraph edgesToInsert = TSnap::LoadEdgeList<PUNGraph>("./Test/AddEdgesTest/newEdges.txt");

    DynamicClique* dynamicClique = new DynamicClique(G, cliques);

    for (TUNGraph::TEdgeI EI = edgesToInsert->BegEI(); EI != edgesToInsert->EndEI(); EI++)
    {
        int u = EI.GetSrcNId(), v = EI.GetDstNId();
        dynamicClique->AddEdge(u, v);
        cout << SFormat("add edge({0}, {1}) finished, new cliques: ", u, v) << endl;
        for (auto& clique : dynamicClique->maxCliques)
        {
            if (clique.empty())
            {
                continue;
            }
            for (auto node : clique)
            {
                cout << node << ", ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

void IMCETest()
{
    // 三个文件，一个原图的文件；一个clique的文件，然后基于clique的文件构建cliquesOfNode；一个是要添加的边的图的文件
    ifstream file("./Test/AddEdgesTest/cliques.txt");
    vector<Clique> cliques;
    string line;
    while (getline(file, line)) // line中不包括每行的换行符
    {
        Clique clique;
        SplitInt(line, clique, ",");
        cliques.push_back(clique);
    }

    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("./Test/AddEdgesTest/graph.txt");
    PUNGraph edgesToInsert = TSnap::LoadEdgeList<PUNGraph>("./Test/AddEdgesTest/newEdges.txt");

    DynamicClique* dynamicClique = new DynamicClique(G, cliques);
    dynamicClique->AddEdges_BInsert(edgesToInsert);
    for (auto& clique : dynamicClique->maxCliques)
    {
        if (clique.empty())
        {
            continue;
        }
        for (auto node : clique)
        {
            cout << node << ", ";
        }
        cout << endl;
    }
    /*for (TUNGraph::TEdgeI EI = edgesToInsert->BegEI(); EI != edgesToInsert->EndEI(); EI++)
    {
        int u = EI.GetSrcNId(), v = EI.GetDstNId();
        dynamicClique->AddEdge(u, v);
        cout << SFormat("add edge({0}, {1}) finished, new cliques: ", u, v) << endl;
        for (auto& clique : dynamicClique->maxCliques)
        {
            if (clique.empty())
            {
                continue;
            }
            for (auto node : clique)
            {
                cout << node << ", ";
            }
            cout << endl;
        }
        cout << endl;
    }*/
}

void TrieTest()
{
    Trie* trie = new Trie();
    ifstream file("./TrieTest/cliques.txt");
    vector<Clique> cliques;
    string line;
    int index = 0;
    while (getline(file, line)) // line中不包括每行的换行符
    {
        Clique clique;
        SplitInt(line, clique, ",");
        UpdateTrie(trie, clique, index);
        cliques.push_back(clique);
        index++;
    }

    Clique cand;
    cand.reserve(cliques.size());
    GetCandsOnTrie(trie, cand);

    DestroyTrie(trie);

    cout << "candidates: " << endl;
    for (auto& index : cand)
    {
        for (auto node : cliques[index])
        {
            cout << node << ", ";
        }
        cout << endl;
    }
}

void BinInsertTest()
{
    ifstream file;
    vector<int> nums;
    file.open("data.txt");
    int cnt;
    file >> cnt;
    for (int i = 0; i < cnt; i++)
    {
        int num;
        file >> num;
        BinInsertNode(nums, num);
        for (auto node : nums)
        {
            cout << node << " ";
        }
        cout << endl;
    }
}

void GetSubGraphTest()
{
    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("./GetSubGraphTest/graph.txt");
    ifstream file("./GetSubGraphTest/testcase.txt");
    string line;
    int i = 1;
    while (getline(file, line)) // line中不包括每行的换行符
    {
        vector<int> nodeSet;
        SplitInt(line, nodeSet, ",");

        vector<Edge> edges, correctEdges;
        // PUNGraph subGraph = GetSubGraph1(G, nodeSet);
        flat_hash_map<int, Clique> subGraph;
        subGraph = GetSubGraph1(G, nodeSet, subGraph);
        for (auto& n_nbr: subGraph)
        {
            int n = n_nbr.first;
            for (auto& nbr : n_nbr.second)
            {
                if (nbr > n)
                {
                    edges.push_back(GetEdge(n, nbr));

                    // 看有没有边(nbr, n)
                    auto& nNbrs = subGraph.at(nbr);
                    if (SearchNode(nNbrs, n) == -1)
                    {
                        cout << "子图有问题！" << endl;
                    }
                }
            }
        }
        sort(edges.begin(), edges.end(), SortEdge);

        PUNGraph correctSubGraph = TSnap::GetSubGraph(G, ToTIntV(nodeSet));
        correctEdges.reserve(correctSubGraph->GetEdges());
        for (TUNGraph::TEdgeI EI = correctSubGraph->BegEI(); EI != correctSubGraph->EndEI(); EI++)
        {
            correctEdges.push_back(GetEdge(EI.GetSrcNId(), EI.GetDstNId()));
        }
        sort(correctEdges.begin(), correctEdges.end(), SortEdge);

        if (edges != correctEdges)
        {
            cout << "第" << i << "个测试用例有问题！" << endl;
        }
        i++;
    }
}

void DeleteNodeTest()
{
    ifstream file("./DeleteNodeTest/testcase.txt");
    string line;
    int testCnt;
    file >> testCnt;
    for (int i = 0; i < testCnt; i++)
    {
        string setStr, delNodeStr;
        file >> setStr >> delNodeStr;
        vector<int> nodeSet, delNodes;
        SplitInt(setStr, nodeSet, ",");
        SplitInt(delNodeStr, delNodes, ",");

        vector<int> nodeSet_c = nodeSet;

        for (auto node : delNodes)
        {
            DeleteNode(nodeSet, node);

            int index = SearchNode(nodeSet_c, node);
            if (index != -1)
            {
                nodeSet_c.erase(nodeSet_c.begin() + index);
            }

            if (nodeSet != nodeSet_c)
            {
                cout << "第" << i + 1 << "个测试用例有问题！" << endl;
                break;
            }
        }
    }
}

void IsInTest()
{
    ifstream file("./Test/IsInTest/testcase.txt");
    string line;
    int testCnt;
    file >> testCnt;
    for (int i = 0; i < testCnt; i++)
    {
        string smallSetStr, bigSetStr;
        file >> smallSetStr >> bigSetStr;
        vector<int> smallSet, bigSet;
        if (smallSetStr != "-1")
        {
            SplitInt(smallSetStr, smallSet, ",");
        }
        if (bigSetStr != "-1")
        {
            SplitInt(bigSetStr, bigSet, ",");
        }

        bool actual = IsIn(smallSet, bigSet);
        
        Counter c;
        set_intersection(smallSet.begin(), smallSet.end(), bigSet.begin(), bigSet.end(), back_inserter(c));
        bool expected = (c.count == smallSet.size());

        if (actual != expected)
        {
            cout << "第" << i + 1 << "个测试用例有问题！" << endl;
        }
    }
}

void GetMaximalCliquesTest()
{
    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("./Test/GetMaximalCliquesTest/graph.txt");
    PUNGraph edgesToInsert = TSnap::LoadEdgeList<PUNGraph>("./Test/GetMaximalCliquesTest/newEdges.txt");

    cout << "未加边：" << endl;
    vector<Clique> cliques;
    // GetMaximalCliques(G, cliques);
    for (auto& clique : cliques)
    {
        for (auto& node : clique)
        {
            cout << node << " ";
        }
        cout << endl;
    }
    cout << endl;

    vector<Edge> edges;
    for (TUNGraph::TEdgeI EI = edgesToInsert->BegEI(); EI != edgesToInsert->EndEI(); EI++)
    {
        int u = EI.GetSrcNId(), v = EI.GetDstNId();
        edges.push_back(GetEdge(u, v));
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(edges.begin(), edges.end(), std::default_random_engine(seed));

    for (auto& e : edges)
    {
        int u = e.first, v = e.second;
        cliques.clear();
        G->AddEdge2(u, v);
        // GetMaximalCliques(G, cliques);
        cout << SFormat("加边({0}, {1}): ", u, v) << endl;
        for (auto& clique : cliques)
        {
            for (auto& node : clique)
            {
                cout << node << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

void MoveTest()
{
    unordered_map<int, Clique> nbrMap;
    for (int i = 0; i < 5; i++)
    {
        vector<int> values;
        values.assign(10, i);
        nbrMap.emplace(i, std::move(values));
    }

    vector<Clique> cliques;
    for (int i = 0; i < 5; i++)
    {
        cliques.push_back(std::move(nbrMap[i]));
    }

    nbrMap[0].push_back(1);
}

void EraseTest()
{
    unordered_set<int> hashset;
    for (int i = 0; i < 10; i++)
    {
        hashset.emplace(i);
    }

    for (auto num : hashset)
    {
        cout << num << " ";
    }
    cout << endl;

    for (auto it = hashset.begin(); it != hashset.end(); )
    {
        it = hashset.erase(it);
        for (auto num : hashset)
        {
            cout << num << " ";
        }
        cout << endl;
    }
}

uint64_t MurmurHash64_Correct(const void* key, int len, uint64_t seed)
{
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

    const uint64_t* data = (const uint64_t*)key;
    const uint64_t* end = data + (len / 8);

    while (data != end)
    {
        uint64_t k = *data++;

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const unsigned char* data2 = (const unsigned char*)data;

    switch (len & 7)
    {
    case 7: h ^= uint64_t(data2[6]) << 48;
    case 6: h ^= uint64_t(data2[5]) << 40;
    case 5: h ^= uint64_t(data2[4]) << 32;
    case 4: h ^= uint64_t(data2[3]) << 24;
    case 3: h ^= uint64_t(data2[2]) << 16;
    case 2: h ^= uint64_t(data2[1]) << 8;
    case 1: h ^= uint64_t(data2[0]);
        h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

void MurmurHashTest()
{
    ifstream file("./Test/MurmurHashTest/testcase.txt");
    string line;
    int testCnt;
    file >> testCnt;
    for (int i = 0; i < testCnt; i++)
    {
        string dataStr;
        file >> dataStr;
        vector<int> data;
        if (dataStr != "-1")
        {
            SplitInt(dataStr, data, ",");
        }
        auto actual = MurmurHash2(data);

        unsigned char* c = (unsigned char*)&(data[0]);
        auto expected = MurmurHash64_Correct(c, data.size() * 4, 0);

        cout << SFormat("actual = {0}, expected = {1} ", actual, expected);

        if (actual != expected)
        {
            cout << "错误！";
        }
        cout << endl;
    }
}