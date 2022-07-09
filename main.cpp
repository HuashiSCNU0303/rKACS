// MaxMinWeightTruss.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "Def.h"
#include "MaxMinWeightTruss.h"
#include "SFormat.h"
#include "Test.h"

#if defined(_WIN64)
#include "Snap.cpp"
#endif

int DEFAULT_K = 0;
int DEFAULT_R = 0;

void QueryOnFile(string dataSet)
{
    string dataset = dataSet;

    // 读取数据
    MaxMinWeightTruss maxWeightTruss;
    maxWeightTruss.G = TSnap::LoadEdgeList<PUNGraph>(SFormat("./DataSets/{0}/graph.txt", dataset).data());
    maxWeightTruss.LoadTrussness(SFormat("./DataSets/{0}/edge_trussness.txt", dataset).data());
    maxWeightTruss.LoadVertexAttribute(SFormat("./DataSets/{0}/user_attributes_int.txt", dataset).data());
    maxWeightTruss.LoadSimilarity(SFormat("./DataSets/{0}/edge_similarity.bin", dataset).data());

    ifstream queryFile(SFormat("./DataSets/{0}/query.txt", dataset));
    int queryCnt;
    queryFile >> queryCnt;
    for (int i = 0; i < queryCnt; i++)
    {
        if (i % 2 == 0)
        {
            cout << "—————第" << i / 2 + 1 << "次—————" << endl;
        }
        ofstream resultFile("./Results/" + dataset + "/exact/" + to_string(i) + ".txt");
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
            results = maxWeightTruss.BottomUpQuery_Maintain(k, r, queryAttrs, IMCE);
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }
        else if (type == 2)
        {
            clock_t startTime = clock();
            results = maxWeightTruss.BottomUpQuery_KRCore(k, r, queryAttrs);
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }
        else if (type == 3)
        {
            clock_t startTime = clock();
            results = maxWeightTruss.BottomUpQuery_Pruning(k, r, queryAttrs); //需要的时候再改
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }
        else if (type == 4)
        {
            clock_t startTime = clock();
            results = maxWeightTruss.BottomUpQuery_Maintain(k, r, queryAttrs, BINSERT);
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }
        else if (type == 5)
        {
            clock_t startTime = clock();
            results = maxWeightTruss.BottomUpQuery_Baseline(k, r, queryAttrs);
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }
        else if (type == 6)
        {
            clock_t startTime = clock();
            results = maxWeightTruss.BottomUpQuery_Maintain(k, r, queryAttrs, BINSERTH);
            printf("查询耗时：%lf\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);
        }
        if (i % 2 == 1)
        {
            cout << "——————————" << endl;
        }

        for (auto weight_coms : results)
        {
            auto weight = weight_coms.first;
            resultFile << weight.numerator << "/" << weight.denominator << " (" << (double)weight.numerator / weight.denominator << "): ";
            sort(weight_coms.second.begin(), weight_coms.second.end(), CompareClique);
            for (auto com : weight_coms.second)
            {
                for (auto node : com)
                {
                    resultFile << node << " ";
                }
                resultFile << " // ";
            }
            resultFile << endl;
        }

    }
}

void RQuery_Time(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    // varying r, from {1, 5, 10, 25, 50, 75, 100}
    const string type = "r", type2 = "time";
    cout << "RQuery_Time() start..." << endl;

    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));
    ofstream avgTimeFile(SFormat(AVG_RESULT_FILE, dataSet, type2, type));

    string rValueStr;
    vector<int> rValues;
    queryFile >> rValueStr;
    SplitInt(rValueStr, rValues, ",");

    string methodValueStr;
    vector<int> methodValues;
    queryFile >> methodValueStr;
    SplitInt(methodValueStr, methodValues, ",");

    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    for (int m : methodValues)
    {
        flat_hash_set<int> failedIndexes;
        for (int r : rValues)
        {
            avgTimeFile << "r = " << r << " :" << endl;

            ofstream optFile(SFormat(RESULT_FILE, dataSet, type2, type, r, METHOD_NAMES[m - 1]));
            double avgTime = 0;
            int failedCnt = 0, successCnt = 0;
            for (int i = 0; i < queryCnt; i++)
            {
                vector<int> queryAttr = queryAttrs[i];
                Time time = failedIndexes.find(i) != failedIndexes.end() ? -1: 
                    maxWeightTruss.Query_TimeWrapper(m, DEFAULT_K, r, queryAttr);
                if (time == -1)
                {
                    failedCnt++;
                    failedIndexes.emplace(i);
                    time = TIME_LIMIT;
                }
                successCnt++;
                optFile << SFormat("time = {0}s", time) << endl;
                avgTime += time;
                if (failedCnt >= queryCnt)
                {
                    break;
                }
            }

            avgTime = avgTime / successCnt;
            avgTimeFile << SFormat("{0}, avg time = {1}s", METHOD_NAMES[m - 1], avgTime) << endl; 
        }
        avgTimeFile << "————————————" << endl;
    }
    cout << "RQuery_Time() finish..." << endl;
}

void KQuery_Time(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    const string type = "k", type2 = "time";
    cout << "KQuery_Time() start..." << endl;

    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));
    ofstream avgTimeFile(SFormat(AVG_RESULT_FILE, dataSet, type2, type));
    
    string kValueStr;
    vector<int> kValues;
    queryFile >> kValueStr;
    SplitInt(kValueStr, kValues, ",");

    string methodValueStr;
    vector<int> methodValues;
    queryFile >> methodValueStr;
    SplitInt(methodValueStr, methodValues, ",");
    
    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    for (int m : methodValues)
    {
        for (int k : kValues)
        {
            avgTimeFile << "k = " << k << " :" << endl;

            ofstream optFile(SFormat(RESULT_FILE, dataSet, type2, type, k, METHOD_NAMES[m - 1]));
            double avgTime = 0;
            int successCnt = 0;
            for (int i = 0; i < queryCnt; i++)
            {
                vector<int> queryAttr = queryAttrs[i];
                Time time = maxWeightTruss.Query_TimeWrapper(m, k, DEFAULT_R, queryAttr);
                if (time == -1)
                {
                    time = TIME_LIMIT;
                }
                successCnt++;
                optFile << SFormat("time = {0}s", time) << endl;
                avgTime += time;
            }

            avgTime = avgTime / successCnt;
            avgTimeFile << SFormat("{0}, avg time = {1}s", METHOD_NAMES[m - 1], avgTime) << endl;
        }
        avgTimeFile << "————————————" << endl;
    }
    cout << "KQuery_Time() finish..." << endl;
}

void AttrCntQuery_Time(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    // varying the size of |queryAttrs|, from {1, 3, 5, 7, 9, 11, 13}
    const string type = "attrCnt", type2 = "time";
    cout << "AttrCntQuery_Time() start..." << endl;
    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));
    ofstream avgTimeFile(SFormat(AVG_RESULT_FILE, dataSet, type2, type));

    string attrCntStr;
    vector<int> attrCnts;
    queryFile >> attrCntStr;
    SplitInt(attrCntStr, attrCnts, ",");

    string methodValueStr;
    vector<int> methodValues;
    queryFile >> methodValueStr;
    SplitInt(methodValueStr, methodValues, ",");

    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    for (int m : methodValues)
    {
        for (int attrCnt : attrCnts)
        {
            avgTimeFile << "attrCnt = " << attrCnt << " :" << endl;

            ofstream optFile(SFormat(RESULT_FILE, dataSet, type2, type, attrCnt, METHOD_NAMES[m - 1]));
            double avgTime = 0;
            int successCnt = 0;
            for (int i = 0; i < queryCnt; i++)
            {
                vector<int> attrs = queryAttrs[i], queryAttr;
                for (int j = 0; j < attrCnt; j++)
                {
                    queryAttr.push_back(attrs[j]);
                }
                Time time = maxWeightTruss.Query_TimeWrapper(m, DEFAULT_K, DEFAULT_R, queryAttr);
                if (time == -1)
                {
                    time = TIME_LIMIT;
                }
                successCnt++;
                optFile << SFormat("time = {0}s", time) << endl;
                avgTime += time;
            }

            avgTime = avgTime / successCnt;
            avgTimeFile << SFormat("{0}, avg time = {1}s", METHOD_NAMES[m - 1], avgTime) << endl;
        }
        avgTimeFile << "————————————" << endl;
    }
    cout << "AttrCntQuery_Time() finish..." << endl;
}

void RQuery_Memory(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    // varying r, from {1, 5, 10, 25, 50, 75, 100}
    const string type = "r", type2 = "memory";
    cout << "RQuery_Memory() start..." << endl;

    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));
    ofstream avgMemoryFile(SFormat(AVG_RESULT_FILE, dataSet, type2, type));

    string rValueStr;
    vector<int> rValues;
    queryFile >> rValueStr;
    SplitInt(rValueStr, rValues, ",");

    string methodValueStr;
    vector<int> methodValues;
    queryFile >> methodValueStr;
    SplitInt(methodValueStr, methodValues, ",");

    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;
    vector<Memory> queryMemories;

    queryMemories.reserve(queryCnt);

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    for (int m : methodValues)
    {
        for (int r : rValues)
        {
            avgMemoryFile << "r = " << r << " :" << endl;

            ofstream optFile(SFormat(RESULT_FILE, dataSet, type2, type, r, METHOD_NAMES[m - 1]));
            for (int i = 0; i < queryCnt; i++)
            {
                vector<int> queryAttr = queryAttrs[i];
                Memory memory = maxWeightTruss.Query_MemoryWrapper(m, DEFAULT_K, r, queryAttr);
                queryMemories.push_back(memory);
            }

            double avgMemory = 0;
            for (auto memory : queryMemories)
            {
                optFile << "memory = " << memory << "s" << endl;
                avgMemory += memory;
            }
            avgMemory = (double)avgMemory / queryCnt;
            avgMemoryFile << METHOD_NAMES[m - 1] << ", avg memory = " << avgMemory << "s" << endl;

            queryMemories.clear();
        }
        avgMemoryFile << "————————————" << endl;
    }
    cout << "RQuery_Memory() finish..." << endl;
}

void KQuery_Memory(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    const string type = "k", type2 = "memory";
    cout << "KQuery_Memory() start..." << endl;

    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));
    ofstream avgMemoryFile(SFormat(AVG_RESULT_FILE, dataSet, type2, type));

    string kValueStr;
    vector<int> kValues;
    queryFile >> kValueStr;
    SplitInt(kValueStr, kValues, ",");

    string methodValueStr;
    vector<int> methodValues;
    queryFile >> methodValueStr;
    SplitInt(methodValueStr, methodValues, ",");

    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;
    vector<Memory> queryMemories;

    queryMemories.reserve(queryCnt);

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    for (int m : methodValues)
    {
        for (int k : kValues)
        {
            avgMemoryFile << "k = " << k << " :" << endl;

            ofstream optFile(SFormat(RESULT_FILE, dataSet, type2, type, k, METHOD_NAMES[m - 1]));
            for (int i = 0; i < queryCnt; i++)
            {
                vector<int> queryAttr = queryAttrs[i];
                Memory memory = maxWeightTruss.Query_MemoryWrapper(m, k, DEFAULT_R, queryAttr);
                queryMemories.push_back(memory);
            }

            double avgMemory = 0;
            for (auto memory : queryMemories)
            {
                optFile << "memory = " << memory << "kb" << endl;
                avgMemory += memory;
            }
            avgMemory = (double)avgMemory / queryCnt;
            avgMemoryFile << METHOD_NAMES[m - 1] << ", avg memory = " << avgMemory << "s" << endl;

            queryMemories.clear();
        }
        avgMemoryFile << "————————————" << endl;
    }
    cout << "KQuery_Memory() finish..." << endl;
}

void AttrCntQuery_Memory(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    // varying the size of |queryAttrs|, from {1, 3, 5, 7, 9, 11, 13}
    const string type = "attrCnt", type2 = "memory";
    cout << "AttrCntQuery_Memory() start..." << endl;
    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));
    ofstream avgMemoryFile(SFormat(AVG_RESULT_FILE, dataSet, type2, type));

    string attrCntStr;
    vector<int> attrCnts;
    queryFile >> attrCntStr;
    SplitInt(attrCntStr, attrCnts, ",");

    string methodValueStr;
    vector<int> methodValues;
    queryFile >> methodValueStr;
    SplitInt(methodValueStr, methodValues, ",");

    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;
    vector<Memory> queryMemories;

    queryMemories.reserve(queryCnt);

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    for (int m : methodValues)
    {
        for (int attrCnt : attrCnts)
        {
            avgMemoryFile << "attrCnt = " << attrCnt << " :" << endl;

            ofstream optFile(SFormat(RESULT_FILE, dataSet, type2, type, attrCnt, METHOD_NAMES[m - 1]));
            for (int i = 0; i < queryCnt; i++)
            {
                vector<int> attrs = queryAttrs[i], queryAttr;
                for (int j = 0; j < attrCnt; j++)
                {
                    queryAttr.push_back(attrs[j]);
                }
                Memory memory = maxWeightTruss.Query_MemoryWrapper(m, DEFAULT_K, DEFAULT_R, queryAttr);
                queryMemories.push_back(memory);
            }

            double avgMemory = 0;
            for (auto memory : queryMemories)
            {
                optFile << "memory = " << memory << "s" << endl;
                avgMemory += memory;
            }
            avgMemory = (double)avgMemory / queryCnt;
            avgMemoryFile << METHOD_NAMES[m - 1] << ", avg memory = " << avgMemory << "s" << endl;

            queryMemories.clear();
        }
        avgMemoryFile << "————————————" << endl;
    }
    cout << "AttrCntQuery_Memory() finish..." << endl;
}

void Query_Quality(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    const string type = "quality";
    cout << "Query_Quality() start..." << endl;

    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));

    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    string methods[3] = { "kc","pruning", "vac" };

    for (int m = 0; m < 3; m++)
    {
        ofstream resultFile(SFormat(QUALITY_RESULT_FILE, dataSet, methods[m]));
        for (int i = 0; i < queryCnt; i++)
        {
            vector<int> queryAttr_ = queryAttrs[i];
            ResultMap results;
            if (m == 2)
            {
                results = maxWeightTruss.VACQuery(DEFAULT_K, queryAttr_);
            }
            else if (m == 1)
            {
                results = maxWeightTruss.BottomUpQuery_Pruning(DEFAULT_K, DEFAULT_R, queryAttr_);
            }
            else
            {
                results = maxWeightTruss.KCQuery(DEFAULT_K, queryAttr_);
            }
            string str = "";
            if (!results.empty())
            {
                auto& resultComs = results.begin()->second;
                Clique com = *resultComs.begin();
                for (auto node : com)
                {
                    str.append(to_string(node) + ",");
                }
                str = str.substr(0, str.length() - 1);
            }
            resultFile << str << endl;
            /*if (!results.empty() && results.begin()->first.numerator != INT_MAX)
            {
                for (auto node : queryAttr_)
                {
                    str.append(to_string(node) + ",");
                }
                str = str.substr(0, str.length() - 1);
                resultFile << str << endl;
            }
            if (i % 100 == 0)
            {
                cout << i << endl;
            }*/
        }
    }

    cout << "Query_Quality() finish..." << endl;
}

void QualityEval(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    const string type = "quality";
    cout << "QualityEval() start..." << endl;

    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));

    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    string queryResultPath = "./Results/{0}/quality/{1}_quality_result.txt";
    string resultMetricDetailPath = "./Results/{0}/quality/{1}_result_metrics.txt";
    string methods[4] = { "pruning", "vac",  "kc", "concs" };
    vector<vector<vector<int>>> queryResults;
    for (int i = 0; i < 4; i++)
    {
        ifstream resultFile(SFormat(queryResultPath.c_str(), dataSet, methods[i]));
        vector<vector<int>> results;
        for (int i = 0; i < queryCnt; i++)
        {
            string resultStr;
            vector<int> result;
            resultFile >> resultStr;
            if (!resultStr.empty())
            {
                SplitInt(resultStr, result, ",");
                results.push_back(result);
            }
        }
        queryResults.push_back(results);
    }

    for (int m = 0; m < 4; m++)
    {
        string method = methods[m];
        auto& results = queryResults[m];
        double totalCQS = 0, totalAPJ = 0, totalAKS = 0, totalACSF = 0, totalECOQUG = 0;
        double totalCSS = 0, totalCAS = 0, totalAScore = 0, totalDensity = 0, totalConDensity = 0, totalKwdC = 0;
        ofstream metricFile(SFormat(resultMetricDetailPath.c_str(), dataSet, method));
        for (int i = 0; i < queryCnt; i++)
        {
            vector<int> queryAttr = queryAttrs[i];
            sort(queryAttr.begin(), queryAttr.end());
            double CQS = maxWeightTruss.ComputeScore2(queryAttr, results[i]);
            double APJ = maxWeightTruss.ComputeScore3(queryAttr, results[i]);
            double AKS = maxWeightTruss.ComputeScore4(queryAttr, results[i]);
            double ACSF = maxWeightTruss.ComputeACSF(queryAttr, results[i]);
            double ECOQUG = maxWeightTruss.ComputeECOQUG(results[i]);
            totalCQS += CQS;
            totalAPJ += APJ;
            totalAKS += AKS;
            totalACSF += ACSF;
            totalECOQUG += ECOQUG;

            double css = maxWeightTruss.ComputeCSS(results[i]);
            double cas = maxWeightTruss.ComputeCAS(results[i]);
            double aScore = maxWeightTruss.ComputeASCore(results[i]);
            double density = maxWeightTruss.ComputeDensity(results[i]);
            double conDensity = maxWeightTruss.ComputeConDensity(queryAttr, results[i]);
            double kwdC = maxWeightTruss.ComputeKwdC(queryAttr, results[i]);

            totalCAS += cas;
            totalCSS += css;
            totalAScore += aScore;
            totalDensity += density;
            totalConDensity += conDensity;
            totalKwdC += kwdC;

            metricFile << SFormat("cqs = {0}, apj = {1}, aks = {2}, css = {3}, cas = {4}, aScore = {5} density = {6}, conDensity = {7}, kwdC = {8}, acsf = {9}, ecoqug = {10}\n", CQS, APJ, AKS, css, cas, aScore, density, conDensity, kwdC, ACSF, ECOQUG);
        }
        printf("%s, cqs = %f, apj = %f, aks = %f, acsf = %f, ecoqug = %f, css = %f, cas = %f, aScore = %f, density = %f, conDensity = %f, kwdC = %f\n", method, totalCQS / queryCnt, totalAPJ / queryCnt,
            totalAKS / queryCnt, totalACSF / queryCnt, totalECOQUG / queryCnt,
            totalCSS / queryCnt, totalCAS / queryCnt, totalAScore / queryCnt, totalDensity / queryCnt, totalConDensity / queryCnt, totalKwdC / queryCnt);
    }

    cout << "QualityEval() finish..." << endl;
}

void DataArrange(PUNGraph G, flat_hash_map<int, vector<Attribute>>& attributes)
{
    for (auto& node_attrs : attributes)
    {
        auto& attrs = node_attrs.second;
        sort(attrs.begin(), attrs.end());
    }

    for (auto NI = G->BegNI(); NI != G->EndNI(); NI++)
    {
        int u = NI.GetId();
        if (attributes.find(u) == attributes.end())
        {
            Clique clique;
            clique.push_back(1500000);
            attributes.emplace(u, clique);
        }
    }

    map<int, Clique> sortedAttributes;
    for (auto& node_attrs : attributes)
    {
        sortedAttributes.emplace(node_attrs.first, node_attrs.second);
    }

    string filePath = "./DataSets/2/user_attributes_int_new.txt";
    ofstream resultFile(filePath);
    for (auto& node_attrs : sortedAttributes)
    {
        string str = "";
        for (auto node : node_attrs.second)
        {
            str.append(to_string(node) + ",");
        }
        str = str.substr(0, str.length() - 1);
        resultFile << node_attrs.first << "\t" << str << endl;
    }
}

void GenQuery(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    vector<Attribute> attrs;
    auto& attrNodeMap = maxWeightTruss.nodesInAttr;
    attrs.reserve(attrNodeMap.size());
    for (auto& attr_nodes : attrNodeMap)
    {
        attrs.push_back(attr_nodes.first);
    }
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(attrs.begin(), attrs.end(), std::default_random_engine(seed));

    ofstream attrFile(SFormat("./DataSets/{0}/qualified_attr.txt", dataSet));

    vector<int> resultAttrs;
    while (true)
    {
        Attribute attr = attrs.back();
        attrs.pop_back();

        vector<int> queryAttrs;
        queryAttrs.push_back(attr);

        Time time = maxWeightTruss.Query_TimeWrapper(BINSERTH, DEFAULT_K, DEFAULT_R, queryAttrs);
        if (time > 0.001)
        {
            resultAttrs.push_back(attr);
            attrFile << SFormat("{0}: {1} nodes\n", attr, 1);
        }

        if (resultAttrs.size() > 1000)
        {
            break;
        }
    }
}

void DataSetStats(MaxMinWeightTruss& maxWeightTruss)
{
    auto& graph = maxWeightTruss.G;
    int nodeCnt = graph->GetNodes(), edgeCnt = graph->GetEdges();
    
    int totalDeg = 0, totalAttrCnt = 0;
    for (auto NI = graph->BegNI(); NI != graph->EndNI(); NI++)
    {
        totalDeg += NI.GetDeg();
        totalAttrCnt += maxWeightTruss.attrsOfNode[NI.GetId()].size();
    }

    int totalAttrNodeCnt = 0, attrCnt = maxWeightTruss.nodesInAttr.size();
    for (auto& attr_nodes : maxWeightTruss.nodesInAttr)
    {
        totalAttrNodeCnt += attr_nodes.second.size();
    }

    cout << SFormat("nodeCnt = {0}, edgeCnt = {1}, avgDeg = {2}, avgAttrCnt = {3}, avgAttrNodeCnt = {4}", 
        nodeCnt, edgeCnt, (double)totalDeg / nodeCnt, (double)totalAttrCnt / nodeCnt, (double)totalAttrNodeCnt / attrCnt);
}

void GenSmallGraph(MaxMinWeightTruss& maxWeightTruss, string dataSet)
{
    const string type = "k", type2 = "time";

    ifstream queryFile(SFormat(QUERY_FILE, dataSet, type));

    /*string attrCntStr;
    vector<int> attrCnts;
    queryFile >> attrCntStr;
    SplitInt(attrCntStr, attrCnts, ",");*/

    int queryCnt;
    queryFile >> queryCnt;
    vector<vector<int>> queryAttrs;

    for (int i = 0; i < queryCnt; i++)
    {
        string queryAttrStr;
        vector<int> queryAttr;
        queryFile >> queryAttrStr;
        SplitInt(queryAttrStr, queryAttr, ",");
        queryAttrs.push_back(queryAttr);
    }

    /*int index = 0;
    for (int i = 0; i < queryCnt; i++)
    {
        if (index >= 100)
        {
            break;
        }
        vector<int> attrs = queryAttrs[i], queryAttr;
        for (int j = 0; j < attrCnts[0]; j++)
        {
            queryAttr.push_back(attrs[j]);
        }
        auto results = maxWeightTruss.BottomUpQuery_Pruning(DEFAULT_K, DEFAULT_R, queryAttr);
        if (!results.empty() && maxWeightTruss.maxKTruss->GetNodes() > 10)
        {
            for (auto attrCnt : attrCnts)
            {
                vector<int> queryAttr;
                for (int j = 0; j < attrCnt; j++)
                {
                    queryAttr.push_back(attrs[j]);
                }
                maxWeightTruss.BottomUpQuery_Pruning(DEFAULT_K, DEFAULT_R, queryAttr);
                ofstream smallGraphFile(SFormat("./DataSets/{0}/SmallGraph/{1}/{2}.edge.txt", dataSet, attrCnt, index));
                ofstream qQueryFile(SFormat("./DataSets/{0}/SmallGraph/{1}/query.txt", dataSet, attrCnt), ios::app);
                auto& smallGraph = maxWeightTruss.maxKTruss;
                smallGraphFile << SFormat("0 {0}\n", smallGraph->GetEdges() * 2);
                for (auto EI = smallGraph->BegEI(); EI != smallGraph->EndEI(); EI++)
                {
                    int u = EI.GetSrcNId(), v = EI.GetDstNId();
                    smallGraphFile << SFormat("{0} {1}\n{1} {0}\n", u, v);
                }
                qQueryFile << queryAttr.size() << " ";
                for (auto attr : queryAttr)
                {
                    qQueryFile << attr << " ";
                }
                qQueryFile << endl;  
            }
            index++;
        }
    }*/
    int index = 0;
    ofstream qQueryFile(SFormat("./DataSets/{0}/SmallGraph/qQuery.txt", dataSet));
    for (int i = 0; i < queryCnt; i++)
    {
        if (index >= 100)
        {
            break;
        }
        vector<int> queryAttr = queryAttrs[i];
        auto results = maxWeightTruss.BottomUpQuery_Pruning(DEFAULT_K, DEFAULT_R, queryAttr);
        if (!results.empty())
        {
            maxWeightTruss.BottomUpQuery_Pruning(DEFAULT_K - 4, DEFAULT_R, queryAttr);
            ofstream smallGraphFile(SFormat("./DataSets/{0}/SmallGraph/{1}.edge.txt", dataSet, index));
            auto& smallGraph = maxWeightTruss.maxKTruss;
            if (smallGraph->GetEdges() < 15 || smallGraph->GetEdges() > 1000)
            {
                continue;
            }
            smallGraphFile << SFormat("0 {0}\n", smallGraph->GetEdges() * 2);
            for (auto EI = smallGraph->BegEI(); EI != smallGraph->EndEI(); EI++)
            {
                int u = EI.GetSrcNId(), v = EI.GetDstNId();
                smallGraphFile << SFormat("{0} {1}\n{1} {0}\n", u, v);
            }
            qQueryFile << queryAttr.size() << " ";
            for (auto attr : queryAttr)
            {
                qQueryFile << attr << " ";
            }
            qQueryFile << endl;
            index++;
        }
    }
}

void Query(int argc, char** argv)
{
    string dataSet = string(argv[1]);
    string type = string(argv[2]);
    string queries = string(argv[3]);
    DEFAULT_K = atoi(argv[4]);
    DEFAULT_R = atoi(argv[5]);

    MaxMinWeightTruss maxWeightTruss;
    maxWeightTruss.G = TSnap::LoadEdgeList<PUNGraph>(SFormat(GRAPH_FILE, dataSet).c_str());
    cout << "Load graph complete..." << endl;
    // maxWeightTruss.G->SortNodeAdjV();
    // TSnap::SaveEdgeList(maxWeightTruss.G, SFormat(GRAPH_FILE, dataSet).c_str());
    maxWeightTruss.LoadTrussness(SFormat(TRUSSNESS_FILE, dataSet).c_str());
    cout << "Load trussness complete..." << endl;
    maxWeightTruss.LoadVertexAttribute(SFormat(ATTRIBUTE_FILE, dataSet).c_str());
    cout << "Load attributes complete..." << endl;
    maxWeightTruss.LoadSimilarity(SFormat(SIMILARITY_FILE, dataSet).c_str());
    cout << "Load similarity complete..." << endl;

    if (type == "t") // time
    {
        for (int i = 0; i < queries.size(); i++)
        {
            int indicator = queries[i];
            if (indicator == 'r')
            {
                RQuery_Time(maxWeightTruss, dataSet);
            }
            else if (indicator == 'k')
            {
                KQuery_Time(maxWeightTruss, dataSet);
            }
            else if (indicator == 'c')
            {
                AttrCntQuery_Time(maxWeightTruss, dataSet);
            }
        }
    }
    else if (type == "m") // memory
    {
        for (int i = 0; i < queries.size(); i++)
        {
            int indicator = queries[i];
            if (indicator == 'r')
            {
                RQuery_Memory(maxWeightTruss, dataSet);
            }
            else if (indicator == 'k')
            {
                KQuery_Memory(maxWeightTruss, dataSet);
            }
            else if (indicator == 'c')
            {
                AttrCntQuery_Memory(maxWeightTruss, dataSet);
            }
        }
    }
    else if (type == "q") // quality
    {
        Query_Quality(maxWeightTruss, dataSet);
        QualityEval(maxWeightTruss, dataSet);
    }
    else if (type == "qe")
    {
        QualityEval(maxWeightTruss, dataSet);
    }
    else if (type == "g")
    {
        // GenQuery(maxWeightTruss, dataSet);
        GenSmallGraph(maxWeightTruss, dataSet);
    }
    else if (type == "s")
    {
        DataSetStats(maxWeightTruss);
    }
    else if (type == "vac")
    {
        Clique queryAttrs;
        auto results = maxWeightTruss.VACQuery(DEFAULT_K, queryAttrs);
        auto& coms = results.begin()->second;
        for (auto& com : coms)
        {
            auto weight = results.begin()->first;
            cout << (double)weight.numerator / weight.denominator << endl;
            auto it = find(com.begin(), com.end(), 170);
            if (it != com.end())
            {
                com.erase(it);
            }
            double score = maxWeightTruss.ComputeASCore(com);
            cout << score << endl;
            for (auto node : com)
            {
                cout << node << " ";
            }
        }
    }
}

int main(int argc, char **argv)
{
    // Query(argc, argv);
    QueryOnFile(string(argv[1]));
    
    
    // MurmurHashTest();
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
