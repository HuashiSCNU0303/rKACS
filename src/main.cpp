#include "Def.h"
#include "MaxMinWeightTruss.h"
#include "SFormat.h"

#if defined(_WIN64)
#include "Snap.cpp"
#endif

int DEFAULT_K = 0;
int DEFAULT_R = 0;

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

    string methods[3] = { "kc", "kac", "vac" };

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
                results = maxWeightTruss.rKACS_Incremental(DEFAULT_K, DEFAULT_R, queryAttr_);
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
        }
    }

    cout << "Query_Quality() finish..." << endl;
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
    if (maxWeightTruss.LoadTrussness(SFormat(TRUSSNESS_FILE, dataSet).c_str()) == 0)
    {
        maxWeightTruss.LoadTrussness(SFormat(TRUSSNESS_FILE, dataSet).c_str());
    }
    cout << "Load trussness complete..." << endl;
    maxWeightTruss.LoadVertexAttribute(SFormat(ATTRIBUTE_FILE, dataSet).c_str());
    cout << "Load attributes complete..." << endl;
    if (maxWeightTruss.LoadSimilarity(SFormat(SIMILARITY_FILE, dataSet).c_str()) == 0)
    {
        maxWeightTruss.LoadSimilarity(SFormat(SIMILARITY_FILE, dataSet).c_str());
    }
    cout << "Load similarity complete..." << endl;

    if (type == "t") // time, ignore the output communities
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
    else if (type == "q") // quality, save the output communities to file
    {
        Query_Quality(maxWeightTruss, dataSet);
    }
}

int main(int argc, char **argv)
{
    Query(argc, argv);
}
