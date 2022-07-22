#pragma once

// the (k, r)-core implementation, 
// slightly modified based on the implementation of Fan Zhang
// https://github.com/fanzhangcs/krCoreEnum/

#include "Def.h"

constexpr auto TimeLimit = 3600;
constexpr auto MxNId = 10000000;

class KRCore
{
	bool maximal = 0, progTimeout = 0;
	long inputK; double inputR;
	long numberInM, numberInE2, numberInC, numberInE, edgesInC;
	long numberViolateInC, numberMC;
	long recursionDepth = 0, recursionDepthCM = 0;
	double algStartTime, algorithmTime, algEndTime;
	double maxCommunityDistance;
	vector<long> degreeInM, degreeInMC, distInC, numberOverLong, degreeInMCE, vnumberInME, degreeMECM, vnumberMECM;
	vector<long> verAttriIndex, verSet0Index, subgraphIndex;
	vector<long> wholeGraphTag, coreSubgraphTag;
	vector<long> subgraphC, subgraphM, subgraphE, subgraphE2;
	vector<long> kcore, checkC, checkM, setTagIDs;
	vector<vector<long> > verSet, verSet0, verAttri, subgraphSet, subgraphDistTag, krcoreStore;
	vector<int> Q; // query attributes
	flat_hash_map<int, int> nodeComAttrCnt;

	void makeSubgraphTage(long tagi, long tag);
	double computeVertexSimilarity(const long vertexID1, const long vertexID2);
	bool checkComputeMaxDistance(const long tag, const vector<long>& setTag);
	double checkMaxDistance(const vector<long>& com, bool tag);
	double checkMaxDistance1(const vector<long>& com, const long id);
	bool computeKcoreU(vector<long>& C, vector<long>& M, vector<long>& numOL, vector<long>& degreeMC, long long& scoreOL, long long& scoreMC);
	bool computeKcore(vector<long>& reC, vector<long>& reE0, vector<long>& E, vector<long>& C, vector<long>& M, const long tag);
	bool computeKcoreInsertion(vector<long>& C, vector<long>& M);
	void computeKcore0();
	bool checkMaximalDoor(long num);
	bool updateCM(vector<long>& reC, vector<long>& reM, const long u);
	long uSelectionInCheckMaximal();
	bool checkMaximal(const long lastID);
	void adaptiveSearchDoor(long tag);
	bool update(vector<long>& reC, vector<long>& reM, vector<long>& reE0, vector<long>& reE, vector<long>& reE2, const long u);
	long uSelection(); 	
	void adaptiveSearch(const long lastID);

public:
	KRCore();
	void dataInput(vector<Edge>& edges);
	void enumerateKRCore();
	vector<Clique> getKRCore();
	void init(int k, double r);
	void initAttrs(const PUNGraph G, const flat_hash_map<int, vector<Attribute>>& attrsOfNode, const vector<int>& queryAttrs);
};

