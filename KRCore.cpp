/*
Copyright (c) of 2018 by Fan Zhang <fanzhang.cs@gmail.com>
*/

#include "KRCore.h"

KRCore::KRCore()
{
	verSet0Index.resize(MxNId);
	verAttriIndex.resize(MxNId);
	subgraphIndex.resize(MxNId);
}

void KRCore::makeSubgraphTage(long tagi, long tag)//find subgraphs of k-core
{
	long neiID, neitagi;
	for (long i = 2; i < verSet0[tagi].size(); i++)
	{
		neiID = verSet0[tagi][i];
		neitagi = verSet0Index[neiID];
		if (wholeGraphTag[neitagi] == 1)
		{
			wholeGraphTag[neitagi] = tag;
			makeSubgraphTage(neitagi, tag);
		}
	}
}

double KRCore::computeVertexSimilarity(const long vertexID1, const long vertexID2) //rely on attribute order
{
	auto& uAttrs = verAttri[verAttriIndex[vertexID1]], vAttrs = verAttri[verAttriIndex[vertexID2]];
	Counter c;
	set_intersection(uAttrs.begin(), uAttrs.end(), vAttrs.begin(), vAttrs.end(), back_inserter(c));

	// normalized the weight into [0, 1]
	int numerator = c.count * (nodeComAttrCnt.at(vertexID1) + nodeComAttrCnt.at(vertexID2));
	int denominator = (uAttrs.size() + vAttrs.size() - c.count) * 2 * Q.size();

	double simi = (double)numerator / denominator;
	simi = -simi;
	if (simi < -1)
	{
		printf("similarity computation wrong!\n");
	}
	return simi;
}

bool KRCore::checkComputeMaxDistance(const long tag, const vector<long> &setTag)
{
	double curDistance;
	for (long i = 0; i < setTag.size(); i++)
	{
		if (setTag[i] == tag)
		{
			for (long j = i + 1; j < setTag.size(); j++)
			{
				if (setTag[j] == tag)
				{
					curDistance = computeVertexSimilarity(verSet[i][0], verSet[j][0]);
					if (curDistance>maxCommunityDistance)
					{
						return 0;
					}
				}
			}
		}
	}
	return 1;
}

double KRCore::checkMaxDistance(const vector<long> &com, bool tag)
{
	for (long i = 0; i < com.size(); i++)
	{
		if (com[i])
		{
			for (long j = 0; j < subgraphDistTag[i].size(); j++)
			{
				long id = subgraphDistTag[i][j];
				if (com[id])
				{
					return 0;
				}
			}
			if (tag)
			{

			}
		}
	}
	return 1;
}

double KRCore::checkMaxDistance1(const vector<long> &com, const long id)
{
	for (long i = 0; i < subgraphDistTag[id].size(); i++)
	{
		long sid = subgraphDistTag[id][i];
		if (com[sid])
		{
			return 0;
		}
	}
	return 1;
}

bool KRCore::computeKcoreU(vector<long> &C, vector<long> &M, vector<long> &numOL, vector<long> &degreeMC, long long &scoreOL, long long &scoreMC)
{
	for (long i = 0; i < kcore.size(); i++)
	{
		if (C[i])
		{
			kcore[i] = 1;
		}
		else
		{
			kcore[i] = 0;
		}
	}
	for (long i = 0; i < kcore.size(); i++)
	{
		if (M[i])
		{
			kcore[i] = 1;
		}
	}

	//update degree
	setTagIDs.clear();

	for (long j = 0; j < kcore.size(); j++)
	{
		if (kcore[j])
		{
			long degree = degreeMC[j];
			subgraphSet[j][1] = degree;

			if (degree < inputK)
			{
				if (M[j])
				{
					return 0;
				}
				setTagIDs.push_back(j);
				kcore[j] = 2;
			}
		}
	}

	for (long i = 0; i < setTagIDs.size(); i++)
	{
		long id = setTagIDs[i];

		C[id] = 0;
		scoreOL += numOL[id];
		scoreMC += degreeMC[id];
		//	numD++;

		for (long j = 0; j < subgraphDistTag[id].size(); j++)
		{
			long jid = subgraphDistTag[id][j];
			if (C[jid])
			{
				numOL[jid]--;
			}
		}
		for (long j = 2; j < subgraphSet[id].size(); j++)
		{
			long jid = subgraphIndex[subgraphSet[id][j]];
			degreeMC[jid]--;
		}

		kcore[id] = 0; //cur vertex deleted
		subgraphSet[id][1] = 0;
		for (long j = 2; j < subgraphSet[id].size(); j++)
		{
			long neiVertexID = subgraphSet[id][j];
			long neiverSetID = subgraphIndex[neiVertexID];
			if (kcore[neiverSetID]) //subgraphSet.size()>neiverSetID && subgraphSet[neiverSetID][0] == neiVertexID &&
			{
				subgraphSet[neiverSetID][1]--; // neighbor degree - 1
				if (kcore[neiverSetID] == 1 && subgraphSet[neiverSetID][1] < inputK) //new candidate for computing
				{
					if (M[neiverSetID])
					{
						return 0;
					}
					setTagIDs.push_back(neiverSetID);
					kcore[neiverSetID] = 2;
				}
			}
		}
	}
	return 1;
}

bool KRCore::computeKcore(vector<long> &reC, vector<long> &reE0, vector<long> &E, vector<long> &C, vector<long> &M, const long tag)
{
	for (long i = 0; i < kcore.size(); i++)
	{
		if (C[i])
		{
			kcore[i] = 1;
		}
		else
		{
			kcore[i] = 0;
		}
	}
	for (long i = 0; i < kcore.size(); i++)
	{
		if (M[i])
		{
			kcore[i] = 1;
		}
	}
	//update degree
	setTagIDs.clear();

	if (tag)
	{
		for (long j = 0; j < kcore.size(); j++)
		{
			if (kcore[j])
			{
				long degree = degreeInMC[j];
				subgraphSet[j][1] = degree;

				if (degree < inputK)
				{
					if (M[j])
					{
						return 0;
					}
					setTagIDs.push_back(j);
					kcore[j] = 2;
				}
			}
		}
	}
	else
	{

		for (long j = 0; j < kcore.size(); j++)
		{
			if (kcore[j])
			{
				subgraphSet[j][1] = degreeMECM[j];
				if (subgraphSet[j][1] < inputK)
				{
					if (M[j])
					{
						return 0;
					}
					setTagIDs.push_back(j);
					kcore[j] = 2;
				}
			}
		}
	}


	for (long i = 0; i < setTagIDs.size(); i++)
	{
		long setTagID = setTagIDs[i];

		C[setTagID] = 0;

		reC.push_back(setTagID);


		if (tag)
		{
			long degreeE = 0;
			for (long k = 2; k < subgraphSet[setTagID].size(); k++)
			{
				long id = subgraphSet[setTagID][k];
				long sid = subgraphIndex[id];
				if (C[sid])
				{
					edgesInC--;
					degreeInMC[sid]--;
					degreeE++;
				}
				if (M[sid])
				{
					degreeInMC[sid]--;
					degreeE++;
				}
				if (E[sid])
				{
					degreeInMC[sid]--;
				}

			}
			numberInC--;
			degreeInMC[setTagID] = degreeE;

			E[setTagID] = 1;
			numberInE++;
			reE0.push_back(setTagID);

			for (long j = 0; j < subgraphDistTag[setTagID].size(); j++)
			{
				long id1 = subgraphDistTag[setTagID][j];
				if (C[id1])
				{
					numberViolateInC--;
					numberOverLong[id1]--;
				}
				vnumberInME[id1]++;
			}
			for (long k = 2; k < subgraphSet[setTagID].size(); k++)
			{
				long nid = subgraphSet[setTagID][k];
				long sid = subgraphIndex[nid];
				if (M[sid])
				{
					numberMC--;
				}
			}
		}
		else
		{
			for (long j = 2; j < subgraphSet[setTagID].size(); j++)
			{
				long id = subgraphSet[setTagID][j];
				long sid = subgraphIndex[id];
				degreeMECM[sid]--;
			}
			for (long j = 0; j < subgraphDistTag[setTagID].size(); j++)
			{
				long id = subgraphDistTag[setTagID][j];
				vnumberMECM[id]--;
			}
		}

		kcore[setTagID] = 0; //cur vertex deleted
		subgraphSet[setTagID][1] = 0;
		for (long j = 2; j < subgraphSet[setTagID].size(); j++)
		{
			long neiVertexID = subgraphSet[setTagID][j];
			long neiverSetID = subgraphIndex[neiVertexID];
			if (kcore[neiverSetID]) //subgraphSet.size()>neiverSetID && subgraphSet[neiverSetID][0] == neiVertexID &&
			{
				subgraphSet[neiverSetID][1]--; // neighbor degree - 1
				if (kcore[neiverSetID] == 1 && subgraphSet[neiverSetID][1] < inputK) //new candidate for computing
				{
					if (M[neiverSetID])
					{
						return 0;
					}
					setTagIDs.push_back(neiverSetID);
					kcore[neiverSetID] = 2;
				}
			}
		}
	}
	return 1;
}

bool KRCore::computeKcoreInsertion(vector<long> &C, vector<long> &M)
{
	for (long i = 0; i < kcore.size(); i++)
	{
		if (C[i])
		{
			kcore[i] = 1;
		}
		else
		{
			kcore[i] = 0;
		}
	}
	for (long i = 0; i < kcore.size(); i++)
	{
		if (M[i])
		{
			kcore[i] = 1;
		}
	}

	//update degree
	setTagIDs.clear();

	for (long j = 0; j < kcore.size(); j++)
	{
		if (kcore[j])
		{
			long degree = 0;
			for (long k = 2; k < subgraphSet[j].size(); k++)
			{
				long nid = subgraphSet[j][k];
				long sid = subgraphIndex[nid];
				if (kcore[sid])
				{
					degree++;
				}
			}
			subgraphSet[j][1] = degree;
			if (degree < inputK)
			{
				if (M[j])
				{
					return 0;
				}
				setTagIDs.push_back(j);
				kcore[j] = 2;
			}
		}
	}

	for (long i = 0; i < setTagIDs.size(); i++)
	{
		long setTagID = setTagIDs[i];

		C[setTagID] = 0;
		kcore[setTagID] = 0; //cur vertex deleted
		subgraphSet[setTagID][1] = 0;
		for (long j = 2; j < subgraphSet[setTagID].size(); j++)
		{
			long neiVertexID = subgraphSet[setTagID][j];
			long neiverSetID = subgraphIndex[neiVertexID];
			if (kcore[neiverSetID]) //subgraphSet.size()>neiverSetID && subgraphSet[neiverSetID][0] == neiVertexID &&
			{
				subgraphSet[neiverSetID][1]--; // neighbor degree - 1
				if (kcore[neiverSetID] == 1 && subgraphSet[neiverSetID][1] < inputK) //new candidate for computing
				{
					if (M[neiverSetID])
					{
						return 0;
					}
					setTagIDs.push_back(neiverSetID);
					kcore[neiverSetID] = 2;
				}
			}
		}
	}
	return 1;
}

void KRCore::computeKcore0()
{
	for (long i = 0; i < setTagIDs.size(); i++)
	{
		long setTagID = setTagIDs[i];
		wholeGraphTag[setTagID] = 0; //cur vertex deleted
		verSet0[setTagID][1] = 0;
		for (long j = 2; j < verSet0[setTagID].size(); j++)
		{
			long neiVertexID = verSet0[setTagID][j];
			long neiverSetID = verSet0Index[neiVertexID];
			if (verSet0[neiverSetID][0] == neiVertexID && wholeGraphTag[neiverSetID] != 0)
			{
				verSet0[neiverSetID][1]--; // neighbor degree - 1
				if (wholeGraphTag[neiverSetID] == 1 && verSet0[neiverSetID][1] < inputK) //new candidate for computing
				{
					setTagIDs.push_back(neiverSetID);
					wholeGraphTag[neiverSetID] = 2;
				}
			}
		}
	}
}

bool KRCore::checkMaximalDoor(long num)
{
	if (!num)
	{
		return 1;
	}

	recursionDepthCM = 1;
	maximal = 0;
	bool result = 0;
	vector<long> useless;
	useless.clear();

	computeKcore(useless, useless, useless, checkC, checkM, 0); //M must in Kcore
	num -= useless.size();

	if (num > 0)
	{
		//build degreeInMCE and vnumberInME
		result = checkMaximal(-1);
	}
	else
	{
		return 1;
	}

	return result;
}

bool KRCore::updateCM(vector<long> &reC, vector<long> &reM, const long u)
{
	if (u >= 0)
	{
		for (long i = 0; i < subgraphDistTag[u].size(); i++)
		{
			long id = subgraphDistTag[u][i];
			if (checkC[id])
			{
				checkC[id] = 0;
				reC.push_back(id);
				for (long j = 2; j < subgraphSet[id].size(); j++)
				{
					long vid = subgraphSet[id][j];
					long sid = subgraphIndex[vid];
					degreeMECM[sid]--;
				}
				for (long j = 0; j < subgraphDistTag[id].size(); j++)
				{
					long vid = subgraphDistTag[id][j];
					vnumberMECM[vid]--;
				}
			}
		}
		//u-degree
		long uDegree = 0;
		for (long i = 2; i < subgraphSet[u].size(); i++)
		{
			long nid = subgraphSet[u][i];
			long sid = subgraphIndex[nid];
			if (checkM[sid] || checkC[sid])
			{
				uDegree++;
			}
		}
		if (uDegree < inputK)
		{
			return 0;
		}
	}

	if (computeKcore(reC, reM, reM, checkC, checkM, 0))
	{
		return 1;
	}
	else
	{
		return 0;
	}

}

long KRCore::uSelectionInCheckMaximal()
{
	long u = -1;
	long size = checkC.size();
	long maxDegree = -1;
	for (long i = 0; i < checkC.size(); i++)
	{
		if (checkC[i] && vnumberInME[i] > 0)
		{
			long degree = degreeMECM[i];
			if (degree > maxDegree)
			{
				maxDegree = degree;
				u = i;
			}
		}
	}
	return u;
}

bool KRCore::checkMaximal(const long lastID)
{
	vector<long> recoverC, recoverM, recoverE;
	recoverC.clear(); recoverM.clear(); recoverE.clear();
	vector<long> degreeMECMstore, vnumberMECMstore;
	if (updateCM(recoverC, recoverM, lastID))
	{
		if (checkMaxDistance(checkC, 0))
		{
			maximal = 1;
			return 0;
		}
		else
		{
			long u = uSelectionInCheckMaximal();

			if (u >= 0)//expansion
			{	//
				checkC[u] = 0;
				checkM[u] = 1;

				recursionDepthCM++;
				checkMaximal(u);
				recursionDepthCM--;
				if (maximal)
				{
					return 0;
				}
				//
				checkM[u] = 0;

				for (long k = 2; k < subgraphSet[u].size(); k++)
				{
					long id = subgraphSet[u][k];
					long sid = subgraphIndex[id];
					degreeMECM[sid]--;
				}
				for (long k = 0; k < subgraphDistTag[u].size(); k++)
				{
					long id = subgraphDistTag[u][k];
					vnumberMECM[id]--;
				}
				recursionDepthCM++;
				checkMaximal(-1);
				recursionDepthCM--;
				if (maximal)
				{
					return 0;
				}
				//
				checkC[u] = 1;
				/*
				subgraphE[u] = 0;
				if (subgraphE2[u])
				{
				subgraphE2[u] = 0;
				numberInE2--;
				}
				*/
			}
		}
	}
	degreeMECM = degreeMECMstore;
	vnumberMECM = vnumberMECMstore;
	//recover C,M,E
	for (long i = 0; i<recoverC.size(); i++)
	{
		long id = recoverC[i];
		checkC[id] = 1;
	}
	for (long i = 0; i<recoverM.size(); i++)
	{
		long id = recoverM[i];
		checkM[id] = 0;
	}
	return 1;
}

void KRCore::adaptiveSearchDoor(long tag)
{
	for (long i = 0; i < coreSubgraphTag.size(); i++)
	{
		vector<long> setInsertion;
		vector<long> tagInsertion;
		vector<long> idInsertion;
		vector<double> distInsertion;
		vector<long> distTableInsertion;
		long id = 0;
		long firstID = i;
		long subgraphID;
		long vertexID = verSet[i][0];
		subgraphC.clear(); subgraphM.clear(); subgraphE.clear(); subgraphE2.clear();
		subgraphSet.clear(); subgraphDistTag.clear();
		checkC.clear(); checkM.clear(); kcore.clear(); numberOverLong.clear(); degreeInMC.clear();
		degreeInM.clear(); degreeInMCE.clear(); vnumberInME.clear();
		numberInC = 0;
		numberInE = 0;
		for (long j = 0; j < coreSubgraphTag.size(); j++)
		{
			long vertexID1 = verSet[j][0];
			double dist = computeVertexSimilarity(vertexID, vertexID1);

			if (dist <= maxCommunityDistance)
			{
				if (j < i)
				{
					subgraphC.push_back(0);
					subgraphM.push_back(0);
					subgraphE.push_back(1);
					numberInE++;
				}
				else if (j == i)
				{
					//printf(" %ld\n", vertexID1);
					subgraphC.push_back(0);
					subgraphM.push_back(1);
					subgraphE.push_back(0);
					numberInM = 1;
					subgraphID = id;
				}
				else
				{
					subgraphC.push_back(1);
					subgraphM.push_back(0);
					subgraphE.push_back(0);
					numberInC++;
				}
				subgraphE2.push_back(0);
				setInsertion.clear(); setInsertion.push_back(vertexID1);
				subgraphSet.push_back(setInsertion); subgraphIndex[vertexID1] = id;
				checkC.push_back(0); checkM.push_back(0);
				kcore.push_back(0);
				id++;
			}
		}

		if (id >= inputK)
		{
			//build set
			edgesInC = 0;
			for (long k = 0; k < subgraphSet.size(); k++)
			{
				long vertexID1 = subgraphSet[k][0];
				long verSet0ID = verSet0Index[vertexID1];
				bool containM = 0;
				long degreeC = 0;
				long degreeE = 0;
				subgraphSet[k].push_back(0);
				for (long l = 2; l < verSet0[verSet0ID].size(); l++)
				{
					long nID = verSet0[verSet0ID][l];
					long subID = subgraphIndex[nID];
					if (subgraphSet.size()>subID && subgraphSet[subID][0] == nID)
					{
						if (subID == subgraphID)
						{
							containM = 1;
						}
						subgraphSet[k].push_back(nID);

						if (subgraphC[subID])
						{
							degreeC++;
							if (subgraphC[k])
							{
								edgesInC++;
							}
						}
						if (subgraphE[subID])
						{
							degreeE++;
						}
					}

				}
				degreeInM.push_back(containM);
				degreeInMCE.push_back(containM + degreeE + degreeC);
				subgraphSet[k][1] = subgraphSet[k].size() - 2;
				degreeInMC.push_back(degreeC + containM);

			}
			edgesInC /= 2;
			//mc number
			numberMC = 0;
			for (long k = 2; k < subgraphSet[subgraphID].size(); k++)
			{
				long sid = subgraphIndex[subgraphSet[subgraphID][k]];
				if (subgraphC[sid])
				{
					numberMC++;
				}
			}

			numberViolateInC = 0;
			for (long k = 0; k < subgraphC.size(); k++)
			{
				long vertexID1 = subgraphSet[k][0];
				distInsertion.clear();
				tagInsertion.clear();
				idInsertion.clear();
				distTableInsertion.clear();
				long id = 0;
				long lid = 0;
				long meNum = 0;
				for (long l = 0; l < subgraphC.size(); l++)
				{
					long vertexID2 = subgraphSet[l][0];
					double distance = computeVertexSimilarity(vertexID1, vertexID2);
					if (distance > maxCommunityDistance)
					{
						distInsertion.push_back(distance);
						distTableInsertion.push_back(0);
						tagInsertion.push_back(l);
						idInsertion.push_back(id);
						id++;
						if (k > subgraphID && l > subgraphID)
						{
							lid++;
							numberViolateInC++;
						}
						if (l < subgraphID)
						{
							meNum++;
						}
					}
					else
					{
						distTableInsertion.push_back(1);
					}
				}
				numberOverLong.push_back(lid);
				vnumberInME.push_back(meNum);
				subgraphDistTag.push_back(tagInsertion);
			}
			numberViolateInC /= 2;

			//build E2
			numberInE2 = 0;
			for (long k = 0; k < subgraphE.size(); k++)
			{
				if (subgraphE[k] && checkMaxDistance1(subgraphC, k))
				{
					subgraphE2[k] = 1;
					numberInE2++;
				}
			}
			recursionDepth = 0;

			adaptiveSearch(-1);

			if (progTimeout)
			{
				return;
			}
		}
	}
}

bool KRCore::update(vector<long> &reC, vector<long> &reM, vector<long> &reE0, vector<long> &reE, vector<long> &reE2, const long u)
{
	if (u >= 0)
	{
		for (long i = 0; i < subgraphDistTag[u].size(); i++)
		{
			long id = subgraphDistTag[u][i];
			if (subgraphC[id])
			{
				subgraphC[id] = 0;
				long degreeE = 0;
				for (long j = 2; j < subgraphSet[id].size(); j++)
				{
					long vid = subgraphSet[id][j];
					long sid = subgraphIndex[vid];
					if (subgraphC[sid])
					{
						edgesInC--;
						degreeInMC[sid]--;
						degreeE++;
					}
					if (subgraphM[sid])
					{
						degreeInMC[sid]--;
						degreeE++;
					}
					if (subgraphE[sid])
					{
						degreeInMC[sid]--;
					}
					degreeInMCE[sid]--;
				}
				//degreeInMC[id] = degreeE;

				reC.push_back(id);
				numberInC--;
				for (long j = 0; j < subgraphDistTag[id].size(); j++)
				{
					long id1 = subgraphDistTag[id][j];
					if (subgraphC[id1])
					{
						numberOverLong[id1]--;
						numberViolateInC--;
					}
					//vnumberInME[id1]--;
				}
				for (long k = 2; k < subgraphSet[id].size(); k++)
				{
					long nid = subgraphSet[id][k];
					long sid = subgraphIndex[nid];
					if (subgraphM[sid])
					{
						numberMC--;
					}
				}
			}
			if (subgraphE[id])
			{
				subgraphE[id] = 0;
				for (long j = 2; j < subgraphSet[id].size(); j++)
				{
					long vid = subgraphSet[id][j];
					long sid = subgraphIndex[vid];
					degreeInMCE[sid]--;
				}
				for (long j = 0; j < subgraphDistTag[id].size(); j++)
				{
					long id1 = subgraphDistTag[id][j];
					vnumberInME[id1]--;
				}
				reE.push_back(id);
				numberInE--;

			}
			if (subgraphE2[id])
			{
				subgraphE2[id] = 0;
				reE2.push_back(id);
				numberInE2--;
			}
		}
	}

	if (computeKcore(reC, reE0, subgraphE, subgraphC, subgraphM, 1))
	{
		if (u >= 0 && (numberInM + numberInC >= inputK))
		{
			for (long i = 0; i < subgraphC.size(); i++)
			{
				if (subgraphE2[i] && degreeInMC[i] >= inputK)
				{
					return 0;
				}
			}
		}
		//early insertion
		if (1) //need update
		{
			//long reCsize = reC.size();
			vector<long> subC;
			subC = subgraphC;
			long num = 0;
			for (long i = 0; i < subgraphC.size(); i++)
			{
				if (subgraphC[i] && numberOverLong[i])
				{
					subC[i] = 0;
					num++;
				}
			}
			if (num && (num + numberInM >= inputK) && computeKcoreInsertion(subC, subgraphM))
			{
				for (long i = 0; i < subC.size(); i++)
				{
					if (subC[i])
					{
						subgraphC[i] = 0;
						subgraphM[i] = 1;
						reC.push_back(i);
						reM.push_back(i);
						for (long k = 2; k < subgraphSet[i].size(); k++)
						{
							long id = subgraphSet[i][k];
							long sid = subgraphIndex[id];
							if (subgraphC[sid])
							{
								edgesInC--;
								numberMC++;
								degreeInM[sid]++;
							}
							if (subgraphM[sid])
							{
								numberMC--;
							}

						}
						numberInM++;
						numberInC--;
						for (long k = 0; k < subgraphDistTag[i].size(); k++)
						{
							long id = subgraphDistTag[i][k];
							if (subgraphC[id])
							{
								numberOverLong[id]--;
								numberViolateInC--;
							}
							if (subgraphE[id])
							{
								subgraphE[id] = 0;
								for (long l = 0; l < subgraphDistTag[id].size(); l++)
								{
									long sid = subgraphDistTag[id][l];
									vnumberInME[sid]--;
								}
								for (long l = 2; l < subgraphSet[id].size(); l++)
								{
									long vid = subgraphSet[id][l];
									long sid = subgraphIndex[vid];
									degreeInMCE[sid]--;
								}
								reE.push_back(id);
								numberInE--;

							}
							if (subgraphE2[id])
							{
								subgraphE2[id] = 0;
								reE2.push_back(id);
								numberInE2--;
							}
							vnumberInME[id]++;
						}

					}
				}
			}
		}
		return 1;
	}
	else
	{
		return 0;
	}

}

long KRCore::uSelection()
{
	long u = -1;
	long size = subgraphC.size();
	long long maxNum = -size;
	long long maxNum1 = -size;
	for (long i = 0; i < size; i++)
	{
		if (subgraphC[i])
		{
			vector<long> subC = subgraphC;
			vector<long> subM = subgraphM;
			vector<long> numOL = numberOverLong;
			vector<long> degreeMC = degreeInMC;
			long long scoreOL1 = numberOverLong[i];
			long long scoreMC1 = 0;//
			subC[i] = 0;
			subM[i] = 1;
			for (long j = 0; j < subgraphDistTag[i].size(); j++)
			{
				long id = subgraphDistTag[i][j];
				if (subC[id])
				{
					numOL[id]--;
				}
			}
			for (long j = 0; j<subgraphDistTag[i].size(); j++)
			{
				long id = subgraphDistTag[i][j];
				if (subC[id])
				{
					subC[id] = 0;
					scoreOL1 += numOL[id];
					scoreMC1 += degreeMC[id];
					for (long k = 0; k < subgraphDistTag[id].size(); k++)
					{
						long kid = subgraphDistTag[id][k];
						if (subC[kid])
						{
							numOL[kid]--;
						}
					}
					for (long k = 2; k < subgraphSet[id].size(); k++)
					{
						long kid = subgraphIndex[subgraphSet[id][k]];
						degreeMC[kid]--;
					}
				}
			}
			bool goon1 = computeKcoreU(subC, subM, numOL, degreeMC, scoreOL1, scoreMC1);

			subC = subgraphC;
			subM = subgraphM;
			numOL = numberOverLong;
			degreeMC = degreeInMC;
			long long scoreOL2 = numberOverLong[i];
			long long scoreMC2 = degreeInMC[i];//
			subC[i] = 0;
			for (long j = 0; j < subgraphDistTag[i].size(); j++)
			{
				long id = subgraphDistTag[i][j];
				if (subC[id])
				{
					numOL[id]--;
				}
			}
			for (long j = 2; j < subgraphSet[i].size(); j++)
			{
				long id = subgraphIndex[subgraphSet[i][j]];
				degreeMC[id]--;
			}
			bool goon2 = computeKcoreU(subC, subM, numOL, degreeMC, scoreOL2, scoreMC2);

			long long score = (scoreOL1 + scoreOL2);
			long long score1 = (scoreMC1 + scoreMC2);
			if (score > maxNum)
			{
				maxNum = score;
				maxNum1 = score1;
				u = i;
			}
			else if (score == maxNum && score1 < maxNum1)
			{
				maxNum1 = score1;
				u = i;
			}
		}
	}
	return u;
}

void KRCore::adaptiveSearch(const long lastID)
{
	double progTime = (double)clock() / CLOCKS_PER_SEC - algStartTime;
	if (progTime > TimeLimit)
	{
		progTimeout = 1;
		return;
	}

	vector<long> recoverC, recoverM, recoverE0, recoverE, recoverE2;
	recoverC.clear(); recoverM.clear(); recoverE0.clear(); recoverE.clear(); recoverE2.clear();
	vector<long> numberLongStore = numberOverLong;
	vector<long> degreeInCStore = degreeInMC;
	vector<long> vnumberInMEStore = vnumberInME;
	vector<long> degreeInMEStore = degreeInMCE;
	long numVio = numberViolateInC;
	long numMC = numberMC;

	if (update(recoverC, recoverM, recoverE0, recoverE, recoverE2, lastID))
	{
		bool goon = 0;
		if (numberInM >= 0 && checkMaxDistance(subgraphC, 1))
		{
			goon = 1;
		}

		if (goon)
		{
			if (1)
			{
				long checkCNum = 0;
				bool checkResult = 1;
				//build checkM
				for (long i = 0; i < checkM.size(); i++)
				{
					if (subgraphC[i] || subgraphM[i])
					{
						checkM[i] = 1;
					}
					else
					{
						checkM[i] = 0;
					}
				}
				goon = 1;

				degreeMECM = degreeInMCE;
				vnumberMECM = vnumberInME;

				//build checkC
				for (long i = 0; i < checkC.size(); i++)
				{
					if (subgraphE[i] && checkMaxDistance1(subgraphC, i))
					{
						checkC[i] = 1;
						checkCNum++;
						if (1)
						{
							if (degreeInMC[i] >= inputK)
							{
								goon = 0;
								break;
							}
						}
					}
					else
					{
						checkC[i] = 0;
						if (subgraphE[i])
						{
							for (long j = 0; j < subgraphDistTag[i].size(); j++)
							{
								long id = subgraphDistTag[i][j];
								vnumberMECM[id]--;
							}
							for (long j = 2; j < subgraphSet[i].size(); j++)
							{
								long vid = subgraphSet[i][j];
								long sid = subgraphIndex[vid];
								degreeMECM[sid]--;
							}
						}
					}
				}
				if (goon)
				{
					checkResult = checkMaximalDoor(checkCNum);

					if (checkResult)
					{
						long coreSize = 0;

						vector<long> temp;
						temp.clear();
						for (long k = 0; k < checkM.size(); k++)
						{
							if (checkM[k])
							{
								temp.push_back(subgraphSet[k][0]);
								coreSize++;
							}
						}
			
						krcoreStore.push_back(temp);
					}
				}
			}
			else
			{
				bool exist = 0;

				vector<long> core;
				core.clear();
				for (long i = 0; i < subgraphC.size(); i++)
				{
					if (subgraphC[i] || subgraphM[i])
					{
						core.push_back(subgraphSet[i][0]);
					}
				}
				krcoreStore.push_back(core);
			}
		}
		else
		{
			long u = uSelection();
			if (u >= 0)//expansion
			{
				subgraphC[u] = 0;
				subgraphM[u] = 1;
				for (long k = 2; k < subgraphSet[u].size(); k++)
				{
					long id = subgraphSet[u][k];
					long sid = subgraphIndex[id];
					if (subgraphC[sid])
					{
						edgesInC--;
					}
					//degreeInMCE[sid]++;
				}

				numberInM++;
				numberInC--;

				for (long k = 0; k < subgraphDistTag[u].size(); k++)
				{
					long id = subgraphDistTag[u][k];
					if (subgraphC[id])
					{
						numberOverLong[id]--;
						numberViolateInC--;
					}
					vnumberInME[id]++;
				}
				for (long k = 2; k < subgraphSet[u].size(); k++)
				{
					long nid = subgraphSet[u][k];
					long sid = subgraphIndex[nid];
					if (subgraphC[sid])
					{
						numberMC++;
						degreeInM[sid]++;
					}
					if (subgraphM[sid])
					{
						numberMC--;
					}
				}
	
				recursionDepth++;
				//C in R(u), no need to update numberOverLong
				adaptiveSearch(u);
				recursionDepth--;
				if (progTimeout)
				{
					return;
				}

				subgraphM[u] = 0;
				numberInM--;
				long degreeE = 0;
				for (long k = 2; k < subgraphSet[u].size(); k++)
				{
					long nid = subgraphSet[u][k];
					long sid = subgraphIndex[nid];
					if (subgraphC[sid])
					{
						numberMC--;
						degreeInM[sid]--;
						degreeInMC[sid]--;
						degreeE++;
					}
					if (subgraphM[sid])
					{
						degreeInMC[sid]--;
						degreeE++;
					}
					if (subgraphE[sid])
					{
						degreeInMC[sid]--;
					}
				}
				degreeInMC[u] = degreeE;
				subgraphE[u] = 1;
				numberInE++;

				if (checkMaxDistance1(subgraphC, u))
				{
					subgraphE2[u] = 1;
					numberInE2++;
				}
	
				recursionDepth++;
				adaptiveSearch(-1);
				recursionDepth--;

				if (progTimeout)
				{
					return;
				}

				subgraphC[u] = 1;
				for (long k = 2; k < subgraphSet[u].size(); k++)
				{
					long id = subgraphSet[u][k];
					long sid = subgraphIndex[id];
					if (subgraphC[sid])
					{
						edgesInC++;
					}
				}
				numberInC++;
				subgraphE[u] = 0;
				numberInE--;

				if (subgraphE2[u])
				{
					subgraphE2[u] = 0;
					numberInE2--;
				}
			}
		}
	}
	//recover C,M,E
	numberViolateInC = numVio;
	numberMC = numMC;
	numberOverLong = numberLongStore;
	degreeInMC = degreeInCStore;
	vnumberInME = vnumberInMEStore;
	degreeInMCE = degreeInMEStore;
	for (long i = 0; i<recoverC.size(); i++)
	{
		long id = recoverC[i];
		subgraphC[id] = 1;
		for (long k = 2; k < subgraphSet[id].size(); k++)
		{
			long vid = subgraphSet[id][k];
			long sid = subgraphIndex[vid];
			if (subgraphC[sid])
			{
				edgesInC++;
			}
		}
		numberInC++;

	}
	for (long i = 0; i<recoverM.size(); i++)
	{
		long id = recoverM[i];
		subgraphM[id] = 0; //!
		numberInM--;
		for (long k = 2; k < subgraphSet[id].size(); k++)
		{
			long nid = subgraphSet[id][k];
			long sid = subgraphIndex[nid];
			degreeInM[sid]--;
		}
	}
	for (long i = 0; i<recoverE.size(); i++)
	{
		long id = recoverE[i];
		subgraphE[id] = 1;
		numberInE++;
	}
	for (long i = 0; i<recoverE0.size(); i++)
	{
		long id = recoverE0[i];
		subgraphE[id] = 0;
		numberInE--;
	}
	for (long i = 0; i<recoverE2.size(); i++)
	{
		long id = recoverE2[i];
		subgraphE2[id] = 1;
		numberInE2++;
	}

}

void KRCore::dataInput(vector<Edge>& edges)
{
	//read edges, build verSet
	long vertexID, neighborID;
	double verSimilarity;
	vector<long> verSetInsertion;
	for (auto& edge: edges)
	{
		vertexID = edge.first;
		neighborID = edge.second;
		// verSimilarity = computeVertexSimilarity(vertexID, neighborID);
		// if (verSimilarity <= maxCommunityDistance) //edge pre-deletion
		if (true)
		{
			if (verSet0.size() == 0)
			{
				verSetInsertion.clear();
				verSetInsertion.push_back(vertexID);
				verSetInsertion.push_back(neighborID);
				verSet0.push_back(verSetInsertion);
			}
			else
			{
				if (vertexID == verSet0[verSet0.size() - 1][0])
				{
					verSet0[verSet0.size() - 1].push_back(neighborID);
				}
				else
				{
					verSetInsertion.clear();
					verSetInsertion.push_back(vertexID);
					verSetInsertion.push_back(neighborID);
					verSet0.push_back(verSetInsertion);
				}
			}
		}
	}

	//change the order of verSet
	sort(verSet0.begin(), verSet0.end(), CompareLongClique);

	long verSetDegree;
	for (long i = 0; i < verSet0.size(); i++)
	{
		wholeGraphTag.push_back(1);
	}

	//build verSet0
	for (long i = 0; i < verSet0.size(); i++)
	{
		verSetDegree = verSet0[i].size() - 1;
		if (verSetDegree < inputK)
		{
			setTagIDs.push_back(i);
			wholeGraphTag[i] = 2;
		}
		verSet0[i].insert(verSet0[i].begin() + 1, verSetDegree);
	}

	//Index verSet0
	for (long i = 0; i < verSet0.size(); i++)
	{
		verSet0Index[verSet0[i][0]] = i;
	}
}

void KRCore::initAttrs(const PUNGraph G, const flat_hash_map<int, vector<Attribute>>& attrsOfNode, const vector<int>& queryAttrs)
{
	int nodeCnt = G->GetNodes();
	verAttri.reserve(nodeCnt);

	Q = queryAttrs;
	sort(Q.begin(), Q.end());

	int i = 0;
	for (auto NI = G->BegNI(); NI != G->EndNI(); NI++, i++)
	{
		int u = NI.GetId();
		auto& uAttrs = attrsOfNode.at(u);

		Counter c;
		set_intersection(Q.begin(), Q.end(), uAttrs.begin(), uAttrs.end(), back_inserter(c));
		nodeComAttrCnt.emplace(u, c.count);

		vector<long> uAttrsL;
		uAttrsL.reserve(uAttrs.size());
		for (auto attr : uAttrs)
		{
			uAttrsL.push_back(attr);
		}
		verAttri.push_back(std::move(uAttrsL));
		verAttriIndex[u] = i;
	}
}

void KRCore::enumerateKRCore()
{
	algStartTime = (double)clock() / CLOCKS_PER_SEC;
	
	//compute k-core
	computeKcore0();

	//find k-core subgraph
	long tag = 1;
	for (long i = 0; i < wholeGraphTag.size(); i++)
	{
		if (wholeGraphTag[i] == 1)
		{
			tag++;
			wholeGraphTag[i] = tag;
			makeSubgraphTage(i, tag);
		}
	}

	//compute on each connected subgraph
	for (long i = 2; i <= tag; i++)
	{
		//build subgraph
		verSet.clear();
		coreSubgraphTag.clear();
		for (long j = 0; j<wholeGraphTag.size(); j++)
		{
			if (wholeGraphTag[j] == i)
			{
				coreSubgraphTag.push_back(1);
				vector<long> temp;
				temp.clear();
				temp.push_back(verSet0[j][0]);
				verSet.push_back(temp);
			}
		}

		if (!checkComputeMaxDistance(1, coreSubgraphTag))
		{
			//enumeration
			adaptiveSearchDoor(i);

			if (progTimeout)
			{
				break;
			}
		}
		else
		{
			vector<long> com;
			com.clear();
			long coreSize = 0;
			for (long j = 0; j < coreSubgraphTag.size(); j++)
			{
				com.push_back(verSet[j][0]);
				coreSize++;
			}
			krcoreStore.push_back(com);
		}
	}
	algEndTime = (double)clock() / CLOCKS_PER_SEC;
	algorithmTime = algEndTime - algStartTime;
}

vector<Clique> KRCore::getKRCore()
{
	// 把krCoreStore输出即可
	vector<Clique> results;
	for (auto& core : krcoreStore)
	{
		Clique clique;
		clique.reserve(core.size());
		for (auto node : core)
		{
			clique.push_back(node);
		}
		results.push_back(clique);
	}
	return results;
}

void KRCore::init(int k, double r)
{
	// 旧数据清理
	maximal = false;
	progTimeout = false;
	numberInM = 0; 
	numberInE2 = 0;
	numberInC = 0;  
	numberInE = 0;
	edgesInC = 0;
	numberViolateInC = 0;
	numberMC = 0;
	recursionDepth = 0;
	recursionDepthCM = 0;
	algStartTime = 0;
	algorithmTime = 0;
	algEndTime = 0;
	degreeInM.clear();
	degreeInMC.clear();
	distInC.clear();
	numberOverLong.clear();
	degreeInMCE.clear();
	vnumberInME.clear();
	degreeMECM.clear();
	vnumberMECM.clear();
	krcoreStore.clear();
	subgraphC.clear();
	subgraphM.clear();
	subgraphE.clear();
	subgraphE2.clear();
	kcore.clear();
	checkC.clear();
	checkM.clear();
	verSet.clear();
	verSet0.clear();
	wholeGraphTag.clear();
	coreSubgraphTag.clear();
	setTagIDs.clear();
	subgraphDistTag.clear();
	subgraphSet.clear();

	// fill(verSet0Index.begin(), verSet0Index.end(), 0);
	fill(subgraphIndex.begin(), subgraphIndex.end(), 0);

	inputK = k;
	inputR = r;
	maxCommunityDistance = -inputR; 
}
