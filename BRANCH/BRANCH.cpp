#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <lemon/list_graph.h>
#include <lemon/bfs.h>
#include <lemon/matching.h>
using namespace std;
using namespace lemon;

//array A = transVec, B and C = contiVec

#define MAX 99999

static double IDENTITY = 0.95;
static double GAPPED_IDENTITY = 0.1;
static double HIGH_IDENTITY = 0.99;

vector<vector<int> > transVec;//key data structure #0

typedef struct structPointerCopy
{
	int pointer;
	int copy;
} PointerCopy;

typedef struct structBC
{
	char base;
	int B;
	int C;
	int isTrans;
	vector<PointerCopy> PC;
	int segVecP;
	int isGene;
	int isSecondary;
	int oriB;
} BC;

vector<vector<BC> > contiVec;//key data structure #1

typedef struct structSegment
{
	int sourceStart;
	int targetStart;
	int size;
	int isTrans;
	int cov;//coverage
	int cpt;//#compatible reads
	int ovl;//#overlapped reads
	int shared;
	vector<int> path;
} Segment;

vector<vector<Segment> > segVec;//key data structure #2

typedef struct structSegJunc
{
	int isTrans;
	int copy;
	double weight;
	int shared;
	int postponed;
	vector<int> path;
} SegJunc;

vector<vector<SegJunc> > segJuncVec;//key data structure #3

vector<vector<vector<int> > > allPathVec;//key data structure #4

typedef struct structPosition
{
	int sourceID;
	int targetID;
	int targetStart;
	int targetEnd;
	vector<Segment> seg;
	int fr;
} Position;

vector<vector<int> > maxMatch;

int PE = 1;

void parse(string buf, int & targetID, int & targetStart, int & targetEnd, int & targetGap, int & sourceID, int & sourceStart, int & sourceEnd, int & sourceGap, int & sourceSize, vector<Segment> & seg, int & fr)
{
	char targetIDBuf[20] = {'\0'}, targetStartBuf[20] = {'\0'}, targetEndBuf[20] = {'\0'}, targetGapBuf[20] = {'\0'}, sourceIDBuf[20] = {'\0'}, sourceStartBuf[20] = {'\0'}, sourceEndBuf[20] = {'\0'}, sourceGapBuf[20] = {'\0'}, sourceSizeBuf[20] = {'\0'}, blockBuf[20] = {'\0'};
	int item = 0, i, j0 = 0, j1 = 0, j2 = 0, j3 = 0, j4 = 0, j5 = 0, j6 = 0, j7 = 0, j8 = 0, j9 = 0, k, tag0 = 1, tag1 = 0, tag2 = 0, tag4 = 1, sp;
	Segment s;

	seg.clear();

	for(i = 0; i < buf.size(); i ++)
	{
		if(buf[i] == '	')
		{
			item ++;
			continue;
		}
		if(item == 13 && tag0 == 1)
		{
			targetIDBuf[j0 ++] = buf[i];
			if(buf[i] == ':')
				tag0 = 0;
		}
		if(item == 15)
			targetStartBuf[j1 ++] = buf[i];
		if(item == 16)
			targetEndBuf[j2 ++] = buf[i];
		if(item == 7)
			targetGapBuf[j3 ++] = buf[i];
		if(item == 9 && tag4 == 1)
		{
			sourceIDBuf[j4 ++] = buf[i];
			if(buf[i] == ':')
				tag4 = 0;
		}
		if(item == 11)
			sourceStartBuf[j5 ++] = buf[i];
		if(item == 12)
			sourceEndBuf[j6 ++] = buf[i];
		if(item == 5)
			sourceGapBuf[j7 ++] = buf[i];
		if(item == 10)
			sourceSizeBuf[j8 ++] = buf[i];
		if(item == 18)
		{
			if(buf[i] == ',')
			{
				s.sourceStart = s.targetStart = -1;
				s.size = atoi(blockBuf);
				seg.push_back(s);
				
				for(k = 0; k < j9; k ++)
					blockBuf[k] = '\0';
				j9 = 0;
			}
			else
				blockBuf[j9 ++] = buf[i];
		}
		if(item == 19)
		{
			if(tag1 == 0)
			{
				sp = 0;
				tag1 = 1;
			}

			if(buf[i] == ',')
			{
				seg[sp ++].sourceStart = atoi(blockBuf);
				for(k = 0; k < j9; k ++)
					blockBuf[k] = '\0';
				j9 = 0;
			}
			else
				blockBuf[j9 ++] = buf[i];
		}
		if(item == 20)
		{
			if(tag2 == 0)
			{
				sp = 0;
				tag2 = 1;
			}

			if(buf[i] == '\0')
				break;
			if(buf[i] == ',')
			{
				seg[sp ++].targetStart = atoi(blockBuf);
				for(k = 0; k < j9; k ++)
					blockBuf[k] = '\0';
				j9 = 0;
			}
			else
				blockBuf[j9 ++] = buf[i];
		}
		if(item == 8)
			fr = buf[i] == '+' ? 1 : 0;
	}
	targetID = atoi(targetIDBuf);
	targetStart = atoi(targetStartBuf);
	targetEnd = atoi(targetEndBuf);
	targetGap = atoi(targetGapBuf);
	sourceID = atoi(sourceIDBuf);
	sourceStart = atoi(sourceStartBuf);
	sourceEnd = atoi(sourceEndBuf);
	sourceGap = atoi(sourceGapBuf);
	sourceSize = atoi(sourceSizeBuf);
}

int parseBLAT(ifstream & a, vector<vector<Position> > & v, double threshold1, double threshold2, int batchStartID)
{
	int sourceIDBak, i, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, vp, fr, tag = 1;
	vector<Segment> seg;
	Position pos;
	vector<Position> posVec;
	string buf;

	v.clear();

	sourceIDBak = batchStartID;

	if(a.is_open())
	{
		while(a.good())
		{
			getline(a, buf);
			if(buf[0] == 0) break;

			parse(buf, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, seg, fr);
//			cout << targetID << ", " << targetStart << ", " << targetEnd << ", " << sourceID << ", " << sourceStart << ", " << sourceEnd << ", " << sourceGap << ", " << sourceSize << ", " << (*seg.begin()).sourceStart << ", " << (*seg.begin()).targetStart << ", " << (*seg.begin()).size << endl;

			if(sourceIDBak == sourceID && tag == 1) continue;
                        tag = 0;

			if((double)(sourceEnd - sourceStart - sourceGap) / sourceSize >= threshold1 && (double)(targetEnd - targetStart - targetGap) / (double)(targetEnd - targetStart) >= threshold2)
			{
				if(sourceIDBak != sourceID)
				{
					for(i = sourceIDBak + 1; i < sourceID; i ++)
					{
						v.push_back(posVec);
						pos.targetID = pos.targetStart = pos.targetEnd = pos.sourceID = -1;
						vp = v.size() - 1;
						v[vp].push_back(pos);
						if(batchStartID != 0 && i == batchStartID + 1000000)
							return 0;
					}
					v.push_back(posVec);
					pos.targetID = targetID;
					pos.targetStart = targetStart;
					pos.targetEnd = targetEnd;
					pos.seg = seg;
					pos.sourceID = sourceID;
					pos.fr = fr;
					vp = v.size() - 1;
					v[vp].push_back(pos);
					if(batchStartID != 0 && sourceID == batchStartID + 1000000)
						return 0;
				}
				else
				{
					pos.targetID = targetID;
					pos.targetStart = targetStart;
					pos.targetEnd = targetEnd;
					pos.seg = seg;
					pos.sourceID = sourceID;
					pos.fr = fr;
					vp = v.size() - 1;
					v[vp].push_back(pos);
				}
				sourceIDBak = sourceID;
			}
		}
		return 1;
	}
	else
	{
		cout << "CANNOT OPEN FILE! (ERROR 1)" << endl;
		exit(-1);
	}
}

int isDonor(char base1, char base2)
{
	if(base1 == 'G' && base2 == 'T' ||
	base1 == 'C' && base2 == 'T' ||
	base1 == 'G' && base2 == 'C' ||
	base1 == 'C' && base2 == 'T' ||
	base1 == 'A' && base2 == 'T' ||
	base1 == 'G' && base2 == 'T')
		return 1;
	else
		return 0;
}

int isAcceptor(char base1, char base2)
{
        if(base1 == 'A' && base2 == 'G' ||
	base1 == 'A' && base2 == 'C' ||
	base1 == 'A' && base2 == 'G' ||
	base1 == 'G' && base2 == 'C' ||
	base1 == 'A' && base2 == 'C' ||
	base1 == 'A' && base2 == 'T')
                return 1;
        else
                return 0;
}

int isDonorAcceptor(char base11, char base12, char base21, char base22)
{
	if(base11 == 'G' && base12 == 'T' && base21 == 'A' && base22 == 'G' ||
	base11 == 'C' && base12 == 'T' && base21 == 'A' && base22 == 'C' ||
	base11 == 'G' && base12 == 'C' && base21 == 'A' && base22 == 'G' ||
	base11 == 'C' && base12 == 'T' && base21 == 'G' && base22 == 'C' ||
	base11 == 'A' && base12 == 'T' && base21 == 'A' && base22 == 'C' ||
	base11 == 'G' && base12 == 'T' && base21 == 'A' && base22 == 'T')
		return 1;
	else
		return 0;
}

void adjustStartEnd(int contiID, vector<Segment> & seg)
{
	int sp, i, j;

	for(sp = 0; sp < seg.size() - 1; sp ++)
	{
con:
		if(seg[sp + 1].targetStart - (seg[sp].targetStart + seg[sp].size) < 10) continue;

		if(seg[sp].size >= 10 && seg[sp + 1].size >= 10 && (!isDonor(contiVec[contiID][seg[sp].targetStart + seg[sp].size].base, contiVec[contiID][seg[sp].targetStart + seg[sp].size + 1].base) || !isAcceptor(contiVec[contiID][seg[sp + 1].targetStart - 2].base, contiVec[contiID][seg[sp + 1].targetStart - 2 + 1].base)))
			for(i = seg[sp].targetStart + seg[sp].size - 3; i <= seg[sp].targetStart + seg[sp].size + 3; i ++)
				for(j = seg[sp + 1].targetStart - 2 - 3; j <= seg[sp + 1].targetStart - 2 + 3; j ++)
					if(isDonorAcceptor(contiVec[contiID][i].base, contiVec[contiID][i + 1].base, contiVec[contiID][j].base, contiVec[contiID][j + 1].base))
					{
						seg[sp].size = i - seg[sp].targetStart;
						seg[sp + 1].sourceStart = seg[sp + 1].sourceStart + (j + 2 - seg[sp + 1].targetStart);
						seg[sp + 1].size = seg[sp + 1].size - (j + 2 - seg[sp + 1].targetStart);
						seg[sp + 1].targetStart = j + 2;
						sp ++;
						if(sp < seg.size() - 1)
							goto con;
						else
							goto end;
					}
	}
end:;
}

void updateContiVec(int contiID, vector<Segment> seg, int transID, int fr, int threshSplit)
{
	int tag, i, j, sp, PCp;
	PointerCopy pc;
	vector<PointerCopy> * PC;

	if(transID == -1)
	{
		for(sp = 0; sp < seg.size(); sp ++)
		{
			for(i = seg[sp].targetStart; i < seg[sp].targetStart + seg[sp].size; i ++)
				contiVec[contiID][i].B ++;
			if(sp != 0)
				contiVec[contiID][seg[sp].targetStart].C ++;
			if(sp != seg.size() - 1)
				contiVec[contiID][seg[sp].targetStart + seg[sp].size - 1].C --;
		}
	}
	else
	{
		for(sp = 0; sp < seg.size(); sp ++)
		{
			for(i = seg[sp].targetStart, j = seg[sp].sourceStart; i < seg[sp].targetStart + seg[sp].size && j < seg[sp].sourceStart + seg[sp].size; i ++, j ++)
			{
				if(fr)
					contiVec[contiID][i].B = transVec[transID][j] == 0 ? 1 : transVec[transID][j];
				else
					contiVec[contiID][i].B = transVec[transID][transVec[transID].size() - 1 - j] == 0 ? 1 : transVec[transID][transVec[transID].size() - 1 - j];
				contiVec[contiID][i].isTrans = 1;
			}
			contiVec[contiID][seg[sp].targetStart].C = - MAX;
			contiVec[contiID][seg[sp].targetStart + seg[sp].size - 1].C = MAX;
		}
	}

	if(transID == -1)
	{
		for(sp = 0; sp < seg.size(); sp ++)
		{	
			for(i = seg[sp].targetStart; i < seg[sp].targetStart + seg[sp].size - 1; i ++)
			{
				tag = 0;
				PC = & contiVec[contiID][i].PC;
				for(PCp = 0; PCp < (*PC).size(); PCp ++)
					if((*PC)[PCp].pointer == i + 1)
					{
						(*PC)[PCp].copy ++;
						tag = 1;
						break;
					}
				if(tag == 0)
				{
					pc.pointer = i + 1;
					pc.copy = 1;
					(*PC).push_back(pc);
				}
			}

			if(sp != seg.size() - 1)
			{
				tag = 0;
				PC = & contiVec[contiID][seg[sp].targetStart + seg[sp].size - 1].PC;
				for(PCp = 0; PCp < (*PC).size(); PCp ++)
					if((*PC)[PCp].pointer == seg[sp + 1].targetStart)
					{
						(*PC)[PCp].copy ++;
						tag = 1;
						break;
					}
				if(tag == 0)
				{
					pc.pointer = seg[sp + 1].targetStart;
					pc.copy = 1;
					(*PC).push_back(pc);
				}
			}
		}
	}
	else
	{
		for(sp = 0; sp < seg.size(); sp ++)
		{
			for(i = seg[sp].targetStart, j = seg[sp].sourceStart; i < seg[sp].targetStart + seg[sp].size - 1 && j < seg[sp].sourceStart + seg[sp].size - 1; i ++, j ++)
			{
				tag = 0;
				PC = & contiVec[contiID][i].PC;
				for(PCp = 0; PCp < (*PC).size(); PCp ++)
				{
					if((*PC)[PCp].pointer == i + 1)
					{
						if(fr)
							(*PC)[PCp].copy = transVec[transID][j] == 0 ? (*PC)[PCp].copy ++ : (*PC)[PCp].copy + transVec[transID][j];
						else
							(*PC)[PCp].copy = transVec[transID][transVec[transID].size() - 1 - j] == 0 ? (*PC)[PCp].copy ++ : (*PC)[PCp].copy + transVec[transID][transVec[transID].size() - 1 - j];
						tag = 1;
						break;
					}
				}
				
				if(tag == 0)
				{
					pc.pointer = i + 1;
					if(fr)
						pc.copy = transVec[transID][j] == 0 ? 1 : transVec[transID][j];
					else
						pc.copy = transVec[transID][transVec[transID].size() - 1 - j] == 0 ? 1 : transVec[transID][transVec[transID].size() - 1 - j];
					(*PC).push_back(pc);
				}
			}

			if(sp != seg.size() - 1)
			{
				tag = 0;
				PC = & contiVec[contiID][seg[sp].targetStart + seg[sp].size - 1].PC;
				for(PCp = 0; PCp < (*PC).size(); PCp ++)
				{
					if((*PC)[PCp].pointer == seg[sp + 1].targetStart)
					{
						if(fr)
							(*PC)[PCp].copy = transVec[transID][seg[sp + 1].sourceStart] == 0 ? (*PC)[PCp].copy ++ : (*PC)[PCp].copy + transVec[transID][seg[sp + 1].sourceStart];
						else
							(*PC)[PCp].copy = transVec[transID][transVec[transID].size() - 1 - seg[sp + 1].sourceStart] == 0 ? (*PC)[PCp].copy ++ : (*PC)[PCp].copy + transVec[transID][transVec[transID].size() - 1 - seg[sp + 1].sourceStart];
						tag = 1;
						break;
					}
				}

				if(tag == 0)
				{
					pc.pointer = seg[sp + 1].targetStart;
					if(fr)
						pc.copy = transVec[transID][seg[sp + 1].sourceStart] == 0 ? 1 : transVec[transID][seg[sp + 1].sourceStart];
					else
						pc.copy = transVec[transID][transVec[transID].size() - 1 - seg[sp + 1].sourceStart] == 0 ? 1 : transVec[transID][transVec[transID].size() - 1 - seg[sp + 1].sourceStart];
					(*PC).push_back(pc);
				}
			}
		}
	}
}

//read alignment files to de novo transfrags + de novo transfrags -> remaining reads + array A
void mapReads2Transfrags(int insertLow, int insertHigh)
{
	ifstream r1, r2, t;
	ifstream ra1, ra2;
	ifstream rai1, rai2;
	ofstream rao1, rao2;
	int i, tp, tbp, vp1, vp2, pp1, pp2, vp, tag, start, end, batchStartID, finish1, finish2;
	string buf;
	vector<Segment> seg;
	Position pos;
	vector<Position> posVec;
	vector<vector<Position> > v1, v2;
	vector<int> v;
	vector<int> transBaseVec;
	int sourceIDBak, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, fr;

	r1.open("tmp/_reads_1.fa");
	r2.open("tmp/_reads_2.fa");
	t.open("tmp/_transfrags.fa");

	ra1.open("tmp/_reads_1_transfrags.psl");
	ra2.open("tmp/_reads_2_transfrags.psl");

        rai1.open("tmp/_reads_1_contigs.psl");
        rai2.open("tmp/_reads_2_contigs.psl");

        rao1.open("tmp/_remaining_reads_1_contigs.psl");
        rao2.open("tmp/_remaining_reads_2_contigs.psl");

//init array A
	if(t.is_open())
	{
		while(t.good())
		{
			getline(t, buf);
			if(buf[0] == 0) break;
			if(buf[0] == '>')
				transVec.push_back(transBaseVec);
			else
			{
				for(i = 0; i < buf.size(); i ++)
				{
					tp = transVec.size() - 1;
					transVec[tp].push_back(1);
				}
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE! (ERROR 2)" << endl;
		exit(-1);
	}

//	for(tp = 0; tp < transVec.size(); tp ++)
//	{
//		for(tbp = 0; tbp < transVec[tp].size(); tbp ++)
//			cout << transVec[tp][tbp];
//		cout << endl;
//	}

//fill in array A
//find mapped reads from pair 1
	batchStartID = -1;
batch:
	finish1 = parseBLAT(ra1, v1, IDENTITY, IDENTITY, batchStartID);
//	for(vp1 = 0; vp1 < v1.size(); vp1 ++)
//		for(pp1 = 0; pp1 < v1[vp1].size(); pp1 ++)
//			cout << v1[vp1][pp1].targetID << ", " << v1[vp1][pp1].targetStart << ", " << v1[vp1][pp1].targetEnd << endl;

//find mapped reads from pair 2
	finish2 = parseBLAT(ra2, v2, IDENTITY, IDENTITY, batchStartID);
//	for(vp2 = 0; vp2 < v2.size(); vp2 ++)
//		for(pp2 = 0; pp2 < v2[vp2].size(); pp2 ++)
//			cout << v2[vp2][pp2].targetID << ", " << v2[vp2][pp2].targetStart << ", " << v2[vp2][pp2].targetEnd << endl;

//find mapped read pairs, and fill in the coverage array using the mapping information
	for(vp1 = 0, vp2 = 0; vp1 < v1.size() && vp2 < v2.size(); vp1 ++, vp2 ++)
	{
		tag = 0;
		for(pp1 = 0; pp1 < v1[vp1].size(); pp1 ++)
			for(pp2 = 0; pp2 < v2[vp2].size(); pp2 ++)
			{
				start = v1[vp1][pp1].targetStart < v2[vp2][pp2].targetStart ? v1[vp1][pp1].targetStart : v2[vp2][pp2].targetStart;
				end = v1[vp1][pp1].targetEnd > v2[vp2][pp2].targetEnd ? v1[vp1][pp1].targetEnd : v2[vp2][pp2].targetEnd;
				if(v1[vp1][pp1].targetID != -1 && v1[vp1][pp1].targetID == v2[vp2][pp2].targetID && 
				(PE && end - start >= insertLow && end - start <= insertHigh ||
				!PE && v1[vp1][pp1].fr == v2[vp2][pp2].fr && v1[vp1][pp1].fr == 1 && v2[vp2][pp2].targetStart - v1[vp1][pp1].targetEnd >= 0 && v2[vp2][pp2].targetStart - v1[vp1][pp1].targetEnd <= 3 && v2[vp2][pp2].targetEnd - v1[vp1][pp1].targetStart >= insertLow - 3 && v2[vp2][pp2].targetEnd - v1[vp1][pp1].targetStart <= insertHigh + 3 ||
				!PE && v1[vp1][pp1].fr == v2[vp2][pp2].fr && v1[vp1][pp1].fr == 0 && v1[vp1][pp1].targetStart - v2[vp2][pp2].targetEnd >= 0 && v1[vp1][pp1].targetStart - v2[vp2][pp2].targetEnd <= 3 && v1[vp1][pp1].targetEnd - v2[vp2][pp2].targetStart >= insertLow - 3 && v1[vp1][pp1].targetEnd - v2[vp2][pp2].targetStart <= insertHigh + 3))
				{
					tag = 1;
					v.push_back(1);
					for(i = v1[vp1][pp1].targetStart; i < v1[vp1][pp1].targetEnd; i ++)
						transVec[v1[vp1][pp1].targetID][i] ++;
					for(i = v2[vp2][pp2].targetStart; i < v2[vp2][pp2].targetEnd; i ++)
						transVec[v2[vp2][pp2].targetID][i] ++;
					goto cont;
				}
			}
cont:
		if(tag == 0)
			v.push_back(0);
	}

	if(!finish1 && !finish2)
	{
		batchStartID = batchStartID + 1000000;
		goto batch;
	}
//	for(vp = 0; vp < v.size(); vp ++)
//		cout << v[vp] << endl;

//generate unmapped read alignment files
        if(rai1.is_open() && rai2.is_open())
        {
                while(rai1.good())
                {
                        getline(rai1, buf);
                        if(buf[0] == 0) break;

                        parse(buf, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, seg, fr);
                        if(sourceID < v.size())
                        {
                                if(v[sourceID] == 0)
                                {
                                        for(i = 0; i < buf.size(); i ++)
                                                rao1 << buf[i];
                                        rao1 << endl;
                                }
                        }
                        else
                        {
                                for(i = 0; i < buf.size(); i ++)
                                        rao1 << buf[i];
                                rao1 << endl;
                        }
                }

                while(rai2.good())
                {
                        getline(rai2, buf);
                        if(buf[0] == 0) break;

                        parse(buf, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, seg, fr);
                        if(sourceID < v.size())
                        {
                                if(v[sourceID] == 0)
                                {
                                        for(i = 0; i < buf.size(); i ++)
                                                rao2 << buf[i];
                                        rao2 << endl;
                                }
                        }
                        else
                        {
                                for(i = 0; i < buf.size(); i ++)
                                        rao2 << buf[i];
                                rao2 << endl;
                        }
                }
        }
        else
        {
                cout << "CANNOT OPEN FILE! (ERROR 3)" << endl;
                exit(-1);
        }
}

//transfrag alignment file to contigs + contigs + array A -> array B + array C
void mapTransfrags2Contigs(string st, int threshSplit)
{
	ifstream c;
	ifstream ta;
	string buf;
	int i, tag, cp, cbp, vp, pp, seqID = -1, batchStartID, finish;
	BC bc;
	vector<vector<Position> > v;
	vector<BC> contiBaseVec;

	c.open("tmp/_contigs.fa");
//init arrays B and C
	ta.open("tmp/_transfrags_contigs.psl");

	if(c.is_open())
	{
		while(c.good())
		{
			getline(c, buf);
			if(buf[0] == 0) break;
			if(buf[0] == '>')
			{
				contiVec.push_back(contiBaseVec);
			}
			else
			{
				for(i = 0; i < buf.size(); i ++)
				{
					cp = contiVec.size() - 1;
					bc.B = bc.C = 0;
					bc.base = buf[i];
					bc.isTrans = 0;
					bc.segVecP = -1;
					bc.isGene = 0;
					bc.isSecondary = 1;
					bc.oriB = 0;
					contiVec[cp].push_back(bc);
				}
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE! (ERROR 5)" << endl;
		exit(-1);
	}

//	for(cp = 0; cp < contiVec.size(); cp ++)
//		cout << contiVec[cp].size() << endl;

//	for(cp = 0; cp < contiVec.size(); cp ++)
//	{
//		i = 0;
//		for(cbp = 0; cbp < contiVec[cp].size(); cbp ++)
//		{
//			cout << contiVec[cp][cbp].base;
//			i ++;
//			if(i % 60 == 0)
//				cout << endl;
//		}
//		cout << endl;
//	}

//fill in arrays B and C
	batchStartID = -1;
batch:
	finish = parseBLAT(ta, v, IDENTITY, GAPPED_IDENTITY, batchStartID);

	for(vp = 0; vp < v.size(); vp ++)
		for(pp = 0; pp < v[vp].size(); pp ++)
		{
			if(v[vp][pp].targetID != -1)
			{
				adjustStartEnd(v[vp][pp].targetID, v[vp][pp].seg);
				updateContiVec(v[vp][pp].targetID, v[vp][pp].seg, vp, v[vp][pp].fr, threshSplit);
			}
		}

	if(!finish)
	{
		batchStartID = batchStartID + 1000000;
		goto batch;
	}

//	for(cp = 0; cp < contiVec.size(); cp ++)
//	{
//		i = 0;
//		for(cbp = 0; cbp < contiVec[cp].size(); cbp ++)
//		{
//			cout << contiVec[cp][cbp].B;
//			i ++;
//			if(i % 60 == 0)
//				cout << endl;
//		}
//		cout << endl;
//	}
}

//remaining read alignment files to contigs + contigs -> array B + array C
void mapReads2Contigs(int insertLow, int insertHigh)
{
	ifstream ra1, ra2;
	int i, vp1, vp2, pp1, pp2, cp, cbp, start, end, batchStartID, finish1, finish2;
	vector<vector<Position> > v1, v2;
	int contiID, tag, PCp;
	PointerCopy pc;
	vector<PointerCopy> * PC;	

	ra1.open("tmp/_remaining_reads_1_contigs.psl");
	ra2.open("tmp/_remaining_reads_2_contigs.psl");

	batchStartID = -1;
batch:
	finish1 = parseBLAT(ra1, v1, IDENTITY, GAPPED_IDENTITY, batchStartID);
	finish2 = parseBLAT(ra2, v2, IDENTITY, GAPPED_IDENTITY, batchStartID);

	for(vp1 = 0, vp2 = 0; vp1 < v1.size() && vp2 < v2.size(); vp1 ++, vp2 ++)
		for(pp1 = 0; pp1 < v1[vp1].size(); pp1 ++)
			for(pp2 = 0; pp2 < v2[vp2].size(); pp2 ++)
			{
				start = v1[vp1][pp1].targetStart < v2[vp2][pp2].targetStart ? v1[vp1][pp1].targetStart : v2[vp2][pp2].targetStart;
				end = v1[vp1][pp1].targetEnd > v2[vp2][pp2].targetEnd ? v1[vp1][pp1].targetEnd : v2[vp2][pp2].targetEnd;
				if(v1[vp1][pp1].targetID == v2[vp2][pp2].targetID && v1[vp1][pp1].targetID != -1 &&
				(PE && end - start >= insertLow && end - start <= insertHigh * 2 ||
                                !PE && v1[vp1][pp1].fr == v2[vp2][pp2].fr && v1[vp1][pp1].fr == 1 && v2[vp2][pp2].targetStart - v1[vp1][pp1].targetEnd >= 0 && v2[vp2][pp2].targetStart - v1[vp1][pp1].targetEnd <= 3 && v2[vp2][pp2].targetEnd - v1[vp1][pp1].targetStart >= insertLow - 3 && v2[vp2][pp2].targetEnd - v1[vp1][pp1].targetStart <= (insertHigh + 3) * 2 ||
                                !PE && v1[vp1][pp1].fr == v2[vp2][pp2].fr && v1[vp1][pp1].fr == 0 && v1[vp1][pp1].targetStart - v2[vp2][pp2].targetEnd >= 0 && v1[vp1][pp1].targetStart - v2[vp2][pp2].targetEnd <= 3 && v1[vp1][pp1].targetEnd - v2[vp2][pp2].targetStart >= insertLow - 3 && v1[vp1][pp1].targetEnd - v2[vp2][pp2].targetStart <= (insertHigh + 3) * 2))
				{
					adjustStartEnd(v1[vp1][pp1].targetID, v1[vp1][pp1].seg);
					adjustStartEnd(v2[vp2][pp2].targetID, v2[vp2][pp2].seg);
					updateContiVec(v1[vp1][pp1].targetID, v1[vp1][pp1].seg, -1, -1, -1);
					updateContiVec(v2[vp2][pp2].targetID, v2[vp2][pp2].seg, -1, -1, -1);
				}
			}

	if(!finish1 && !finish2)
	{
		batchStartID = batchStartID + 1000000;
		goto batch;
	}
//	for(cp = 0; cp < contiVec.size(); cp ++)
//	{
//		i = 0;
//		for(cbp = 0; cbp < contiVec[cp].size(); cbp ++)
//		{
//			cout << contiVec[cp][cbp].PC[0].pointer << "  ";
//			i ++;
//			if(i % 60 == 0)
//				cout << endl;
//		}
//		cout << endl;
//	}
}

void closeCoverageGaps(int insertLow, int insertHigh)
{
	ifstream ra1, ra2;
        int i, vp1, vp2, pp1, pp2, cp, cbp, start, end, s, e, PCp, tag, tag1, bp, contiID, batchStartID, finish1, finish2;
        vector<vector<Position> > v1, v2;
	PointerCopy pc;
	vector<PointerCopy> * PC;

	for(cp = 0; cp < contiVec.size(); cp ++)
		for(bp = 0; bp < contiVec[cp].size(); bp ++)
			if(contiVec[cp][bp].B > 0)
				contiVec[cp][bp].oriB = 1;

        ra1.open("tmp/_remaining_reads_1_contigs.psl");
       	ra2.open("tmp/_remaining_reads_2_contigs.psl");

	batchStartID = -1;
batch:
        finish1 = parseBLAT(ra1, v1, IDENTITY, GAPPED_IDENTITY, batchStartID);
        finish2 = parseBLAT(ra2, v2, IDENTITY, GAPPED_IDENTITY, batchStartID);

        for(vp1 = 0, vp2 = 0; vp1 < v1.size() && vp2 < v2.size(); vp1 ++, vp2 ++)
                for(pp1 = 0; pp1 < v1[vp1].size(); pp1 ++)
                        for(pp2 = 0; pp2 < v2[vp2].size(); pp2 ++)
			{
				start = v1[vp1][pp1].targetStart < v2[vp2][pp2].targetStart ? v1[vp1][pp1].targetStart : v2[vp2][pp2].targetStart;
                                end = v1[vp1][pp1].targetEnd > v2[vp2][pp2].targetEnd ? v1[vp1][pp1].targetEnd : v2[vp2][pp2].targetEnd;
				s = v1[vp1][pp1].targetEnd < v2[vp2][pp2].targetEnd ? v1[vp1][pp1].targetEnd + 1 : v2[vp2][pp2].targetEnd + 1;
				e = v1[vp1][pp1].targetStart > v2[vp2][pp2].targetStart ? v1[vp1][pp1].targetStart - 1: v2[vp2][pp2].targetStart - 1;
				if(v1[vp1][pp1].targetID == v2[vp2][pp2].targetID && v1[vp1][pp1].targetID != -1 && end - start >= insertLow && end - start <= insertHigh)
				{
					contiID =v1[vp1][pp1].targetID;
					tag = 0;
					for(i = s; i <= e + 1; i ++)
					{
						if(contiVec[contiID][i].oriB == 0 && tag == 0)
						{
							contiVec[contiID][i].B ++;
							contiVec[contiID][i].C ++;
							
							tag1 = 0;
							PC = & contiVec[contiID][i - 1].PC;
							for(PCp = 0; PCp < (*PC).size(); PCp ++)
								if((*PC)[PCp].pointer == i)
								{
									(*PC)[PCp].copy ++;
									tag1 = 1;
									break;
								}
							if(tag1 == 0)
							{
								pc.pointer = i;
								pc.copy = 1;
								contiVec[contiID][i - 1].PC.push_back(pc);
							}

							tag = 1;
						}
						else if(contiVec[contiID][i].oriB  == 0 && tag == 1)
						{
							contiVec[contiID][i].B ++;
						}
						else if(contiVec[contiID][i].oriB > 0 && tag == 1)
						{
							contiVec[contiID][i - 1].C ++;

							tag1 = 0;
							PC = & contiVec[contiID][i - 1].PC;
							for(PCp = 0; PCp < (*PC).size(); PCp ++)
								if((*PC)[PCp].pointer == i)
								{
									(*PC)[PCp].copy ++;
									tag1 = 1;
									break;
								}
							if(tag1 == 0)
							{
								pc.pointer = i;
								pc.copy = 1;
								contiVec[contiID][i - 1].PC.push_back(pc);
							}

							tag = 0;
						}
					}

				}
			}
	if(!finish1 && !finish2)
	{
		batchStartID = batchStartID + 1000000;
		goto batch;
	}
}

void filterSegs(int threshSize, int threshCov)
{
	int sp, ssp;

	for(sp = 0; sp < segVec.size(); sp ++)
	{
		for(ssp = 0; ssp < segVec[sp].size();)
		{
			if((segVec[sp][ssp].size < threshSize || segVec[sp][ssp].cov < threshCov) && segVec[sp][ssp].isTrans == 0)
			{
				segVec[sp].erase(segVec[sp].begin() + ssp);
				continue;
			}
			ssp ++;
		}
	}
}

void filterSegJuncs(int contiID, int threshConn)
{
	int sp, jp;

	for(sp = 0; sp < segVec[contiID].size(); sp ++)
		for(jp = 0; jp < segVec[contiID].size(); jp ++)
			if(segJuncVec[sp][jp].copy < threshConn && segJuncVec[sp][jp].isTrans == 0)
			{
				segJuncVec[sp][jp].copy = 0;
				segJuncVec[sp][jp].isTrans = 0;
				segJuncVec[sp][jp].weight = -1;
				segJuncVec[sp][jp].shared = 0;
			}
}

int isOverlapped(int x1, int y1, int x2, int y2)
{
	if(x1 <= x2 && y1 >= x2 || x1 <= y2 && y1 >= y2 || x1 >= x2 && y1 <= y2)
		return 1;
	else
		return 0;
}

int isContained(int x1, int y1, int x2, int y2)
{
	if(x1 <= x2 && y1 >= y2) 
		return 1;
	else
		return 0;
}

void countCmpOvlReads(int insertLow, int insertHigh)
{
	ifstream ra1, ra2;
	int i, vp, vp1, vp2, pp1, pp2, sp, spBak, spInit, vsp, contiID, start, end, batchStartID, finish1, finish2;
	vector<vector<Position> > v1, v2;

	ra1.open("tmp/_reads_1_contigs.psl");
	ra2.open("tmp/_reads_2_contigs.psl");

	batchStartID = -1;
batch:
	finish1 = parseBLAT(ra1, v1, IDENTITY, GAPPED_IDENTITY, batchStartID);
	finish2 = parseBLAT(ra2, v2, IDENTITY, GAPPED_IDENTITY, batchStartID);

	for(vp1 = 0, vp2 = 0; vp1 < v1.size() && vp2 < v2.size(); vp1 ++, vp2 ++)
	{
		for(pp1 = 0; pp1 < v1[vp1].size(); pp1 ++)
			for(pp2 = 0; pp2 < v2[vp2].size(); pp2 ++)
			{
				start = v1[vp1][pp1].targetStart < v2[vp2][pp2].targetStart ? v1[vp1][pp1].targetStart : v2[vp2][pp2].targetStart;
                                end = v1[vp1][pp1].targetEnd > v2[vp2][pp2].targetEnd ? v1[vp1][pp1].targetEnd : v2[vp2][pp2].targetEnd;

				if(v1[vp1][pp1].targetID == v2[vp2][pp2].targetID && v1[vp1][pp1].targetID != -1 && 
				(PE && end - start >= insertLow && end - start <= insertHigh * 2 ||
                                !PE && v1[vp1][pp1].fr == v2[vp2][pp2].fr && v1[vp1][pp1].fr == 1 && v2[vp2][pp2].targetStart - v1[vp1][pp1].targetEnd >= 0 && v2[vp2][pp2].targetStart - v1[vp1][pp1].targetEnd <= 3 && v2[vp2][pp2].targetEnd - v1[vp1][pp1].targetStart >= insertLow - 3 && v2[vp2][pp2].targetEnd - v1[vp1][pp1].targetStart <= (insertHigh + 3) * 2 ||
                                !PE && v1[vp1][pp1].fr == v2[vp2][pp2].fr && v1[vp1][pp1].fr == 0 && v1[vp1][pp1].targetStart - v2[vp2][pp2].targetEnd >= 0 && v1[vp1][pp1].targetStart - v2[vp2][pp2].targetEnd <= 3 && v1[vp1][pp1].targetEnd - v2[vp2][pp2].targetStart >= insertLow - 3 && v1[vp1][pp1].targetEnd - v2[vp2][pp2].targetStart <= (insertHigh + 3) * 2))
				{
					adjustStartEnd(v1[vp1][pp1].targetID, v1[vp1][pp1].seg);
					adjustStartEnd(v2[vp2][pp2].targetID, v2[vp2][pp2].seg);
					contiID = v1[vp1][pp1].targetID;
					spInit = spBak = -1;
					for(i = v1[vp1][pp1].targetStart; i < v1[vp1][pp1].targetEnd; i ++)
					{
						sp = contiVec[contiID][i].segVecP;
						if(sp != spBak && sp != -1)
						{				
							segVec[contiID][sp].ovl ++;
							for(vsp = 0; vsp < v1[vp1][pp1].seg.size(); vsp ++)
							{
								if(isOverlapped(segVec[contiID][sp].targetStart, segVec[contiID][sp].targetStart + segVec[contiID][sp].size - 1, v1[vp1][pp1].seg[vsp].targetStart, v1[vp1][pp1].seg[vsp].targetStart + v1[vp1][pp1].seg[vsp].size - 1))
								{	
									segVec[contiID][sp].cpt ++;
									break;
								}
							}
							if(spInit == -1) spInit = sp;
							spBak = sp;
						}

					}
					if(PE || !PE && v1[vp1][pp1].fr == 0) spBak = -1;// keep spBak when v1[vp1][pp1].fr == 1
					for(i = v2[vp2][pp2].targetStart; i < v2[vp2][pp2].targetEnd; i ++)
					{
						sp = contiVec[contiID][i].segVecP;
						if((PE || !PE && v1[vp1][pp1].fr == 1) && sp != spBak && sp != -1 ||
						(!PE && v1[vp1][pp1].fr == 0) && sp != spBak && sp != -1 && sp != spInit)
						{
							segVec[contiID][sp].ovl ++;
							for(vsp = 0; vsp < v2[vp2][pp2].seg.size(); vsp ++)
							{
								if(isOverlapped(segVec[contiID][sp].targetStart, segVec[contiID][sp].targetStart + segVec[contiID][sp].size - 1, v2[vp2][pp2].seg[vsp].targetStart, v2[vp2][pp2].seg[vsp].targetStart + v2[vp2][pp2].seg[vsp].size - 1))
								{
									segVec[contiID][sp].cpt ++;
									break;
								}
							}
							spBak = sp;
						}
					}
				}
			}
	}
	if(!finish1 && !finish2)
	{
		batchStartID = batchStartID + 1000000;
		goto batch;
	}
}

double computeWeight(double x, double y)
{
	if(x > 1 || y > 1)
	{
		cout << "EDGE WEIGHT ERROR!" << endl;
		exit(-1);
	}
	if(x > y) return - log(1 - (x - y));
	else return - log(1 - (y - x));
}

void weightSegs(int contiID)
{
	int sp, jp;
	double x, y;

	for(sp = 0; sp < segJuncVec.size(); sp ++)
		for(jp = 0; jp < segJuncVec.size(); jp ++)
			if(segJuncVec[sp][jp].copy != 0)
			{
				if(segVec[contiID][sp].ovl == 0)
					x = 0;
				else
					x = (double)segVec[contiID][sp].cpt / (segVec[contiID][sp].ovl * segVec[contiID][sp].size);
				if(segVec[contiID][jp].ovl == 0)
					y = 0;
				else
					y = (double)segVec[contiID][jp].cpt / (segVec[contiID][jp].ovl * segVec[contiID][jp].size);
				segJuncVec[sp][jp].weight = computeWeight(x, y);				
			}
}

void initPath()
{
	ifstream ta;
	vector<vector<Position> > tVec;
	vector<vector<int> > pathVec;
	vector<int> pVec;
	int cp, tp, pp, sp, ssp, ap, ip, ppp, tag, contiID, batchStartID, finish;

	ta.open("tmp/_transfrags_contigs.psl");

	for(cp = 0; cp < contiVec.size(); cp ++)
		allPathVec.push_back(pathVec);

	batchStartID = -1;
batch:
	finish = parseBLAT(ta, tVec, IDENTITY, GAPPED_IDENTITY, batchStartID);

	for(tp = 0; tp < tVec.size(); tp ++)
		for(pp = 0; pp < tVec[tp].size(); pp ++)
		{
			contiID = tVec[tp][pp].targetID;
			if(contiID != -1)
			{
				adjustStartEnd(tVec[tp][pp].targetID, tVec[tp][pp].seg);
				allPathVec[contiID].push_back(pVec);
				for(sp = 0; sp < tVec[tp][pp].seg.size(); sp ++)
				{
					for(ssp = 0; ssp < segVec[contiID].size(); ssp ++)
						if(isContained(tVec[tp][pp].seg[sp].targetStart, tVec[tp][pp].seg[sp].targetStart + tVec[tp][pp].seg[sp].size - 1, segVec[contiID][ssp].targetStart, segVec[contiID][ssp].targetStart + segVec[contiID][ssp].size - 1))
						{
							tag = 0;
							for(ppp = 0; ppp < allPathVec[contiID][(int)allPathVec[contiID].size() - 1].size(); ppp ++)
								if(allPathVec[contiID][(int)allPathVec[contiID].size() - 1][ppp] == ssp)
									tag = 1;
							if(tag == 0)
								allPathVec[contiID][(int)allPathVec[contiID].size() - 1].push_back(ssp);
						}
				}
				if(allPathVec[contiID][(int)allPathVec[contiID].size() - 1].size() == 0 || allPathVec[contiID][(int)allPathVec[contiID].size() - 1].size() == 1)
					allPathVec[contiID].erase(allPathVec[contiID].begin() + (int)allPathVec[contiID].size() - 1);
			}
		}

	if(!finish)
	{
		batchStartID = batchStartID + 1000000;
		goto batch;
	}

	for(ap = 0; ap < allPathVec.size(); ap ++)
	{
		for(pp = 1; pp < allPathVec[ap].size();)
		{
			tag = 0;
			for(ppp = 0; ppp < pp; ppp ++)
			{
				if(allPathVec[ap][pp] == allPathVec[ap][ppp])
				{
					tag = 1;
					allPathVec[ap].erase(allPathVec[ap].begin() + pp);
					break;
				}
			}
			if(tag == 0)
				pp ++;
		}
	}

	for(ap = 0; ap < allPathVec.size(); ap ++)
		for(pp = 0; pp < allPathVec[ap].size(); pp ++)
			for(ip = 0; ip < allPathVec[ap][pp].size(); ip ++)
				segVec[ap][allPathVec[ap][pp][ip]].shared ++;

//      for(ap = 0; ap < allPathVec.size(); ap ++)
//      {
//              cout << "all path id: " << ap << endl;
//              for(pp = 0; pp < allPathVec[ap].size(); pp ++)
//              {
//                      for(ip = 0; ip < allPathVec[ap][pp].size(); ip ++)
//                              cout << allPathVec[ap][pp][ip];
//                      cout << endl;
//              }
//              cout << endl;
//      }
}

void initSegVec(int threshSize, int threshCov, int threshSplit, int insertLow, int insertHigh)
{
	int tag = 0, cp, cbp, sp, ssp, size, isTrans, cov = 0;
	Segment s;
	vector<Segment> sVec;

	for(cp = 0; cp < contiVec.size(); cp ++)
	{
		tag = 0;
		segVec.push_back(sVec);
		for(cbp = 0; cbp < contiVec[cp].size(); cbp ++)
		{
			if(contiVec[cp][cbp].B > 0 && tag == 0)// seg start
			{
				s.targetStart = cbp;
				s.size = -1;
				s.cov = s.cpt = s.ovl = s.shared = 0;
				s.isTrans = -1;
				segVec[segVec.size() - 1].push_back(s);
				tag = 1;
				if(contiVec[cp][cbp].isTrans == 1)
					isTrans = 1;
				else
					isTrans = 0;
				cov = contiVec[cp][cbp].B;
				contiVec[cp][cbp].segVecP = segVec[segVec.size() - 1].size() - 1;
				size = 1;
			}
			else if(contiVec[cp][cbp].B == 0 && tag == 1 || cbp == contiVec[cp].size() - 1 && tag == 1)// seg end
			{
				segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].size = cbp - segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].targetStart;
				tag = 0;
				if(isTrans > segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].size / 2)
					segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].isTrans = 1;
				else
					segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].isTrans = 0;
				segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].cov = cov / segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].size;
				size = 0;
			}
			else if(contiVec[cp][cbp].B > 0 && tag == 1)
			{
				if(cbp < contiVec[cp].size() - 1 && contiVec[cp][cbp + 1].B > 0 && (contiVec[cp][cbp].C > threshSplit || contiVec[cp][cbp + 1].C < - threshSplit))//end and start
				{
					if(contiVec[cp][cbp].isTrans == 1) isTrans ++;
					segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].size = cbp - segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].targetStart + 1;
					tag = 0;
					if(isTrans > segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].size / 2)
						segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].isTrans = 1;
					else
						segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].isTrans = 0;
					cov = cov + contiVec[cp][cbp].B;
					segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].cov = cov / segVec[segVec.size() - 1][segVec[segVec.size() - 1].size() - 1].size;
					contiVec[cp][cbp].segVecP = segVec[segVec.size() - 1].size() - 1;
					s.targetStart = cbp + 1;
					s.size = -1;
					s.cov = s.cpt = s.ovl = s.shared = 0;
					segVec[segVec.size() - 1].push_back(s);
					tag = 1;
					if(contiVec[cp][cbp + 1].isTrans == 1)
						isTrans = 1;
					else
						isTrans = 0;
					cov = contiVec[cp][cbp + 1].B;
					contiVec[cp][cbp + 1].segVecP = segVec[segVec.size() - 1].size() - 1;
					size = 1;
					cbp ++;
				}

				else// in seg
				{
					if(contiVec[cp][cbp].isTrans == 1)
						isTrans ++;
					cov = cov + contiVec[cp][cbp].B;
					contiVec[cp][cbp].segVecP = segVec[segVec.size() - 1].size() - 1;
					size ++;
				}
			}
		}
	}
//	for(sp = 0; sp < segVec.size(); sp ++)
//		for(ssp = 0; ssp < segVec[sp].size(); ssp ++)
//			cout << segVec[sp][ssp].targetStart << ", " << segVec[sp][ssp].size << ", " << segVec[sp][ssp].cov << ", " << segVec[sp][ssp].cpt << ", " << segVec[sp][ssp].ovl << ", " << segVec[sp][ssp].shared << endl;
//	cout << endl;

	countCmpOvlReads(insertLow, insertHigh);

//	for(sp = 0; sp < segVec.size(); sp ++)
//		for(ssp = 0; ssp < segVec[sp].size(); ssp ++)
//			cout << segVec[sp][ssp].targetStart << ", " << segVec[sp][ssp].size << ", " << segVec[sp][ssp].cov << ", " << segVec[sp][ssp].cpt << ", " << segVec[sp][ssp].ovl << ", " << segVec[sp][ssp].shared << endl;
//	cout << endl;

	filterSegs(threshSize, threshCov);

//	for(sp = 0; sp < segVec.size(); sp ++)
//		for(ssp = 0; ssp < segVec[sp].size(); ssp ++)
//			cout << segVec[sp][ssp].targetStart << ", " << segVec[sp][ssp].size << ", " << segVec[sp][ssp].cov << ", " << segVec[sp][ssp].cpt << ", " << segVec[sp][ssp].ovl << ", " << segVec[sp][ssp].shared << endl;
//	cout << endl;

	initPath();

//	for(sp = 0; sp < segVec.size(); sp ++)
//		for(ssp = 0; ssp < segVec[sp].size(); ssp ++)
//			cout << sp << ", " << segVec[sp][ssp].targetStart << ", " << segVec[sp][ssp].size << ", " << segVec[sp][ssp].cov << ", " << segVec[sp][ssp].cpt << ", " << segVec[sp][ssp].ovl << ", " << segVec[sp][ssp].shared << endl;
//	cout << endl;

//	for(cp = 0; cp < contiVec.size(); cp ++)
//	{
//		for(cbp = 0; cbp < contiVec[cp].size(); cbp ++)
//		{
//			cout << contiVec[cp][cbp].B;
//			if((cbp + 1 )% 60 == 0)
//				cout << endl;
//		}
//		cout << endl;
//	}
//	cout << endl;
}

void initSegJuncVec(int contiID, int threshConn)
{
	int i, sp, jp, PCp, ssp, pp, ip;
	vector<PointerCopy> * PC;
	vector<SegJunc> sjVec;
	SegJunc sj;

	sj.copy = 0;
	sj.isTrans = 0;
	sj.weight = -1;
	sj.shared = 0;
	sj.postponed = 0;
	for(sp = 0; sp < segVec[contiID].size(); sp ++)
	{
		sjVec.clear();
		for(jp = 0; jp < segVec[contiID].size(); jp ++)
			sjVec.push_back(sj);
		segJuncVec.push_back(sjVec);
	}

//	for(sp = 0; sp < segJuncVec.size(); sp ++)
//	{
//		for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//			cout << segJuncVec[sp][jp].copy;
//		cout << endl;
//	}
//	cout << endl;

	for(sp = 0; sp < (int)segVec[contiID].size() - 1; sp ++)
	{
		for(i = segVec[contiID][sp].targetStart; i < segVec[contiID][sp].targetStart + segVec[contiID][sp].size; i ++)
		{
			PC = & contiVec[contiID][i].PC;
			for(PCp = 0; PCp < (*PC).size(); PCp ++)
			{
				for(jp = sp + 1; jp < segVec[contiID].size(); jp ++)
				{
					if((*PC)[PCp].pointer >= segVec[contiID][jp].targetStart && (*PC)[PCp].pointer < segVec[contiID][jp].targetStart + segVec[contiID][jp].size)
					{
						segJuncVec[sp][jp].copy = (*PC)[PCp].copy == 0 ? segJuncVec[sp][jp].copy + 1 : segJuncVec[sp][jp].copy + (*PC)[PCp].copy;
						for(pp = 0; pp < allPathVec[contiID].size(); pp ++)
						{
							for(ip = 0; ip < (int)allPathVec[contiID][pp].size() - 1; ip ++)
							{
								if(allPathVec[contiID][pp][ip] == sp && allPathVec[contiID][pp][ip + 1] == jp)
								{
									segJuncVec[sp][jp].isTrans = 1;
									segJuncVec[sp][jp].shared = segVec[contiID][sp].shared < segVec[contiID][jp].shared ? segVec[contiID][sp].shared : segVec[contiID][jp].shared;
								}
							}
						}
					}
				}
			}
		}
	}

//	for(sp = 0; sp < segVec[contiID].size(); sp ++)
//	{
//		for(jp = 0; jp < segVec[contiID].size(); jp ++)
//			cout << segJuncVec[sp][jp].copy << "  ";
//		cout << endl;
//	}
//	cout << endl;

	weightSegs(contiID);

//	for(sp = 0; sp < segVec[contiID].size(); sp ++)
//	{
//		for(jp = 0; jp < segVec[contiID].size(); jp ++)
//			cout << segJuncVec[sp][jp].copy << "  ";
//		cout << endl;
//	}
//	cout << endl;

	filterSegJuncs(contiID, threshConn);

//	for(sp = 0; sp < segVec[contiID].size(); sp ++)
//	{
//		for(jp = 0; jp < segVec[contiID].size(); jp ++)
//			cout << segJuncVec[sp][jp].copy << "  ";
//		cout << endl;
//	}
//	cout << endl;
}

void generateSingleTrans(ofstream & it, int contiID, int & seqID)
{
        int start, end, i;

        start = segVec[contiID][0].targetStart;
        end = segVec[contiID][0].targetStart + segVec[contiID][0].size;
        it << ">" << seqID ++ << ": " << contiID << ".0: 0" << endl;
        for(i = start; i < end; i ++)
        {
                it << contiVec[contiID][i].base;
                if((i - start + 1) % 60 == 0 && i != end - 1)
                        it << endl;
        }
        it << endl;
}

void checkError()
{
	int sp, jp;

	for(sp = 0; sp < segJuncVec.size(); sp ++)
		for(jp = 0; jp <= sp; jp ++)
			if(segJuncVec[sp][jp].copy > 0)
			{
				cout << "JUNCTION GRAPH ERROR!" << endl;
				for(sp = 0; sp < segJuncVec.size(); sp ++)
				{
					for(jp = 0; jp < segJuncVec.size(); jp ++)
						cout << segJuncVec[sp][jp].copy << "  ";
					cout << endl;
				}
				exit(-1);
			}
}

//----------------------------------------------------------------------------------------------------------
//FUNCTIONS BELOW ARE TO FIND MMPP
//----------------------------------------------------------------------------------------------------------

void connectSink(int x, int y)
{
	int i, tag;

	tag = 0;
	for(i = 0; i < segJuncVec[y].size() - 2; i ++)
		if(segJuncVec[y][i].copy > 0)
			tag = 1;
	if(tag == 0)
	{
		segJuncVec[y][segJuncVec.size() - 1].copy = 1;
		segJuncVec[y][segJuncVec.size() - 1].weight = 0;
	}
}

void connectSource(int x, int y)
{
	int i, tag;

	tag = 0;
	for(i = 0; i < segJuncVec.size() - 2; i ++)
		if(segJuncVec[i][x].copy > 0)
			tag = 1;
	if(tag == 0)
	{
		segJuncVec[segJuncVec.size() - 2][x].copy = 1;
		segJuncVec[segJuncVec.size() - 2][x].weight = 0;
	}

}

void getExistingPath(vector<int> & pathVec, int contiID, int pp)
{
	int ip;

	pathVec.clear();
	for(ip = 0; ip < allPathVec[contiID][pp].size(); ip ++)
	{
		pathVec.push_back(allPathVec[contiID][pp][ip]);
		segVec[contiID][allPathVec[contiID][pp][ip]].shared --;
		if(ip != allPathVec[contiID][pp].size() - 1)
			segJuncVec[allPathVec[contiID][pp][ip]][allPathVec[contiID][pp][ip + 1]].shared --;
	}
}

void addNewNode(vector<int> pathVec, int contiID)
{
	int sp, jp, pp, spp;
	vector<SegJunc> sjVec;
	SegJunc sj;
	vector<Segment> sVec;
	Segment s;

	s.sourceStart = s.targetStart = s.size = s.isTrans = s.cov = s.cpt = s.ovl = -1;
	s.shared = 0;

	for(pp = 0; pp < (int)pathVec.size() - 1; pp ++)
	{
		if(segVec[contiID][pathVec[pp]].path.begin() != segVec[contiID][pathVec[pp]].path.end())
			for(spp = 0; spp < segVec[contiID][pathVec[pp]].path.size(); spp ++)
				s.path.push_back(segVec[contiID][pathVec[pp]].path[spp]);
		else
			s.path.push_back(pathVec[pp]);
		if(segJuncVec[pathVec[pp]][pathVec[pp + 1]].path.begin() != segJuncVec[pathVec[pp]][pathVec[pp + 1]].path.end())
			for(spp = 0; spp < segJuncVec[pathVec[pp]][pathVec[pp + 1]].path.size(); spp ++)
				s.path.push_back(segJuncVec[pathVec[pp]][pathVec[pp + 1]].path[spp]);
	}
	if(pp > 0)
		if(segVec[contiID][pathVec[pp]].path.begin() != segVec[contiID][pathVec[pp]].path.end())
			for(spp = 0; spp < segVec[contiID][pathVec[pp]].path.size(); spp ++)
				s.path.push_back(segVec[contiID][pathVec[pp]].path[spp]);
		else
			s.path.push_back(pathVec[pp]);

	segVec[contiID].push_back(s);

	sj.copy = sj.isTrans = sj.shared = 0;
	sj.weight = -1;
	for(sp = 0; sp < segJuncVec.size(); sp ++)
		segJuncVec[sp].push_back(sj);
	for(jp = 0; jp < segJuncVec[0].size(); jp ++)
		sjVec.push_back(sj);
	segJuncVec.push_back(sjVec);
}

void getPathForward(int x, int y, vector<int> & pathVec, int sinkp, int contiID)
{
        int i, tag = 0;

	maxMatch[x][y] = 0;
        for(i = 0; i < maxMatch[y].size(); i ++)
                if(maxMatch[y][i] == 1 && i != sinkp)
                {
                        maxMatch[y][i] == 0;
			segJuncVec[y][i].shared = 1;
			segVec[contiID][i].shared = 1;
                        pathVec.push_back(i);
                        getPathForward(y, i, pathVec, sinkp, contiID);
                        tag = 1;
                }
        if(tag == 0 && pathVec[pathVec.size() - 1] != y)
                pathVec.push_back(y);
}

void getPathReverse(int x, int y, vector<int> & pathVec, int sourcep, int contiID)
{
        int i, tag = 0;

	maxMatch[x][y] = 0;
        for(i = 0; i < maxMatch[x].size(); i ++)
                if(maxMatch[i][x] == 1 && i != sourcep)
                {
                        maxMatch[i][x] = 0;
			segJuncVec[i][x].shared = 1;
			segVec[contiID][i].shared = 1;
                        pathVec.insert(pathVec.begin(), i);
                        getPathReverse(i, x, pathVec, sourcep, contiID);
                        tag = 1;
                }
        if(tag == 0 && pathVec[0] != x)
                pathVec.insert(pathVec.begin(), x);
}

void addSourceSink(int contiID, int & sourcep, int & sinkp)
{
	vector<int> pathVec;
	int sp, jp;

	addNewNode(pathVec, contiID);
	sourcep = segJuncVec.size() - 1;
        addNewNode(pathVec, contiID);
        sinkp = segJuncVec.size() - 1;
        for(sp = 0; sp < segJuncVec.size() - 2; sp ++)
                for(jp = 0; jp < segJuncVec[sp].size() - 2; jp ++)
                        if(segJuncVec[sp][jp].copy > 0)
                        {
                                connectSource(sp, jp);
                                connectSink(sp, jp);
                        }
//      for(sp = 0; sp < segJuncVec.size(); sp ++)
//              for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//                      if(segJuncVec[sp][jp].copy < 0)
//                              segJuncVec[sp][jp].copy = - segJuncVec[sp][jp].copy;

}

void initIO(vector<int> pathVec, vector<int> & I, vector<int> & O)
{
	int pp, tag, sp;

	I.clear();
	for(pp = 0; pp < pathVec.size(); pp ++)
	{
		tag = 0;
		for(sp = 0; sp < segJuncVec.size(); sp ++)
			if(segJuncVec[sp][pathVec[pp]].copy > 0)
				if(pp > 0)
				{
					if(sp != pathVec[pp - 1])
					{
						I.push_back(pathVec[pp]);
						tag = 1;
						break;
					}
				}
				else// pp == 0
				{
					I.push_back(pathVec[pp]);
					tag = 1;
					break;
				}
		if(tag == 0)
			I.push_back(-1);
	}

	O.clear();
	for(pp = 0; pp < pathVec.size(); pp ++)
	{
		tag = 0;
		for(sp = 0; sp < segJuncVec.size(); sp ++)
			if(segJuncVec[pathVec[pp]][sp].copy > 0)
				if(pp < (int)pathVec.size() - 1)
				{
					if(sp != pathVec[pp + 1])
					{
						O.push_back(pathVec[pp]);
						tag = 1;
						break;
					}
				}
				else// pp == (int)pathVec.size() - 1
				{
					O.push_back(pathVec[pp]);
					tag = 1;
					break;
				}
		if(tag == 0)
			O.push_back(-1);
	}
}

void keepConnection(vector<int> pathVec, vector<int> I, vector<int> O, int contiID)
{
	int ip, op, pp, spp;

	for(ip = 0; ip < (int)I.size() - 1; ip ++)
		for(op = ip + 1; op < O.size(); op ++)
			if(I[ip] != -1 && O[op] != -1)
				if(ip == 0 && op == (int)O.size() - 1)
				{
					segJuncVec[I[ip]][segJuncVec.size() - 1].copy ++;
					segJuncVec[I[ip]][segJuncVec.size() - 1].shared = segJuncVec[I[ip]][segJuncVec.size() - 1].isTrans = segJuncVec[I[ip]][segJuncVec.size() - 1].weight = 0;

					segJuncVec[segJuncVec.size() - 1][O[op]].copy ++;
					segJuncVec[segJuncVec.size() - 1][O[op]].shared = segJuncVec[segJuncVec.size() - 1][O[op]].isTrans;
					for(pp = ip; pp < op; pp ++)
						segJuncVec[segJuncVec.size() - 1][O[op]].weight = segJuncVec[segJuncVec.size() - 1][O[op]].weight + segJuncVec[pathVec[pp]][pathVec[pp + 1]].weight;
				}
				else
				{
					for(pp = ip; pp < op; pp ++)
						if(segJuncVec[pathVec[pp]][pathVec[pp + 1]].shared == 0 && segJuncVec[I[ip]][O[op]].shared == 0)
						{
							segJuncVec[I[ip]][O[op]].copy ++;
							segJuncVec[I[ip]][O[op]].shared = segJuncVec[I[ip]][O[op]].isTrans = segJuncVec[I[ip]][O[op]].weight = 0;
							for(pp = ip; pp < op; pp ++)
								segJuncVec[I[ip]][O[op]].weight = segJuncVec[I[ip]][O[op]].weight + segJuncVec[pathVec[pp]][pathVec[pp + 1]].weight;
							segJuncVec[I[ip]][O[op]].path.clear();
							for(pp = ip; pp < op; pp ++)// could reuse pointers since will break
							{
								if(segJuncVec[pathVec[pp]][pathVec[pp + 1]].path.begin() != segJuncVec[pathVec[pp]][pathVec[pp + 1]].path.end() && !(I[ip] == pathVec[pp] && O[op] == pathVec[pp + 1]))// avoid dead loop assigning a path to itself
									for(spp = 0; spp < segJuncVec[pathVec[pp]][pathVec[pp + 1]].path.size(); spp ++)
									{
										segJuncVec[I[ip]][O[op]].path.push_back(segJuncVec[pathVec[pp]][pathVec[pp + 1]].path[spp]);
									}
								if(pp + 1 != op)
									if(segVec[contiID][pathVec[pp + 1]].path.begin() != segVec[contiID][pathVec[pp + 1]].path.end())
										for(spp = 0; spp < segVec[contiID][pathVec[pp + 1]].path.size(); spp ++)
											segJuncVec[I[ip]][O[op]].path.push_back(segVec[contiID][pathVec[pp + 1]].path[spp]);
									else
										segJuncVec[I[ip]][O[op]].path.push_back(pathVec[pp + 1]);
							}
							break;
						}
				}
}

void deleteExistingPath(int contiID, vector<int> pathVec, vector<int> & I, vector<int> & O)
{
	int pp;

	for(pp = 0; pp < pathVec.size(); pp ++)
	{
		if(segVec[contiID][pathVec[pp]].shared > 0)
		{
			if(pp == 0)
				segJuncVec[pathVec[pp]][segJuncVec.size() - 1].postponed = 1;
			if(pp == pathVec.size() - 1)
				segJuncVec[segJuncVec.size() - 1][pathVec[pp]].postponed = 1;
			I[pp] = O[pp] = -1;
		}
		if(pp != pathVec.size() - 1)
			if(segJuncVec[pathVec[pp]][pathVec[pp + 1]].shared == 0)
			{
				segJuncVec[pathVec[pp]][pathVec[pp + 1]].copy --; 
				segJuncVec[pathVec[pp]][pathVec[pp + 1]].isTrans = 0;
				if(segJuncVec[pathVec[pp]][pathVec[pp + 1]].copy == 0)
					segJuncVec[pathVec[pp]][pathVec[pp + 1]].weight = -1;
			}
	}
}

void initNodeIO(vector<int> N, vector<int> & nI, vector<int> & nO, int np)
{
	int sp;

	nI.clear();
	nO.clear();
	for(sp = 0; sp < segJuncVec.size(); sp ++)
	{
		if(segJuncVec[sp][N[np]].copy > 0)
			nI.push_back(sp);
		if(segJuncVec[N[np]][sp].copy > 0)
			nO.push_back(sp);
	}
}

void keepNodeConnection(int contiID, vector<int> N, vector<int> nI, vector<int> nO, int np, int boundp)
{
	int inp, onp, pp, sp;

	for(inp = 0; inp < nI.size(); inp ++)
		for(onp = 0; onp < nO.size(); onp ++)
		{
			segJuncVec[nI[inp]][nO[onp]].copy = 1;
			segJuncVec[nI[inp]][nO[onp]].shared = segJuncVec[nI[inp]][nO[onp]].isTrans = 0;
			segJuncVec[nI[inp]][nO[onp]].weight = segJuncVec[nI[inp]][N[np]].weight + segJuncVec[N[np]][nO[onp]].weight;
			segJuncVec[nI[inp]][nO[onp]].path.clear();

			if(nO[onp] == segJuncVec.size() - 1 || segJuncVec[N[np]][nO[onp]].postponed)
				for(pp = 0; pp < segJuncVec[nI[inp]][N[np]].path.size(); pp ++)
					segJuncVec[nI[inp]][nO[onp]].path.push_back(segJuncVec[nI[inp]][N[np]].path[pp]);
			else if(nI[inp] == segJuncVec.size() - 1 || segJuncVec[nI[inp]][N[np]].postponed)
				for(pp = 0; pp < segJuncVec[N[np]][nO[onp]].path.size(); pp ++)
					segJuncVec[nI[inp]][nO[onp]].path.push_back(segJuncVec[N[np]][nO[onp]].path[pp]);
			else
			{
				for(pp = 0; pp < segJuncVec[nI[inp]][N[np]].path.size(); pp ++)
					segJuncVec[nI[inp]][nO[onp]].path.push_back(segJuncVec[nI[inp]][N[np]].path[pp]);
				if(segVec[contiID][N[np]].path.begin() == segVec[contiID][N[np]].path.end())
					segJuncVec[nI[inp]][nO[onp]].path.push_back(N[np]);
				else
					for(pp = 0; pp < segVec[contiID][N[np]].path.size(); pp ++)
						segJuncVec[nI[inp]][nO[onp]].path.push_back(segVec[contiID][N[np]].path[pp]);
				for(pp = 0; pp < segJuncVec[N[np]][nO[onp]].path.size(); pp ++)
					segJuncVec[nI[inp]][nO[onp]].path.push_back(segJuncVec[N[np]][nO[onp]].path[pp]);
			}
		}

	for(sp = 0; sp < segJuncVec.size(); sp ++)
	{
		if(segJuncVec[sp][N[np]].copy > 0)
		{
			segJuncVec[sp][N[np]].copy = segJuncVec[sp][N[np]].shared = segJuncVec[sp][N[np]].isTrans = segJuncVec[sp][N[np]].postponed = 0;
			segJuncVec[sp][N[np]].weight = -1;
		}
		if(segJuncVec[N[np]][sp].copy > 0)
		{
			segJuncVec[N[np]][sp].copy = segJuncVec[N[np]][sp].shared = segJuncVec[N[np]][sp].isTrans = segJuncVec[N[np]][sp].postponed = 0;
			segJuncVec[N[np]][sp].weight = -1;
		}
	}
}

void findMaxMatch()
{
	int sp, jp;

	ListGraph bp;
	ListGraph::EdgeMap<double> wt(bp);
	ListGraph::Node nd;
	ListGraph::Edge eg;
	vector<ListGraph::Node> u, v;
	vector<vector<ListGraph::Edge> > e;
	vector<ListGraph::Edge> egVec;
	vector<int> intVec;

	maxMatch.clear();

	for(sp = 0; sp < segJuncVec.size(); sp ++)
	{
		u.push_back(nd);
		v.push_back(nd);
		u[sp] = bp.addNode();
		v[sp] = bp.addNode();
	}

	for(sp = 0; sp < segJuncVec.size(); sp ++)
	{
		maxMatch.push_back(intVec);
		e.push_back(egVec);
		for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
		{
			maxMatch[sp].push_back(0);
			e[sp].push_back(eg);
			if(segJuncVec[sp][jp].copy > 0)
			{
				e[sp][jp] = bp.addEdge(u[sp], v[jp]);
				wt[e[sp][jp]] = - segJuncVec[sp][jp].weight;
			}
			else
			{
				e[sp][jp] = bp.addEdge(u[sp], v[jp]);
				wt[e[sp][jp]] = - MAX;
			}
		}
	}

	MaxWeightedPerfectMatching<ListGraph, ListGraph::EdgeMap<double> > mwpm(bp, wt);
	mwpm.run();

	ListGraph::EdgeIt ei(bp);
	for(sp = segJuncVec.size() - 1; sp >= 0; sp --)
		for(jp = segJuncVec[sp].size() - 1; jp >= 0; jp --)
		{	
			maxMatch[sp][jp] = segJuncVec[sp][jp].copy == 0 ? 0 : mwpm.matching(ei);
			++ei;
		}
}

void updateExistingPath(int contiID, int sourcep, int sinkp)
{
	int sp, jp;
	vector<int> pathVec;

	allPathVec[contiID].clear();
	for(sp = 0; sp < maxMatch.size(); sp ++)
		for(jp = 0; jp < maxMatch[sp].size(); jp ++)
			if(maxMatch[sp][jp] == 1 && sp != sourcep && jp != sinkp)
			{
				pathVec.clear();
				pathVec.push_back(sp);
				pathVec.push_back(jp);
				getPathForward(sp, jp, pathVec, sinkp, contiID);
				getPathReverse(sp, jp, pathVec, sourcep, contiID);
				allPathVec[contiID].push_back(pathVec);
				segJuncVec[sp][jp].shared = 1;
				segVec[contiID][sp].shared = segVec[contiID][jp].shared = 1;
			}
}

int finish(int contiID)
{
	int pp;

	for(pp = 0; pp < allPathVec[contiID].size(); pp ++)
		if(allPathVec[contiID][pp].size() > 1)
			return 0;
	return 1;
}

void generateTranscripts(ofstream & it, int contiID, int sourcep, int sinkp, int & seqID)
{
	int i, transID = 0, ssp, sp, bp, counter, spBak, tp, isTrans;
	vector<vector<int> > transVec;
	vector<int> trans;

	for(i = 0; i < segJuncVec[sourcep].size(); i ++)
	{
		if(segJuncVec[sourcep][i].copy > 0 && i != sinkp)
		{
			spBak = -1;
			trans.clear();
			isTrans = 0;
			for(ssp = 0; ssp < segJuncVec[sourcep][i].path.size(); ssp ++)
			{
				if(segJuncVec[sourcep][i].path[ssp] != spBak)
				{
					trans.push_back(segJuncVec[sourcep][i].path[ssp]);
					spBak = segJuncVec[sourcep][i].path[ssp];
				}
				if(segJuncVec[sourcep][i].isTrans) isTrans = 1;
			}
			if(segVec[contiID][i].path.begin() != segVec[contiID][i].path.end())
			{
				for(ssp = 0; ssp < segVec[contiID][i].path.size(); ssp ++)
				{
					if(segVec[contiID][i].path[ssp] != spBak)
					{
						trans.push_back(segVec[contiID][i].path[ssp]);
						spBak = segVec[contiID][i].path[ssp];
					}
					if(segVec[contiID][i].isTrans) isTrans = 1;
				}
			}
			else
			{
				if(i != spBak)
				{
					trans.push_back(i);
					spBak = i;
				}
				if(segVec[contiID][i].isTrans) isTrans = 1;
			}
			for(ssp = 0; ssp < segJuncVec[i][sinkp].path.size(); ssp ++)
			{
				if(segJuncVec[i][sinkp].path[ssp] != spBak)
				{
					trans.push_back(segJuncVec[i][sinkp].path[ssp]);
					spBak = segJuncVec[i][sinkp].path[ssp];
				}
				if(segVec[contiID][i].isTrans) isTrans = 1;
			}
			for(tp = 0; tp < transVec.size(); tp ++)
				if(transVec[tp] == trans)
					goto conti;
			if(isTrans == 0) goto conti;
			transVec.push_back(trans);
//check duplication
			spBak = -1;
			it << ">" << seqID ++ << ": " << contiID << "." << transID ++ << ": ";
			for(ssp = 0; ssp < segJuncVec[sourcep][i].path.size(); ssp ++)
				if(segJuncVec[sourcep][i].path[ssp] != spBak)
				{
					it << segJuncVec[sourcep][i].path[ssp];
					spBak = segJuncVec[sourcep][i].path[ssp];
				}
			if(segVec[contiID][i].path.begin() != segVec[contiID][i].path.end())
			{
				for(ssp = 0; ssp < segVec[contiID][i].path.size(); ssp ++)
					if(segVec[contiID][i].path[ssp] != spBak)
					{
						it << segVec[contiID][i].path[ssp];
						spBak = segVec[contiID][i].path[ssp];
					}
			}
			else
				if(i != spBak)
				{
					it << i;
					spBak = i;
				}
			for(ssp = 0; ssp < segJuncVec[i][sinkp].path.size(); ssp ++)
				if(segJuncVec[i][sinkp].path[ssp] != spBak)
				{
					it << segJuncVec[i][sinkp].path[ssp];
					spBak = segJuncVec[i][sinkp].path[ssp];
				}
			it << endl;
//output title
			counter = 0;
			spBak = -1;
			for(ssp = 0; ssp < segJuncVec[sourcep][i].path.size(); ssp ++)
			{
                                sp = segJuncVec[sourcep][i].path[ssp];
				if(sp != spBak)
				{
                                	for(bp = segVec[contiID][sp].targetStart; bp < segVec[contiID][sp].targetStart + segVec[contiID][sp].size - 1; bp ++)
                                	{
						if(counter % 60 == 0 && counter != 0)
							it << endl;
						counter ++;
                                	        it << contiVec[contiID][bp].base;
                                	}
					spBak = sp;
				}
			}
			if(segVec[contiID][i].path.begin() != segVec[contiID][i].path.end())
				for(ssp = 0; ssp < segVec[contiID][i].path.size(); ssp ++)
				{
					sp = segVec[contiID][i].path[ssp];
					if(sp != spBak)
					{
						for(bp = segVec[contiID][sp].targetStart; bp < segVec[contiID][sp].targetStart + segVec[contiID][sp].size - 1; bp ++)
						{
							if(counter % 60 == 0 && counter != 0)
								it << endl;
							counter ++;
							it << contiVec[contiID][bp].base;
						}
						spBak = sp;
					}
				}
			else
				if(i != spBak)
				{
					for(bp = segVec[contiID][i].targetStart; bp < segVec[contiID][i].targetStart + segVec[contiID][i].size - 1; bp ++)
					{
						if(counter % 60 == 0 && counter != 0)
							it << endl;
						counter ++;
						it << contiVec[contiID][i].base;
					}
					spBak = i;
				}
                        for(ssp = 0; ssp < segJuncVec[i][sinkp].path.size(); ssp ++)
                        {
                                sp = segJuncVec[i][sinkp].path[ssp];
				if(sp != spBak)
				{
                                	for(bp = segVec[contiID][sp].targetStart; bp < segVec[contiID][sp].targetStart + segVec[contiID][sp].size - 1; bp ++)
                                	{
						if(counter % 60 == 0 && counter != 0)
							it << endl;
						counter ++;
                                	        it << contiVec[contiID][bp].base;
                                	}
					spBak = sp;
				}
                        }
			it << endl;
//output sequence
conti:;
		}		
	}
}

void setCopy()
{
	int sp, jp;

	for(sp = 0; sp < segJuncVec.size(); sp ++)
		for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
			if(segJuncVec[sp][jp].copy > 0)
				segJuncVec[sp][jp].copy = 1;
}

void transformSegJuncVec(ofstream & it, int contiID, int & seqID)
{
	int pp, ip, op, sourcep, sinkp, sp, jp, spp;
	vector<int> pathVec;
	vector<int> I, O, nI, nO;

	if(segJuncVec.size() == 0)
		return;

	setCopy();

	addSourceSink(contiID, sourcep, sinkp);

cont:
//	cout << endl;
//	cout << "******cont******" << endl;
//	cout << endl;

	for(pp = 0; pp < allPathVec[contiID].size(); pp ++)
	{
		getExistingPath(pathVec, contiID, pp);

//		cout << "======itrt======" << endl;
//		cout << endl;

//		cout << "------sjv0------" << endl;
//		for(sp = 0; sp < segJuncVec.size(); sp ++)
//		{
//			for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//				cout << segJuncVec[sp][jp].copy << "  ";
//			cout << endl;
//		}
//		cout << "----------------" << endl;
//		cout << endl;

		addNewNode(pathVec, contiID);

//		cout << "------sjv1------" << endl;
//		for(sp = 0; sp < segJuncVec.size(); sp ++)
//		{
//			for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//				cout << segJuncVec[sp][jp].copy << "  ";
//			cout << endl;
//		}
//		cout << "----------------" << endl;
//		cout << endl;

		initIO(pathVec, I, O);
		keepConnection(pathVec, I, O, contiID);

//		cout << "------sjv2------" << endl;
//		for(sp = 0; sp < segJuncVec.size(); sp ++)
//		{
//			for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//				cout << segJuncVec[sp][jp].copy << "  ";
//			cout << endl;
//		}
//		cout << "----------------" << endl;
//		cout << endl;

		deleteExistingPath(contiID, pathVec, I, O);

//		cout << "------sjv3------" << endl;
//		for(sp = 0; sp < segJuncVec.size(); sp ++)
//		{
//			for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//				cout << segJuncVec[sp][jp].copy << "  ";
//			cout << endl;
//		}
//		cout << "----------------" << endl;
//		cout << endl;

		for(ip = 0; ip < I.size(); ip ++)
		{
			if(I[ip] != -1)
			{
				initNodeIO(I, nI, nO, ip);
				keepNodeConnection(contiID, I, nI, nO, ip, 0);
			}
		}
		for(op = 0; op < O.size(); op ++)
		{
			if(O[op] != -1)
			{
				initNodeIO(O, nI, nO, op);
				keepNodeConnection(contiID, O, nI, nO, op, 1);
			}
		}

//		cout << "------sjv4------" << endl;
//		for(sp = 0; sp < segJuncVec.size(); sp ++)
//		{
//			for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//				cout << segJuncVec[sp][jp].copy << "  ";
//			cout << endl;
//		}
//		cout << "----------------" << endl;
//		cout << endl;

//		cout << "================" << endl;
//		cout << endl;
	}

	findMaxMatch();

//	cout << "-------weit------" << endl;
//	for(sp = 0; sp < segJuncVec.size(); sp ++)
//	{
//		for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//			cout << segJuncVec[sp][jp].weight << "  ";
//		cout << endl;
//	}
//	cout << "----------------" << endl;
//	cout << endl;

//	cout << "------mmc------" << endl;
//	for(sp = 0; sp < maxMatch.size(); sp ++)
//	{
//		for(jp = 0; jp < maxMatch[sp].size(); jp ++)
//			cout << maxMatch[sp][jp] << "  ";
//		cout << endl;
//	}
//	cout << "----------------" << endl;
//	cout << endl;

//	cout << "------cvtd------" << endl;
//	for(sp = 0; sp < segJuncVec.size(); sp ++)
//	{
//		for(jp = 0; jp < segJuncVec[sp].size(); jp ++)
//		{
//			if(/*segJuncVec[sp][jp].copy == 1 &&*/ segJuncVec[sp][jp].path.begin() != segJuncVec[sp][jp].path.end())
//			{
//				cout << sp << "." << jp << ": ";
//				for(spp = 0; spp < segJuncVec[sp][jp].path.size(); spp ++)
//					cout << segJuncVec[sp][jp].path[spp] << "  ";
//				cout << endl;
//			}
//		}
//	}

//	cout << "+" << endl;

//	for(sp = 0; sp < segVec[contiID].size(); sp ++)
//	{
//		if(segVec[contiID][sp].path.begin() != segVec[contiID][sp].path.end())
//		{
//			cout << sp << ": ";
//			for(spp = 0; spp < segVec[contiID][sp].path.size(); spp ++)
//				cout << segVec[contiID][sp].path[spp];
//			cout << endl;
//		}
//	}
//	cout << "----------------" << endl;
//	cout << endl;

	updateExistingPath(contiID, sourcep, sinkp);

//	cout << "****************" << endl;
//	cout << endl;

	if(!finish(contiID)) goto cont;

	generateTranscripts(it, contiID, sourcep, sinkp, seqID);
	segJuncVec.clear();
}

void formalizeInput(ifstream & in, string file)
{
	string buf;
	int i, sp, bp;
	unsigned long seqID = 0;
	ofstream out;
	vector<vector<char> > seqVec;
	vector<char> sv;

	out.open(file.c_str());

	if(file == "tmp/_transfrags.fa" || file == "tmp/_contigs.fa")
	{
		if(in.is_open())
		{
			while(in.good())
			{
				getline(in, buf);
				if(buf[0] == 0)
					break;
				if(buf[0] == '>')
					seqVec.push_back(sv);
				else
					for(i = 0; i < buf.size(); i ++)
						seqVec[seqVec.size() - 1].push_back(buf[i]);
			}

			for(sp = 0; sp < seqVec.size(); sp ++)
			{
				if(seqVec[sp].size() >= 60)
				{
					out << ">" << seqID ++ << endl;
					for(bp = 0; bp < seqVec[sp].size(); bp ++)
					{
						out << seqVec[sp][bp];
						if((bp + 1) % 60 == 0 || bp == seqVec[sp].size() - 1)
							out << endl;
					}
				}
			}
		}
		else
		{
			cout << "CANNOT OPEN FILE! (ERROR 6)" << endl;
			exit(-1);
		}
	}
	else
	{
		if(in.is_open())
		{
			while(in.good())
			{
				getline(in, buf);
				if(buf[0] == 0)
					break;
				if(buf[0] == '>')
				{
					getline(in, buf);
					out << ">" << seqID ++ << endl;
					for(i = 0; i < buf.size(); i ++)
						out << buf[i];
					out << endl;
				}
				else
				{
					for(i = 0; i < buf.size(); i ++)
						out << buf[i];
					out << endl;
				}		
			}
		}
		else
		{
			cout << "CANNOT OPEN FILE! (ERROR 6)" << endl;
			exit(-1);
		}
	}
}

void formalizeSingleReads(ifstream & in, string file1, string file2, int & insertLow, int & insertHigh)
{
	string buf;
	int i, tag = 0;
	unsigned long seqID = 0;
	ofstream out1;
	ofstream out2;

	out1.open(file1.c_str());
	out2.open(file2.c_str());

	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0)
				break;
			if(buf[0] == '>')
			{
				out1 << ">" << seqID << endl;
				out2 << ">" << seqID ++ << endl;
			}
			else
			{
				for(i = 0; i < buf.size() / 2; i ++)
					out1 << buf[i];
				out1 << endl;
				for(i = buf.size() / 2; i < buf.size(); i ++)
					out2 << buf[i];
				out2 << endl;
				if(tag == 0)
				{
					insertLow = insertHigh = buf.size();
					tag = 1;
				}
				else if(tag == 1 && insertLow != buf.size())
				{
					cout << "INVALID INPUT FILE!" << endl;
					exit(-1);
				}
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE! (ERROR 7)" << endl;
		exit(-1);
	}
}

typedef struct SeqMappedStruct
{
	vector<char> seq;
	int mapped;
} SeqMapped;

void generateFinalTranscripts(ofstream & it)
{
	ifstream tf, tc, tfa, tca;
	string buf;
	vector<SeqMapped> transfrags;
	vector<SeqMapped> transcripts;
	SeqMapped sm;

	vector<vector<Position> > v, u;
	int vp, up, pp, seqID = 0, counter, i, tvp, tcp, tfp, pointerBak, tag;

	system("blat tmp/_transcripts.fa tmp/_transfrags.fa -noHead tmp/_transfrags_transcripts.psl >> tmp/blat_doc.txt");
	tf.open("tmp/_transfrags.fa");
	tc.open("tmp/_transcripts.fa");
	tfa.open("tmp/_transfrags_transcripts.psl");

	sm.mapped = 0;
	if(tf.is_open())
	{
		while(tf.good())
		{
			getline(tf, buf);
			if(buf[0] == 0)
				break;
			if(buf[0] == '>')
				transfrags.push_back(sm);
			else
				for(i = 0; i < buf.size(); i ++)
					transfrags[transfrags.size() - 1].seq.push_back(buf[i]);
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE! (ERROR 8)" << endl;
		exit(-1);
	}

        if(tc.is_open())
        {
                while(tc.good())
                {
                        getline(tc, buf);
                        if(buf[0] == 0)
                                break;
                        if(buf[0] == '>')
                                transcripts.push_back(sm);
                        else
                                for(i = 0; i < buf.size(); i ++)
                                        transcripts[transcripts.size() - 1].seq.push_back(buf[i]);
                }
        }
        else
        {
                cout << "CANNOT OPEN FILE! (ERROR 9)" << endl;
                exit(-1);
        }


	parseBLAT(tfa, v, IDENTITY, IDENTITY, 0);
	tfa.clear(); tfa.seekg(0);
	parseBLAT(tfa, u, HIGH_IDENTITY, HIGH_IDENTITY, 0);

	for(vp = 0; vp < v.size(); vp ++)
		for(pp = 0; pp < v[vp].size(); pp ++)
			if(v[vp][pp].targetID != -1)
				transcripts[v[vp][pp].targetID].mapped = 1;

	for(up = 0; up < u.size(); up ++)
		for(pp = 0; pp < u[up].size(); pp ++)
			if(u[up][pp].targetID != -1)
				transfrags[u[up][pp].sourceID].mapped = 1;

	for(tvp = 0; tvp < transcripts.size(); tvp ++)
	{
		pointerBak = -1;
		counter = 0;
		it << ">" << seqID ++ << endl;
		for(tcp = 0; tcp < transcripts[tvp].seq.size(); tcp ++)
		{
			it << transcripts[tvp].seq[tcp];
			counter ++;
			if(counter % 60 == 0 && counter != 0)
				it << endl;
		}
		if(counter % 60 != 0)
			it << endl;
	}

	for(tvp = 0; tvp < transfrags.size(); tvp ++)
		if(transfrags[tvp].mapped == 0)
		{
			counter = 0;
			it << ">" << seqID ++ << endl;
			for(tfp = 0; tfp < transfrags[tvp].seq.size(); tfp ++)
			{
				it << transfrags[tvp].seq[tfp];
				counter ++;
				if(counter % 60 == 0 && counter != 0)
					it << endl;
			}
			if(counter % 60 != 0)
				it << endl;
		}
}

//----------------------------------------------------------------------------------------------------------
//FUNCTIONS TO FIND MMPP END
//----------------------------------------------------------------------------------------------------------

void clearUp()
{
	contiVec.clear();
	segVec.clear();
	allPathVec.clear();
}

void * task0(void * arg)
{
	system("blat tmp/_transfrags.fa tmp/_reads_1.fa -noHead -minScore=90 -fastMap tmp/_reads_1_transfrags.psl > tmp/blat_doc.txt");
}

void * task1(void * arg)
{
	system("blat tmp/_transfrags.fa tmp/_reads_2.fa -noHead -minScore=90 -fastMap tmp/_reads_2_transfrags.psl >> tmp/blat_doc.txt");
}

void * task2(void * arg)
{
	system("blat tmp/_contigs.fa tmp/_transfrags.fa -noHead tmp/_transfrags_contigs.psl >> tmp/blat_doc.txt");
}

void * task3(void * arg)
{
	system("blat tmp/_contigs.fa tmp/_reads_1.fa -noHead tmp/_reads_1_contigs.psl >> tmp/blat_doc.txt");
}

void * task4(void * arg)
{
	system("blat tmp/_contigs.fa tmp/_reads_2.fa -noHead tmp/_reads_2_contigs.psl  >> tmp/blat_doc.txt");
}

void nonParallelMap()
{
	system("blat tmp/_transfrags.fa tmp/_reads_1.fa -noHead -minScore=90 -fastMap tmp/_reads_1_transfrags.psl > tmp/blat_doc.txt");
	system("blat tmp/_transfrags.fa tmp/_reads_2.fa -noHead -minScore=90 -fastMap tmp/_reads_2_transfrags.psl >> tmp/blat_doc.txt");
	system("blat tmp/_contigs.fa tmp/_transfrags.fa -noHead tmp/_transfrags_contigs.psl >> tmp/blat_doc.txt");
	system("blat tmp/_contigs.fa tmp/_reads_1.fa -noHead tmp/_reads_1_contigs.psl >> tmp/blat_doc.txt");
	system("blat tmp/_contigs.fa tmp/_reads_2.fa -noHead tmp/_reads_2_contigs.psl >> tmp/blat_doc.txt");
}

void parallelMap()
{
	pthread_t t0, t1, t2, t3, t4;

        if(pthread_create(&t0, NULL, task0, NULL) != 0) {nonParallelMap(); return;}
	if(pthread_create(&t1, NULL, task1, NULL) != 0) {nonParallelMap(); return;}
	if(pthread_create(&t2, NULL, task2, NULL) != 0) {nonParallelMap(); return;}
	if(pthread_create(&t3, NULL, task3, NULL) != 0) {nonParallelMap(); return;}
	if(pthread_create(&t4, NULL, task4, NULL) != 0) {nonParallelMap(); return;}

	if(pthread_join(t0, NULL) != 0) {perror("Thread join faied"); exit(EXIT_FAILURE);}
	if(pthread_join(t1, NULL) != 0) {perror("Thread join faied"); exit(EXIT_FAILURE);}
	if(pthread_join(t2, NULL) != 0) {perror("Thread join faied"); exit(EXIT_FAILURE);}
	if(pthread_join(t3, NULL) != 0) {perror("Thread join faied"); exit(EXIT_FAILURE);}
	if(pthread_join(t4, NULL) != 0) {perror("Thread join faied"); exit(EXIT_FAILURE);}
}

void print()
{
	cout << "BRANCH --read1 reads_1.fa --read2 reads_2.fa --transfrag transfrags.fa --contig contigs.fa --transcript transcripts.fa [--insertLow insertLow --insertHigh insertHigh --threshSize threshSize --threshCov threshCov --threshSplit threshSplit --threshConn threshConn --closeGap --noAlignment --lowEukaryote]" << endl;
	cout << "Inputs: " << endl;
	cout << "--read1 is the first pair of PE RNA reads or single-end RNA reads in fasta format" << endl;
	cout << "--read2 is the second pair of PE RNA reads in fasta format" << endl;
	cout << "--transfrag is the de novo RNA transfrags to be extended" << endl;
	cout << "--contig is the DNA contigs or the genes of close related species" << endl;
	cout << "Output: " << endl;
	cout << "--transcript is the extended de novo transfrags" << endl;
	cout << "Options: " << endl;
	cout << "--insertLow is the lower bound of insert length (highly recommended; default: 0)" << endl;
	cout << "--insertHigh is the upper bound of insert length (highly recommended; default: 99999)" << endl;
	cout << "--threshSize is the minimum size of a genome region that could be identified as an exon (default: 2 bp)" << endl;
	cout << "--threshCov is the minimum coverage of a genome region that could be identified as an exon (default: 2)" << endl;
	cout << "--threshSplit is the minimum upstream and downstream junction coverages to split a genome region into more than one exons (default: 2)" << endl;
	cout << "--threshConn is the minimum connectivity of two exons that could be identified as a junction (default: 2)" << endl;
	cout << "--closeGap closes sequencing gaps using PE read information (default: none)" << endl;
	cout << "--noAlignment skips the initial time-consuming alignment step, if all the alignment files have been provided in tmp directory (default: none)" << endl;
	cout << "--lowEukaryote runs in a different mode for low eukaryotes with rare splice variants (default: none)" << endl;
}

int main(int argc, char * argv[])
{
	int contiID, i, seqID, tagRead1 = 0, tagRead2 = 0, tagTransfrag = 0, tagContig = 0, threshSize = 2, threshCov = 2, threshConn = 2, threshSplit = 2, tagThreshSize = 0, tagThreshCov = 0, tagThreshConn = 0, tagTranscript = 0, tagThreshSplit = 0, tagInsertLow = 0, tagInsertHigh = 0, insertLow = 0, insertHigh = MAX, tagGene = 0, tagCloseGap = 0, tagNoAlignment = 0, tagLowEukaryote = 0;
	ifstream r1, r2, t, ci;
	ofstream g, it, iit, co;
	string s, st, command, sCheck;
	stringstream ssCheck;

	cout << "BRANCH: boosting RNA-Seq assemblies with partial or related genome sequences" << endl;
	cout << "By Ergude Bao, CS Department, UC-Riverside. All Rights Reserved" << endl << endl;

	for(i = 1; i < argc; i ++)
	{
		s = argv[i];
		if(s == "--read1")
		{
			if(tagRead1 == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			r1.open(argv[++ i]);
			if(!r1.is_open())
			{
				cout << "CANNOT OPEN FILE! (ERROR 10)" << endl;
				print();
				return 0;
			}
			tagRead1 = 1;
		}
		else if(s == "--read2")
		{
			if(tagRead2 == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			r2.open(argv[++ i]);
			if(!r2.is_open())
			{
				cout << "CANNOT OPEN FILE! (ERROR 11)" << endl;
				print();
				return 0;
			}
			tagRead2 = 1;
		}
		else if(s == "--transfrag")
		{
			if(tagTransfrag == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			t.open(argv[++ i]);
			if(!t.is_open())
			{
				cout << "CANNOT OPEN FILE! (ERROR 12)" << endl;
				print();
				return 0;
			}
			st = argv[i];
			tagTransfrag = 1;
		}
		else if(s == "--contig")
		{
			if(tagContig == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			ci.open(argv[++ i]);
			if(!ci.is_open())
			{
				cout << "CANNOT OPEN FILE! (ERROR 13)" << endl;
				print();
				return 0;
			}
			tagContig = 1;
		}
                else if(s == "--transcript")
                {
                        if(tagTranscript == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        it.open(argv[++ i]);
                        if(!it.is_open())
                        {
                                cout << "CANNOT OPEN FILE! (ERROR 14)" << endl;
                                print();
                                return 0;
                        }
                        tagTranscript = 1;
                }
                else if(s == "--gene")
                {
                        if(tagGene == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        g.open(argv[++ i]);
                        if(!g.is_open())
                        {
                                cout << "CANNOT OPEN FILE! (ERROR 15)" << endl;
                                print();
                                return 0;
                        }
                        tagGene = 1;
                }
		else if(s == "--threshSize")
		{
			if(tagThreshSize == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			threshSize = atoi(argv[++ i]);
			ssCheck.str("");
			ssCheck << threshSize;
			sCheck = ssCheck.str();
			if(sCheck != argv[i])
			{
				print();
				return 0;
			}
			tagThreshSize = 1;
		}
		else if(s == "--threshCov")
		{
			if(tagThreshCov == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			threshCov = atoi(argv[++ i]);
			ssCheck.str("");
			ssCheck << threshCov;
			sCheck = ssCheck.str();
			if(sCheck != argv[i])
			{
				print();
				return 0;
			}
			tagThreshCov = 1;
		}
		else if(s == "--threshConn")
		{
			if(tagThreshConn == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			threshConn = atoi(argv[++ i]);
			ssCheck.str("");
			ssCheck << threshConn;
			sCheck = ssCheck.str();
			if(sCheck != argv[i])
			{
				print();
				return 0;
			}
			tagThreshConn = 1;
		}
		else if(s == "--threshSplit")
		{
			if(tagThreshSplit == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			threshSplit = atoi(argv[++ i]);
			ssCheck.str("");
			ssCheck << threshSplit;
			sCheck = ssCheck.str();
			if(sCheck != argv[i])
			{
				print();
				return 0;
			}
			tagThreshSplit = 1;
		}
                else if(s == "--closeGap")
                {
                        if(tagCloseGap == 1)
                        {
                                print();
                                return 0;
                        }
                        tagCloseGap = 1;
		}
		else if(s == "--noAlignment")
		{
			if(tagNoAlignment == 1)
			{
				print();
				return 0;
			}
			tagNoAlignment = 1;
		}
                else if(s == "--lowEukaryote")
                {
                        if(tagLowEukaryote == 1)
                        {
                                print();
                                return 0;
                        }
                        tagLowEukaryote = 1;
                }
		else if(s == "--insertLow")
		{
			if(tagInsertLow == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			insertLow = atoi(argv[++ i]);
			ssCheck.str("");
			ssCheck << insertLow;
			sCheck = ssCheck.str();
			if(sCheck != argv[i])
			{
				print();
				return 0;
			}
			tagInsertLow = 1;
		}
                else if(s == "--insertHigh")
                {
                        if(tagInsertHigh == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        insertHigh = atoi(argv[++ i]);
                        ssCheck.str("");
                        ssCheck << insertHigh;
                        sCheck = ssCheck.str();
                        if(sCheck != argv[i])
                        {
                                print();
                                return 0;
                        }
                        tagInsertHigh = 1;
                }
		else
		{
			print();
			return 0;
		}
	}

	if(tagRead1 == 0 && tagRead2 == 0 || tagRead1 == 0 && tagRead2 == 1 || tagRead1 == 1 && tagRead2 == 0 && (insertLow != 0 || insertHigh != MAX || tagCloseGap != 0) || 
	tagTransfrag == 0 || tagTranscript == 0 || insertLow > insertHigh || insertLow < 0)
	{
		print();
		return 0;
	}

	system("test -d \"tmp\"; t=$?; if [ $t -eq 1 ]; then mkdir tmp; fi");

	if(tagNoAlignment == 0)
	{
		if(tagRead1 == 1 && tagRead2 == 0)//single-end
		{
			formalizeSingleReads(r1, "tmp/_reads_1.fa", "tmp/_reads_2.fa", insertLow, insertHigh);
			PE = 0;
		}
		else//paired-end
		{
			formalizeInput(r1, "tmp/_reads_1.fa");
			formalizeInput(r2, "tmp/_reads_2.fa");
		}
		formalizeInput(t, "tmp/_transfrags.fa");
		formalizeInput(ci, "tmp/_contigs.fa");
	}

	if(tagLowEukaryote == 1)
	{
		IDENTITY = 1;
		GAPPED_IDENTITY = 1;
		HIGH_IDENTITY = 1;		
	}

	if(tagNoAlignment == 0)
	{
		parallelMap();
		cout << "(1) Alignment to de novo transfrags and contigs finished" << endl;
	}
	else
		cout << "(1) Alignment to de novo transfrags and contigs skipped" << endl;

	mapReads2Transfrags(insertLow, insertHigh);
	cout << "(2) Read alignment to de novo transfrags loaded" << endl;

	mapTransfrags2Contigs(st, threshSplit);
	cout << "(3) De novo transfrag alignment to contigs loaded" << endl;

	mapReads2Contigs(insertLow, insertHigh);
	cout << "(4) Remaining read alignment to contigs loaded" << endl;

	if(tagCloseGap)
		closeCoverageGaps(insertLow, insertHigh);
	initSegVec(threshSize, threshCov, threshSplit, insertLow, insertHigh);
	cout << "(5) Exons identified" << endl;

	iit.open("tmp/_transcripts.fa");
	for(contiID = 0, seqID = 0; contiID < contiVec.size(); contiID ++)
	{
		if(segVec[contiID].size() == 0)
			continue;
		if(segVec[contiID].size() == 1)
		{
			generateSingleTrans(iit, contiID, seqID);
			continue;
		}
		initSegJuncVec(contiID, threshConn);
		transformSegJuncVec(iit, contiID, seqID);
	}

	generateFinalTranscripts(it);
	cout << "(6) Transcripts generated" << endl;
	return 1;
}


