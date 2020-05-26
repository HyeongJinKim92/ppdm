#include <iostream>
#include <fstream>
#include <gmp.h>
#include <string>
#include <sstream>
#include <thread>
#include <chrono>
#include <vector>
#include <iostream>
#include "parallel_protocol.h"

using namespace std;


void DP_GSRO_MultithreadforTopk(paillier_ciphertext_t **qLL, paillier_ciphertext_t **qRR, std::vector<int> input, boundary *node, paillier_ciphertext_t **alpha, protocol proto, bool verify)
{
	cout << "DP_GSRO_Multithread" << endl;
	for(int i = 0; i < input.size(); i++)
	{
		int num = input[i];
		alpha[num] = proto.DP_GSRO_inThread(qLL, qRR, node[num].LL, node[num].RR, 0);
		for( int j = 0 ; j < proto.dim ; j++ )
		{
			node[num].LL[j] = proto.SM_p1(node[num].LL[j], proto.SBN(alpha[num]));
			node[num].RR[j] = proto.SM_p1(node[num].RR[j], proto.SBN(alpha[num]));
		}
	}
}



void DP_GSRO_Multithread(paillier_ciphertext_t **qLL, paillier_ciphertext_t **qRR, std::vector<int> input, boundary *node, paillier_ciphertext_t **alpha, protocol proto, int idx)
{
	//cout << "DP_GSRO_Multithread" << endl;
	for(int i = 0; i < input.size(); i++)
	{
		int num = input[i];
		alpha[num] = proto.DP_GSRO_inThread(qLL, qRR, node[num].LL, node[num].RR, 0);
	}
}

void DP_GSRO_dataquery_Multithread(paillier_ciphertext_t *** cand, boundary q, std::vector<int> input, paillier_ciphertext_t **alpha, protocol proto, int idx, paillier_ciphertext_t * cipher_rand)
{
	for(int i = 0; i < input.size(); i++)
	{
		int num = input[i];
		alpha[num] = proto.DP_GSRO_inThread(cand[num], cand[num], q.LL, q.RR, idx);
		for( int j = 0 ; j < proto.dim ; j++ ){
			paillier_mul( proto.pubkey, cand[num][j], cand[num][j], cipher_rand);
		}
	}
}

void sNodeRetrievalforDatainThread(std::vector<int> inputJ, std::vector<int *> inputNodeGroup, int &remained,	std::vector<int> inputCnt, boundary *node, paillier_ciphertext_t ***data, paillier_ciphertext_t **alpha, paillier_ciphertext_t ***cand, protocol proto, int idx)
{
	cout << "sNodeRetrievalforDatainThread" << endl;

	int end_dim = proto.dim;
	paillier_ciphertext_t * TMP_VALUE = new paillier_ciphertext_t;
	TMP_VALUE = proto.ciper_MAX;

	if (proto.query == RANGE){}
	else if (proto.query == TOPK){TMP_VALUE = proto.ciper_zero;}
	else if (proto.query == KNN){}
	else if(proto.query == CLASSIFICATION){end_dim = proto.dim+1;}

	paillier_ciphertext_t **tmp = new paillier_ciphertext_t *[end_dim];
	for(int i = 0; i < inputJ.size(); i++)
	{
		if(remained == 0)
			break;

		for(int j = 1; j <= inputNodeGroup[i][0]; j++)
		{
			int nodeId = inputNodeGroup[i][j];
			if(node[nodeId].NumData >= inputJ[i] + 1)
			{
				int dataId = node[nodeId].indata_id[inputJ[i]];
				for(int k = 0; k < end_dim ; k++)
				{
					if(j == 1)
					{
						cand[inputCnt[i]][k] = proto.SM_p1(data[dataId][k], alpha[nodeId], idx);
					}
					else
					{
						tmp[k] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
						mpz_init(tmp[k]->c);
						tmp[k] = proto.SM_p1(data[dataId][k], alpha[nodeId], idx);
						paillier_mul(proto.pubkey, cand[inputCnt[i]][k], cand[inputCnt[i]][k], tmp[k]);
					}
				}
				if(node[nodeId].NumData == j + 1)
					remained--;
			}
			else
			{
				for(int k = 0; k < end_dim ; k++)
				{
					if(j == 1)
					{
						cand[inputCnt[i]][k] = proto.SM_p1(TMP_VALUE, alpha[nodeId], idx);
					}
					else
					{
						tmp[k] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
						mpz_init(tmp[k]->c);
						tmp[k] = proto.SM_p1(TMP_VALUE, alpha[nodeId], idx);
						paillier_mul(proto.pubkey, cand[inputCnt[i]][k], cand[inputCnt[i]][k], tmp[k]);
					}
					
				}
			}
		}
	}
	delete[] tmp;
}

void ComputeScore_Multithread(paillier_ciphertext_t *** cand, paillier_ciphertext_t ** q, paillier_ciphertext_t ** SCORE, std::vector<int> input, protocol proto, bool idx)
{
	for( int i = 0 ; i < input.size() ; i++ )
	{
		int num = input[i];
		SCORE[num] = proto.computeScore(q, cand[num]);
	}
}

void MAXnInthread2(std::vector<int> inputCnt, paillier_ciphertext_t **SCORE_MINUS_MAX, paillier_ciphertext_t **SCORE,	paillier_ciphertext_t *MAX, paillier_ciphertext_t *C_RAND, protocol proto, int idx)
{
	for(int i = 0; i < inputCnt.size(); i++)
	{
		int num = inputCnt[i];
		paillier_subtract(proto.pubkey, SCORE_MINUS_MAX[num], SCORE[num], MAX);
		SCORE_MINUS_MAX[num] = proto.SM_p1(SCORE_MINUS_MAX[num], C_RAND, idx);
	}
}

void MAXnInThread3_1(int index, int s, std::vector<int> inputCnt, paillier_ciphertext_t **V, paillier_ciphertext_t **SCORE,	paillier_ciphertext_t ***cand, paillier_ciphertext_t *MIN, protocol proto)
{
	paillier_ciphertext_t *TMP_alpha = new paillier_ciphertext_t;
	mpz_init(TMP_alpha->c);

	for(int i = 0; i < inputCnt.size(); i++)
	{
		int c = inputCnt[i];
		//std::cout << c << std::endl;
		paillier_subtract(proto.pubkey, TMP_alpha, proto.ciper_one, V[c]);
		paillier_mul(proto.pubkey, SCORE[c], proto.SM_p1(V[c], MIN), proto.SM_p1(TMP_alpha, SCORE[c]));
	}
	delete TMP_alpha;
}

void MAXnInThread3_2_v2(int index, int offset, int record_size, int dim, paillier_ciphertext_t** sub_result, paillier_ciphertext_t **V, paillier_ciphertext_t ***V2, paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, protocol proto)
{
	int size = offset+record_size;
	for(int i = offset; i < size; i++ )
	{
		for(int j = 0; j < proto.dim; j++)
		{
			V2[i][j] = proto.SM_p1(V[i], cand[i][j]);
			paillier_mul(proto.pubkey, sub_result[j], V2[i][j], sub_result[j]);
		}
	}
}


void SMSnInthread2(std::vector<int> inputCnt, paillier_ciphertext_t **DIST_MINUS_MIN, paillier_ciphertext_t **DIST,	paillier_ciphertext_t *MIN, paillier_ciphertext_t *C_RAND, protocol proto, int idx)
{
	for(int i = 0; i < inputCnt.size(); i++)
	{
		int num = inputCnt[i];
		paillier_subtract(proto.pubkey, DIST_MINUS_MIN[num], DIST[num], MIN);
		DIST_MINUS_MIN[num] = proto.SM_p1(DIST_MINUS_MIN[num], C_RAND, idx);
	}
}

void SSEDInThread(std::vector<int> inputCnt, paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t **q, int cnt, protocol proto)
{
	for(int i = 0; i < inputCnt.size(); i++)
	{
		int num = inputCnt[i];
		DIST[num] = proto.DP_SSED(cand[num], q, proto.dim);
	}
}


void SMSnInThread3_1(int index, int s, std::vector<int> inputCnt, paillier_ciphertext_t **V, paillier_ciphertext_t **DIST,	paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, protocol proto)
{
	paillier_ciphertext_t *TMP_alpha = new paillier_ciphertext_t;
	mpz_init(TMP_alpha->c);

	for(int i = 0; i < inputCnt.size(); i++)
	{
		int c = inputCnt[i];
		//std::cout << c << std::endl;
		paillier_subtract(proto.pubkey, TMP_alpha, proto.ciper_one, V[c]);
		paillier_mul(proto.pubkey, DIST[c], proto.SM_p1(V[c], MAX, index), proto.SM_p1(TMP_alpha, DIST[c], index));
	}

	delete TMP_alpha;
}


void SMSnInThread3_2_v2(int index, int offset, int record_size, int dim, paillier_ciphertext_t** sub_result, paillier_ciphertext_t **V, paillier_ciphertext_t ***V2, paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, protocol proto)
{
	int size = offset+record_size;
	for(int i = offset; i < size; i++ )
	{
		for(int j = 0; j < dim+1; j++)
		{
			V2[i][j] = proto.SM_p1(V[i], cand[i][j], index);
			paillier_mul(proto.pubkey, sub_result[j], V2[i][j], sub_result[j]);
			//paillier_mul(proto.pubkey, Result[s][j], V2[cnt][j], Result[s][j]);
		}
	}
}

void sNodeRetrievalforClassificationinThread(std::vector<int> inputJ, std::vector<int *> inputNodeGroup, int &remained,	std::vector<int> inputCnt, boundary *node, paillier_ciphertext_t ***data, paillier_ciphertext_t **alpha,	paillier_ciphertext_t ***cand, protocol proto, int idx)
{
	paillier_ciphertext_t **tmp = new paillier_ciphertext_t *[proto.dim + 1];
	//std::cout << "inputJ.size() = " << inputJ.size() << std::endl;
	for(int i = 0; i < inputJ.size(); i++)
	{
		//std::cout << index << "th i = " << i << std::endl;

		if(remained == 0)
			break;

		for(int j = 1; j <= inputNodeGroup[i][0]; j++)
		{
			int nodeId = inputNodeGroup[i][j];
			//printf("selected node ID : %d\n", nodeId);

			if(node[nodeId].NumData >= inputJ[i] + 1)
			{
				int dataId = node[nodeId].indata_id[inputJ[i]];
				//printf("selected data ID : %d\n", dataId);

				//paillier_ciphertext_t** tmp = new paillier_ciphertext_t *[proto.dim + 1];
				for(int k = 0; k < proto.dim + 1; k++)
				{
					if(j == 1)
					{
						//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey,data[dataId][k]), paillier_dec(0, proto.pubkey, proto.prvkey, alpha[nodeId]));
						//printf("cnt : %d,  dim : %d\n", inputCnt[i], k);
						//cand[inputCnt[i]][k] = proto.SM_p1_multiThread(data[dataId][k], alpha[nodeId]);
						cand[inputCnt[i]][k] = proto.SM_p1(data[dataId][k], alpha[nodeId], idx);
						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, cand[inputCnt[i]][k]));
						//std::cout << "aaaa " << index << std::endl;
					}
					else
					{
						//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey,data[dataId][k]), paillier_dec(0, proto.pubkey, proto.prvkey, alpha[nodeId]));
						//printf("cnt : %d,  dim : %d\n", inputCnt[i], k);
						
						tmp[k] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
						mpz_init(tmp[k]->c);
						//tmp[k] = proto.SM_p1_multiThread(data[dataId][k], alpha[nodeId]);
						tmp[k] = proto.SM_p1(data[dataId][k], alpha[nodeId], idx);
						paillier_mul(proto.pubkey, cand[inputCnt[i]][k], cand[inputCnt[i]][k], tmp[k]);

						//paillier_ciphertext_t *tmp = new paillier_ciphertext_t;
						//mpz_init(tmp->c);
						//tmp = proto.SM_p1_multiThread(data[dataId][k], alpha[nodeId]);

						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, cand[inputCnt[i]][k]));
						
						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, tmp[k]));

						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, tmp));
						//paillier_mul(proto.pubkey, cand[inputCnt[i]][k], cand[inputCnt[i]][k], tmp);

						//gmp_printf("cnt : %d, (%Zd)\n", inputCnt[i], paillier_dec(0, proto.pubkey, proto.prvkey, cand[inputCnt[i]][k]));
					}
				}

				if(node[nodeId].NumData == j + 1)
					remained--;
			}
			else
			{
				for(int k = 0; k < proto.dim + 1; k++)
				{
					if(j == 1)
					{
						//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, proto.ciper_MAX), paillier_dec(0, proto.pubkey, proto.prvkey, alpha[nodeId]));
						//printf("cnt : %d,  dim : %d\n", inputCnt[i], k);
						//cand[inputCnt[i]][k] = proto.SM_p1_multiThread(proto.ciper_MAX, alpha[nodeId]);
						cand[inputCnt[i]][k] = proto.SM_p1(proto.ciper_MAX, alpha[nodeId], idx);
						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, cand[inputCnt[i]][k]));
					}
					else
					{
						//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, proto.ciper_MAX), paillier_dec(0, proto.pubkey, proto.prvkey, alpha[nodeId]));
						//printf("cnt : %d,  dim : %d\n", inputCnt[i], k);

						tmp[k] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
						mpz_init(tmp[k]->c);
						//tmp[k] = proto.SM_p1_multiThread(proto.ciper_MAX, alpha[nodeId]);
						tmp[k] = proto.SM_p1(proto.ciper_MAX, alpha[nodeId], idx);
						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, cand[inputCnt[i]][k]));
						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, tmp[k]));

						//paillier_ciphertext_t *tmp = new paillier_ciphertext_t;
						//mpz_init(tmp->c);
						//tmp = proto.SM_p1_multiThread(proto.ciper_MAX, alpha[nodeId]);
						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, cand[inputCnt[i]][k]));
						//gmp_printf("%Zd \n", paillier_dec(0, proto.pubkey, proto.prvkey, tmp));
						
						paillier_mul(proto.pubkey, cand[inputCnt[i]][k], cand[inputCnt[i]][k], tmp[k]);
						//paillier_mul(proto.pubkey, cand[inputCnt[i]][k], cand[inputCnt[i]][k], tmp);
						//gmp_printf("(%Zd)\n", paillier_dec(0, proto.pubkey, proto.prvkey, cand[inputCnt[i]][k]));
					}
					
				}
			}
		}
		//*(cnt)++;
	}
	delete[] tmp;
}