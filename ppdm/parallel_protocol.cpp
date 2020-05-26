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
	//cout << "DP_GSRO_Multithread" << endl;
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
	//cout << "sNodeRetrievalforDatainThread" << endl;

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
		for(int j = 0; j < dim; j++)
		{
			V2[i][j] = proto.SM_p1(V[i], cand[i][j], index);
			paillier_mul(proto.pubkey, sub_result[j], V2[i][j], sub_result[j]);
			//paillier_mul(proto.pubkey, Result[s][j], V2[cnt][j], Result[s][j]);
		}
	}
}

void sNodeRetrievalforClassificationinThread(std::vector<int> inputJ, std::vector<int *> inputNodeGroup, int &remained,	std::vector<int> inputCnt, boundary *node, paillier_ciphertext_t ***data, paillier_ciphertext_t **alpha,	paillier_ciphertext_t ***cand, protocol proto, int idx)
{
	int dim = proto.dim;
	if( proto.query == CLASSIFICATION) dim = proto.dim+1;



	paillier_ciphertext_t **tmp = new paillier_ciphertext_t *[dim];
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
				for(int k = 0; k < dim; k++)
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
				for(int k = 0; k < dim ; k++)
				{
					if(j == 1)
					{
						cand[inputCnt[i]][k] = proto.SM_p1(proto.ciper_MAX, alpha[nodeId], idx);
					}
					else
					{
						tmp[k] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
						mpz_init(tmp[k]->c);
						tmp[k] = proto.SM_p1(proto.ciper_MAX, alpha[nodeId], idx);
						paillier_mul(proto.pubkey, cand[inputCnt[i]][k], cand[inputCnt[i]][k], tmp[k]);
					}
				}
			}
		}
	}
	delete[] tmp;
}


void NodeExpansion_Multithread(paillier_ciphertext_t **q, paillier_ciphertext_t *K_DIST, std::vector<int> input, boundary *node, paillier_ciphertext_t **alpha, protocol proto, int idx)
{
	
	//std::cout << "idx = " << idx << std::endl;
	paillier_ciphertext_t **psi = new paillier_ciphertext_t *[3];
	paillier_ciphertext_t *cipher_Max = paillier_create_enc(10000);
	paillier_ciphertext_t *temp_coord1;
	paillier_ciphertext_t *temp_coord2;
	paillier_ciphertext_t *temp_coord3;
	paillier_ciphertext_t *temp_alpha = new paillier_ciphertext_t;
	mpz_init(temp_alpha->c);

	paillier_ciphertext_t **shortestPoint = new paillier_ciphertext_t *[proto.dim];
	for(int i = 0; i < proto.dim; i++)
	{
		shortestPoint[i] = new paillier_ciphertext_t;
		mpz_init(shortestPoint[i]->c);
	}

	paillier_ciphertext_t *temp_nodedist = new paillier_ciphertext_t;
	mpz_init(temp_nodedist->c);

	paillier_ciphertext_t **cihper_nodedist = new paillier_ciphertext_t*[input.size()];
	for(int i = 0; i < input.size(); i++)
	{
		cihper_nodedist[i] = new paillier_ciphertext_t;
		mpz_init(cihper_nodedist[i]->c);
	}

	for(int i = 0; i < input.size(); i++)
	{
		int num = input[i];
		for(int j = 0; j < proto.dim; j++)
		{
			psi[0] = proto.GSCMP(q[j], node[num].LL[j]);
			psi[1] = proto.GSCMP(q[j], node[num].RR[j]);
			psi[2] = proto.SBXOR(psi[0], psi[1], idx);
			temp_coord1 = proto.SM_p1(psi[2], q[j], idx);
			temp_coord2 = proto.SM_p1(psi[0], node[num].LL[j], idx);
			temp_coord3 = proto.SM_p1(proto.SBN(psi[0]), node[num].RR[j], idx);
			paillier_mul(proto.pubkey, temp_coord3, temp_coord2, temp_coord3);
			temp_coord3 = proto.SM_p1(proto.SBN(psi[2]), temp_coord3, idx);
			paillier_mul(proto.pubkey, shortestPoint[j], temp_coord1, temp_coord3);
		}
		temp_nodedist = proto.DP_SSED(q, shortestPoint, proto.dim);
		paillier_subtract(proto.pubkey, temp_alpha, proto.ciper_one, alpha[num]);
		paillier_mul(proto.pubkey, cihper_nodedist[i], proto.SM_p1(alpha[num], cipher_Max, idx), proto.SM_p1(temp_alpha, temp_nodedist, idx));
		alpha[num] = proto.GSCMP(cihper_nodedist[i], K_DIST);
	}
	delete[] psi;

	delete cipher_Max;
	delete temp_coord1;
	delete temp_coord2;
	delete temp_coord3;

	delete temp_alpha;
	delete[] shortestPoint;
	delete temp_nodedist;
	delete[] cihper_nodedist;
}

void SSED_SBD_inThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t **query, paillier_ciphertext_t ***data,
	int dim, paillier_ciphertext_t **cipher_distance, paillier_ciphertext_t ***cipher_SBD_distance, protocol proto)
{
	int len = inputIdx.size();
	
	for(int i = 0; i < len; i++)
	{
		cipher_distance[inputIdx[i]] = proto.SSEDm(query, data[inputIdx[i]], proto.dim);
		cipher_SBD_distance[inputIdx[i]] = proto.SBD(cipher_distance[inputIdx[i]]);
		/*
		paillier_print("dist : ", cipher_distance[inputIdx[i]]);	
		for(int j = 0; j < proto.size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, proto.pubkey, proto.prvkey, cipher_SBD_distance[inputIdx[i]][j]));
		}
		printf("\n");
		*/
	}
}
void SSED_inThread(std::vector<int> inputIdx, paillier_ciphertext_t **query, paillier_ciphertext_t ***data, paillier_ciphertext_t **cipher_distance,
	protocol proto)
{
	int len = inputIdx.size();
	
	for(int i = 0; i < len; i++)
	{
		cipher_distance[inputIdx[i]] = proto.SSEDm(query, data[inputIdx[i]], proto.dim);
		//cipher_SBD_distance[inputIdx[i]] = proto.SBD(cipher_distance[inputIdx[i]]);
		/*
		paillier_print("dist : ", cipher_distance[inputIdx[i]]);	
		for(int j = 0; j < proto.size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, proto.pubkey, proto.prvkey, cipher_SBD_distance[inputIdx[i]][j]));
		}
		printf("\n");
		*/
		
	}
}

void SBD_inThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t **query, paillier_ciphertext_t ***data, paillier_ciphertext_t **cipher_distance,
	paillier_ciphertext_t ***cipher_SBD_distance, protocol proto)
{
	int len = inputIdx.size();
	
	for(int i = 0; i < len; i++)
	{
		//cipher_distance[inputIdx[i]] = proto.SSEDm(query, data[inputIdx[i]], dim);
		cipher_SBD_distance[inputIdx[i]] = proto.SBD(cipher_distance[inputIdx[i]]);
		/*
		paillier_print("dist : ", cipher_distance[inputIdx[i]]);	
		for(int j = 0; j < proto.size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, proto.pubkey, proto.prvkey, cipher_SBD_distance[inputIdx[i]][j]));
		}
		printf("\n");
		*/
	}
}
void Smin_n_InThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t ***copy_cipher, paillier_ciphertext_t ***min, protocol proto)
{
	paillier_ciphertext_t ***ciper_SBD_distance = new paillier_ciphertext_t **[inputIdx.size()];
	//paillier_ciphertext_t **smin;
	//std::cout << "i = " << inputIdx.size() << std::endl;
	for(int i = 0; i < inputIdx.size(); i++)
	{
		ciper_SBD_distance[i] = copy_cipher[inputIdx[i]];
	}
	/*
	for(int i = 0; i < inputIdx.size(); i++)
	{
		std::cout << "i = " << i;
		for(int j = 0; j < proto.size; j++ ) {
			gmp_printf(" %Zd", paillier_dec(0, proto.pubkey, proto.prvkey, ciper_SBD_distance[i][j]));
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/
	min[index] = proto.Smin_n(ciper_SBD_distance, inputIdx.size());
	/*
	for(int i = 0; i < proto.size; i++)
	{
		gmp_printf("%Zd ", paillier_dec(0, proto.pubkey, proto.prvkey, smin[i]));
		
	}
	std::cout << std::endl;
	*/
	delete[] ciper_SBD_distance;

	//min[index] = smin;
	/*
	for(int i = 0; i < proto.size; i++)
	{
		gmp_printf("%Zd ", paillier_dec(0, proto.pubkey, proto.prvkey, min[index][i]));
		
	}
	std::cout << std::endl;
	*/
}

void caldistanceInThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t ***cipher_SBD_distance, paillier_ciphertext_t **cipher_distance,
	protocol proto)
{
	paillier_ciphertext_t *cipher_dist	= paillier_create_enc_zero();
	int len = inputIdx.size();

	for(int i = 0; i < len; i++)
	{
		for(int j = proto.size; j > 0; j--)
		{
			int t = (int)pow(2, j - 1);
			paillier_ciphertext_t *cipher_binary = paillier_create_enc(t);
			cipher_binary = proto.SM_p1(cipher_binary, cipher_SBD_distance[inputIdx[i]][proto.size - j]);
			paillier_mul(proto.pubkey, cipher_dist, cipher_binary, cipher_dist);
		}
		cipher_distance[inputIdx[i]] = cipher_dist;
		cipher_dist = paillier_enc(0, proto.pubkey, proto.plain_zero, paillier_get_rand_devurandom);
	}
	/*
	if( s != 0 )
	{
		for(int i = 0; i < n; i++ )
		{
			for(int j = size; j > 0; j--)
			{
				t = (int)pow(2, j-1);
				ciper_binary = paillier_create_enc(t);
				ciper_binary = SM_p1(ciper_binary, ciper_SBD_distance[i][size-j]);
				paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
			}
			ciper_distance[i] = ciper_dist;
			ciper_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
		}
	}
	*/
}

void classipmInThread1(int index, std::vector<int> inputIdx, paillier_ciphertext_t **cipher_distance, paillier_ciphertext_t *cipher_min,
	paillier_ciphertext_t *cipher_rand, paillier_ciphertext_t **cipher_mid, protocol proto)
{
	paillier_ciphertext_t *temp_dist = paillier_create_enc_zero();
	int len = inputIdx.size();
	for(int i = 0; i < len; i++)
	{
		paillier_subtract(proto.pubkey, temp_dist, cipher_distance[inputIdx[i]], cipher_min);
		cipher_mid[inputIdx[i]] = proto.SM_p1(temp_dist, cipher_rand, index);
		//paillier_print("dist - dist : ", cipher_mid[inputIdx[i]]);
	}
}

void classipmInThread2(int index, std::vector<int> inputIdx, paillier_ciphertext_t **cipher_V, paillier_ciphertext_t **cipher_V2,
	paillier_ciphertext_t ***data, paillier_ciphertext_t **cipher_rabel_result, int s, protocol proto)
{
	paillier_ciphertext_t *temp_dist = paillier_create_enc_zero();
	int len = inputIdx.size();
	for(int i = 0; i < len; i++)
	{
		cipher_V2[inputIdx[i]] = proto.SM_p1(cipher_V[inputIdx[i]], data[inputIdx[i]][proto.dim]);
		//paillier_print("ciper_V2 : ", cipher_V2[inputIdx[i]]);
		paillier_mul(proto.pubkey, cipher_rabel_result[s], cipher_V2[inputIdx[i]], cipher_rabel_result[s]);
	}
}

void SBOR_inThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t **cipher_V, paillier_ciphertext_t ***cipher_SBD_distance, protocol proto)
{
	int len = inputIdx.size();
	
	for(int i = 0; i < len; i++)
	{
		for(int j = 0; j < proto.size; j++)
		{
			cipher_SBD_distance[inputIdx[i]][j] = proto.SBOR(cipher_V[inputIdx[i]], cipher_SBD_distance[inputIdx[i]][j]);
		}
	}
}

void Smin_basic1_InThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2, bool func,
	paillier_ciphertext_t* ciper_Rand_value, paillier_plaintext_t* Rand_value,
	paillier_ciphertext_t** ciper_W, paillier_ciphertext_t** ciper_R, paillier_ciphertext_t** ciper_H, paillier_ciphertext_t** ciper_L,
	paillier_ciphertext_t** ciper_M, paillier_ciphertext_t** ciper_lambda, paillier_ciphertext_t** ciper_min, paillier_ciphertext_t** ciper_O,
	paillier_ciphertext_t** ciper_G, paillier_ciphertext_t** temp, paillier_ciphertext_t** temp2, paillier_ciphertext_t** temp3, protocol proto)
{
	int len = inputIdx.size();

	for(int i = 0; i < len; i++)
	{
		if(func){	// true :  F : u>v	
			paillier_subtract(proto.pubkey, ciper_W[inputIdx[i]], ciper1[inputIdx[i]], proto.SM_p1(ciper1[inputIdx[i]], ciper2[inputIdx[i]], index));	// W
			paillier_subtract(proto.pubkey, temp[inputIdx[i]], ciper2[inputIdx[i]], ciper1[inputIdx[i]]);	
			paillier_mul(proto.pubkey, ciper_R[inputIdx[i]], temp[inputIdx[i]], ciper_Rand_value);	// Gamma
		}else{
			paillier_subtract(proto.pubkey, ciper_W[inputIdx[i]], ciper2[inputIdx[i]], proto.SM_p1(ciper1[inputIdx[i]], ciper2[inputIdx[i]], index));
			paillier_subtract(proto.pubkey, temp[inputIdx[i]], ciper1[inputIdx[i]], ciper2[inputIdx[i]]);
			paillier_mul(proto.pubkey, ciper_R[inputIdx[i]], temp[inputIdx[i]], ciper_Rand_value);
		}
		ciper_G[inputIdx[i]] = proto.SBXOR(ciper1[inputIdx[i]], ciper2[inputIdx[i]]);

		if(inputIdx[i]==0){
			paillier_exp(proto.pubkey,ciper_H[inputIdx[i]], proto.ciper_zero, Rand_value);
			paillier_mul(proto.pubkey,ciper_H[inputIdx[i]], ciper_H[inputIdx[i]], ciper_G[inputIdx[i]]);		
		}else{
			paillier_exp(proto.pubkey,temp[inputIdx[i]], ciper_H[inputIdx[i]-1], Rand_value);
			paillier_mul(proto.pubkey,ciper_H[inputIdx[i]], temp[inputIdx[i]], ciper_G[inputIdx[i]]);
		}

		paillier_mul(proto.pubkey, ciper_O[inputIdx[i]], ciper_H[inputIdx[i]], proto.ciper_minus);	// PI
		paillier_exp(proto.pubkey, ciper_O[inputIdx[i]], ciper_O[inputIdx[i]], Rand_value);
		
		paillier_mul(proto.pubkey, ciper_L[inputIdx[i]], ciper_O[inputIdx[i]], ciper_W[inputIdx[i]]);
	}
}