#include <iostream>
#include <fstream>
#include <gmp.h>
#include <string>
#include <sstream>
#include <thread>
#include <chrono>
#include <vector>
#include <mutex>
#include <iostream>
#include "parallel_protocol.h"

using namespace std;

static mutex mtx;


// PARALLEL SBD_SRO 2020.05.27
void SRO_SBD_inThread(std::vector<int> input, paillier_ciphertext_t *** ciper_qLL_bit, paillier_ciphertext_t *** ciper_qRR_bit, paillier_ciphertext_t *** data, paillier_ciphertext_t ** alpha, protocol proto)
{
	cout << "SRO_SBD_inThread" << endl;
	paillier_ciphertext_t*** candLL_bit = new paillier_ciphertext_t**[proto.dim];
	paillier_ciphertext_t*** candRR_bit = new paillier_ciphertext_t**[proto.dim];

	for( int i = 0 ; i < input.size() ; i++ ) 
	{
		int num = input[i];
		for( int j = 0 ; j < proto.dim ; j++ ) 
		{
			candLL_bit[j] = proto.SBD_for_SRO(data[num][j], 0);			
			candRR_bit[j] = proto.SBD_for_SRO(data[num][j], 1);		
		}		
		alpha[num] = proto.SRO(candLL_bit, candRR_bit, ciper_qLL_bit, ciper_qRR_bit);
	}
	
	delete [] candLL_bit;
	delete [] candRR_bit;
}

// PARALLEL SBD_COMPUTE_SCORE
void CS_SBD_inThread(std::vector<int>input, paillier_ciphertext_t ** query, paillier_ciphertext_t *** data, paillier_ciphertext_t ** cipher_distance, paillier_ciphertext_t *** cipher_SBD_distance, protocol proto)
{
	for( int i = 0 ; i < input.size() ; i++ ){
		int num = input[i];
		cipher_distance[num] = proto.computeScore(query, data[num]);
		cipher_SBD_distance[num] = proto.SBD(cipher_distance[num]);
	}
}

void Smax_n_InThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t ***copy_cipher, paillier_ciphertext_t ***min, protocol proto)
{
	paillier_ciphertext_t ***ciper_SBD_distance = new paillier_ciphertext_t **[inputIdx.size()];
	for(int i = 0; i < inputIdx.size(); i++)
	{
		ciper_SBD_distance[i] = copy_cipher[inputIdx[i]];
	}
	min[index] = proto.Smax_n(ciper_SBD_distance, inputIdx.size());
	delete[] ciper_SBD_distance;
}

void Smax_basic1_InThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2, bool func, paillier_ciphertext_t* ciper_Rand_value, paillier_plaintext_t* Rand_value, paillier_ciphertext_t** ciper_W, paillier_ciphertext_t** ciper_R, paillier_ciphertext_t** ciper_H, paillier_ciphertext_t** ciper_L,	paillier_ciphertext_t** ciper_M, paillier_ciphertext_t** ciper_lambda, paillier_ciphertext_t** ciper_min, paillier_ciphertext_t** ciper_O, paillier_ciphertext_t** ciper_G, paillier_ciphertext_t** temp, paillier_ciphertext_t** temp2, paillier_ciphertext_t** temp3, protocol proto)
{
	
}


void recalculate_dist_inThread(std::vector<int>input, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_distance, protocol proto)
{
	int t = 0;
	paillier_ciphertext_t * cipher_binary = new paillier_ciphertext_t;
	mpz_init(cipher_binary->c);
	paillier_ciphertext_t * cipher_dist = paillier_create_enc(0);

	for( int i = 0 ; i < input.size() ; i++ ){
		int num = input[i];
		for( int j = proto.size ; j > 0 ; j--){
			t = (int)pow(2, j-1);
			cipher_binary = paillier_create_enc(t);
			cipher_binary = proto.SM_p1(cipher_binary, cipher_SBD_distance[num][(proto.size)-j]);
			paillier_mul(proto.pubkey, cipher_dist, cipher_binary, cipher_dist);
		}
		cipher_distance[num] = cipher_dist;
		cipher_dist = paillier_create_enc(0);
	}
}


void DuringProcessedInTOPK_PB_inThread( std::vector<int>input, paillier_ciphertext_t** cipher_distance, paillier_ciphertext_t** cipher_mid, paillier_ciphertext_t* cipher_min, paillier_ciphertext_t* cipher_rand, protocol proto)
{
	paillier_ciphertext_t * temp_dist = new paillier_ciphertext_t;
	mpz_init(temp_dist->c);
	
	for( int i = 0 ; i < input.size() ; i++ ){
		int num = input[i];
		paillier_subtract(proto.pubkey, temp_dist, cipher_distance[num], cipher_min);
		cipher_mid[num] = proto.SM_p1(temp_dist, cipher_rand);
	}
}

void ExtractTOPKInTOPK_PB_inThread( std::vector<int>input, paillier_ciphertext_t*** data, paillier_ciphertext_t** cipher_V, paillier_ciphertext_t ** cipher_thread_result, protocol proto)
{
	paillier_ciphertext_t * temp_value = new paillier_ciphertext_t;
	mpz_init(temp_value->c);

	for( int i = 0 ; i < input.size() ; i++ ){
		int num = input[i];
		for( int j = 0 ; j < proto.dim; j++ ){
			temp_value = proto.SM_p1(cipher_V[num], data[num][j]);
			paillier_mul(proto.pubkey, cipher_thread_result[j], temp_value, cipher_thread_result[j]);
		}
	}
}

void UPDATE_SBD_SCORE_InTOPK_PB_inThread(std::vector<int>input, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_V, protocol proto)
{
	//cout << "UPDATE_SBD_SCORE_InTOPK_PB_inThread" << endl;
	for( int i = 0 ; i < input.size() ; i++ ){
		int num = input[i];
		for( int j = 0 ; j < proto.size ; j++ ){
			cipher_SBD_distance[num][j] = proto.SM_p1(proto.SBN(cipher_V[num]), cipher_SBD_distance[num][j]);
		}
	}
}



void recalculate_DISTforkNN_inThread(std::vector<int> input, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_distance, protocol proto)
{
	int t = 0;
	paillier_ciphertext_t * cipher_binary = new paillier_ciphertext_t;
	mpz_init(cipher_binary->c);
	paillier_ciphertext_t * cipher_dist = paillier_create_enc(0);

	for( int i = 0 ; i < input.size() ; i++ ){
		int num = input[i];
		for( int j = proto.size ; j > 0 ; j--){
			t = (int)pow(2, j-1);
			cipher_binary = paillier_create_enc(t);
			cipher_binary = proto.SM_p1(cipher_binary, cipher_SBD_distance[num][(proto.size)-j]);
			paillier_mul(proto.pubkey, cipher_dist, cipher_binary, cipher_dist);
		}
		cipher_distance[num] = cipher_dist;
		cipher_dist = paillier_create_enc(0);
	}
}

void UPDATE_SBD_SCORE_InKNN_PB_inThread(std::vector<int>input, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_V, protocol proto)
{
	//cout << "UPDATE_SBD_SCORE_InKNN_PB_inThread" << endl;
	for( int i = 0 ; i < input.size() ; i++ ){
		int num = input[i];
		for( int j = 0 ; j < proto.size ; j++ ){
			cipher_SBD_distance[num][j] = proto.SBOR(cipher_V[num], cipher_SBD_distance[num][j]);
		}
	}
}













void GSCMPforTopk_inThread(std::vector<int> input, boundary* node, paillier_ciphertext_t** psi, paillier_ciphertext_t* kth_score, paillier_ciphertext_t** q, paillier_ciphertext_t** alpha, protocol proto)
{
	cout << "GSCMPforTopk_inThread" << endl;
	paillier_ciphertext_t** ShortestPoint = new paillier_ciphertext_t*[proto.dim];
	for ( int i = 0 ; i < proto.dim ; i++ )
	{
		ShortestPoint[i] = new paillier_ciphertext_t;
		mpz_init(ShortestPoint[i]->c);
	}

	for( int i = 0 ; i < input.size() ; i++ )
	{	
		int num = input[i];
		for( int j = 0 ; j < proto.dim ; j++ )
		{
			paillier_mul(proto.pubkey, ShortestPoint[j], proto.SM_p1(node[num].RR[j], psi[j]), proto.SM_p1(node[num].LL[j], proto.SBN(psi[j])));	
		}
		alpha[num] = proto.GSCMP(kth_score , proto.computeScore(q, ShortestPoint));	// k번째 score보다 높을 가능성이 있으면 alpha=1
	}
}

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
	}
}
void SSED_inThread(std::vector<int> inputIdx, paillier_ciphertext_t **query, paillier_ciphertext_t ***data, paillier_ciphertext_t **cipher_distance,
	protocol proto)
{
	int len = inputIdx.size();
	
	for(int i = 0; i < len; i++)
	{
		cipher_distance[inputIdx[i]] = proto.SSEDm(query, data[inputIdx[i]], proto.dim);	
	}
}

void SBD_inThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t **query, paillier_ciphertext_t ***data, paillier_ciphertext_t **cipher_distance,
	paillier_ciphertext_t ***cipher_SBD_distance, protocol proto)
{
	int len = inputIdx.size();
	
	for(int i = 0; i < len; i++)
	{
		cipher_SBD_distance[inputIdx[i]] = proto.SBD(cipher_distance[inputIdx[i]]);
	}
}
void Smin_n_InThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t ***copy_cipher, paillier_ciphertext_t ***min, protocol proto)
{
	paillier_ciphertext_t ***ciper_SBD_distance = new paillier_ciphertext_t **[inputIdx.size()];
	for(int i = 0; i < inputIdx.size(); i++)
	{
		ciper_SBD_distance[i] = copy_cipher[inputIdx[i]];
	}	
	min[index] = proto.Smin_n(ciper_SBD_distance, inputIdx.size());
	delete[] ciper_SBD_distance;
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



void Comp_Cluster_inThread(paillier_ciphertext_t*** cipher, paillier_ciphertext_t*** former_Center, std::vector<int> inputIdx, paillier_ciphertext_t*** NewSumCluster, paillier_ciphertext_t** NewSumCntCluster, protocol proto, int i)
{
	paillier_ciphertext_t** Distance_center_data = new paillier_ciphertext_t*[proto.k];
	paillier_ciphertext_t** idx_arr = new paillier_ciphertext_t*[proto.k];
	paillier_ciphertext_t* tmp_data = new paillier_ciphertext_t;

	for( int j = 0 ; j < proto.k ; j++){
		Distance_center_data[j] = new paillier_ciphertext_t;
		idx_arr[j] = new paillier_ciphertext_t;

		mpz_init(Distance_center_data[j]->c);
		mpz_init(idx_arr[j]->c);
	}
	mpz_init(tmp_data->c);

	for(int i = 0 ; i < inputIdx.size() ; i++){
		int num = inputIdx[i];
		for( int j = 0 ; j < proto.k ; j++){
			Distance_center_data[j] = proto.SSEDm(cipher[num], former_Center[j], proto.dim);
		}
		
		idx_arr = proto.Smin_bool(Distance_center_data, proto.k);

		mtx.lock();
		for( int j = 0 ; j < proto.k ; j++ ){
			paillier_mul(proto.pubkey, NewSumCntCluster[j], idx_arr[j], NewSumCntCluster[j]);
			for( int e = 0 ; e < proto.dim ; e++ ){
				tmp_data = proto.SM_p1(idx_arr[j], cipher[num][e]);
				paillier_mul(proto.pubkey, NewSumCluster[j][e], NewSumCluster[j][e], tmp_data);
			}

		}
		mtx.unlock();
	}
	delete [] Distance_center_data;
	delete [] idx_arr;
	delete [] tmp_data;
}