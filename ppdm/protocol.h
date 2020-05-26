#ifndef _PROTOCOL_H_
#define	 _PROTOCOL_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <string>
#include <sstream>
#include <thread>
#include <chrono>
#include <vector>
#include <mutex>
#include <iostream>
#include <math.h>
#include <time.h>
#include <gmp.h>

#include "paillier.h"
#include "types.h"


/*
#define BOOL int
#define true 1
#define false 0
*/
class protocol
{

	private :
		static std::mutex mtx;

	public :
		paillier_pubkey_t* pubkey;
		paillier_prvkey_t* prvkey;

		paillier_plaintext_t* plain_zero;
		paillier_plaintext_t* plain_one;
		paillier_plaintext_t* plain_two;
		paillier_plaintext_t* plain_minus;
		paillier_plaintext_t* plain_minus_two;
		paillier_plaintext_t* plain_MAX;
		paillier_plaintext_t* rand1;
		paillier_plaintext_t* rand2;	



		paillier_ciphertext_t* ciper_zero;
		paillier_ciphertext_t* ciper_one;
		paillier_ciphertext_t* ciper_minus;
		paillier_ciphertext_t* ciper_MAX;
		paillier_ciphertext_t* l;		
		paillier_ciphertext_t** temp;
		paillier_ciphertext_t* result;
		

	// generate noise (random numbers)

	// encrypt noise 
		paillier_ciphertext_t* ciper_rand1;
		paillier_ciphertext_t* ciper_rand2;




		bool Print;
		
		//for SMSn 
		paillier_plaintext_t* SMSn_value;
		bool smsn_flag;
		int smsn_cnt;
		int smsn_result_idx;
		int smsn_tmp_val;
		int smsn_tmp_rand;
		//for MAXn
		bool MAXn_flag;
		int MAXn_cnt;
		int MAXn_result_idx;
		int MAXn_tmp_val;
		int MAXn_tmp_rand;
		// parameters
		QUERY query;
		RANGE_TYPE rquery_t;
		TOPK_TYPE tquery_t;
		KNN_TYPE kquery_t;
		CLASSIFICATION_TYPE cquery_t;
		KMEANS_TYPE mquery_t;		
		int modul;
		int NumData;
		int tree_level;
		int thread_num;
		int dim;
		int size;
		int DP;
		int k;
		int FanOut;
		int totalNumOfRetrievedNodes;		// 제안하는 기법에서 탐색한 노드의 총 수를 의미

		// time variables
		float total_time;
		float node_SBD_time;
		float node_SRO_time;
		float node_expansion_time;
		float data_extract_first_time;
		float data_extract_second_time;
		float data_SSED_SBD_time;
		float sMINn_first_time;
		float sMINn_second_time;
		float data_SBOR_time;
		float data_SPE_time;

		//pub
		void protocol_setkey(paillier_pubkey_t* pubkey, paillier_prvkey_t* prvkey, int modulus);
		void protocol_free();

		void protocol_1_free(paillier_ciphertext_t* cipher, bool print=false);
		void protocol_2_free(paillier_ciphertext_t** cipher, int x, bool print= false);
		void protocol_3_free(paillier_ciphertext_t*** cipher, int x, int y, bool print= false);
		void protocol_4_free(paillier_ciphertext_t**** cipher, int x, int y, int z, bool print= false);

		void protocol_1_printf(paillier_ciphertext_t* cipher, char* str, bool print);
		void protocol_2_printf(paillier_ciphertext_t** cipher, int x, char* str, bool print);
		void protocol_3_printf(paillier_ciphertext_t*** cipher, int x, int y, char* str, bool print);
		
		void set_Param(int datanum, int dimension, int bitsize, int requiredK, int tree_level);
		void set_Param(parsed_query* user_query);


		paillier_ciphertext_t* SM_p2(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2);
		paillier_ciphertext_t* SM_p1(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2);
		paillier_ciphertext_t* SM_p1(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2, int idx);

		paillier_ciphertext_t* SSED(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2);
		paillier_ciphertext_t* SSEDm(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2, int col_num);
		
		paillier_ciphertext_t* SBOR(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2);
		
		paillier_ciphertext_t* SBXOR(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2);
		paillier_ciphertext_t* SBXOR(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2, int idx);
				

		paillier_ciphertext_t* SBN(paillier_ciphertext_t* ciper);
		
		paillier_ciphertext_t** SBD(paillier_ciphertext_t* ciper1);
		paillier_ciphertext_t*  SBD_underBob(paillier_ciphertext_t* ciper1, int round);
		paillier_ciphertext_t*  SBD_underAlice(paillier_ciphertext_t* ciper);
		paillier_ciphertext_t** SBD_for_SRO(paillier_ciphertext_t* ciper1, int extra);

		paillier_ciphertext_t*  SSR(paillier_ciphertext_t* ciper);
		paillier_ciphertext_t*  SSR_sub(paillier_ciphertext_t* ciper);
		
		paillier_ciphertext_t** Smin_basic1(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2);
		paillier_ciphertext_t* Smin_basic2(paillier_ciphertext_t** ciper_R, paillier_ciphertext_t** ciper_L, paillier_ciphertext_t** ciper_M, paillier_ciphertext_t* alpha);
		paillier_ciphertext_t* Smin_for_alpha(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2);
		paillier_ciphertext_t** Smin_n(paillier_ciphertext_t*** ciper, int number);
		paillier_ciphertext_t** Smin_bool_n(paillier_ciphertext_t*** ciper, paillier_ciphertext_t** data, int number);
		paillier_ciphertext_t** Smin_bool(paillier_ciphertext_t** data, int number);
		paillier_ciphertext_t** Smin_bool_n_sub(paillier_ciphertext_t** ciper, int number);
		
		paillier_ciphertext_t* SRO(paillier_ciphertext_t*** qLL, paillier_ciphertext_t*** qRR, paillier_ciphertext_t*** nodeLL, paillier_ciphertext_t*** nodeRR);
		paillier_ciphertext_t* SRO2(paillier_ciphertext_t** ciper_L, paillier_ciphertext_t* alpha);


		//garbled
		void makeGate();
		int G_CMP(int mr1, int m1r1, int mr2, int m2r2);


		////////////////////***************range***************////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		int** sRange_B(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int *result_num); //리턴 타입 변경!!
		int** sRange_I(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num);
		int** sRange_G(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num);
		int** sRange_PB(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num);
		int** sRange_PGI(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num);
		int** sRange_PAI(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num);


		//range_i		
		int ** sRange_sub(paillier_ciphertext_t** alpha, int node_num, int * set_num);
		int ** FsRange_Bob(paillier_ciphertext_t*** ciper_result, int rand, int k, int col_num);
		paillier_ciphertext_t*** sNodeRetrievalforRange(paillier_ciphertext_t*** data, paillier_ciphertext_t*** ciper_qLL_bit, paillier_ciphertext_t*** ciper_qRR_bit, boundary* node, int NumData, int NumNode, int* cnt, int* NumNodeGroup, int type);

		//range_g
		paillier_ciphertext_t* GSRO(paillier_ciphertext_t** qLL, paillier_ciphertext_t** qRR, paillier_ciphertext_t** nodeLL, paillier_ciphertext_t** nodeRR);
		paillier_ciphertext_t** GSRO_sub(paillier_ciphertext_t** A, paillier_ciphertext_t** B, int * R1, int * R2);
		paillier_ciphertext_t* DP_GSRO(paillier_ciphertext_t** qLL, paillier_ciphertext_t** qRR, paillier_ciphertext_t** nodeLL, paillier_ciphertext_t** nodeRR);
		paillier_ciphertext_t** unDP_GSRO(paillier_ciphertext_t* A, paillier_ciphertext_t* B, int * R1, int * R2);
		paillier_ciphertext_t*** GSRO_sNodeRetrievalforRange(paillier_ciphertext_t*** data, paillier_ciphertext_t** ciper_qLL, paillier_ciphertext_t** ciper_qRR, boundary* node, int NumData, int NumNode, int* cnt, int* NumNodeGroup, int type);
		paillier_ciphertext_t*** sRange_result(paillier_ciphertext_t** alpha, paillier_ciphertext_t*** cand, int cnt, int NumNodeGroup, int * result_num);

		//parallel
		paillier_ciphertext_t* DP_GSRO_inThread(paillier_ciphertext_t** qLL, paillier_ciphertext_t** qRR,	paillier_ciphertext_t** nodeLL,	paillier_ciphertext_t** nodeRR, int idx);
		paillier_ciphertext_t*** Parallel_GSRO_sNodeRetrievalforRange(paillier_ciphertext_t*** data, paillier_ciphertext_t** cipher_qLL, paillier_ciphertext_t** cipher_qRR, boundary* node, int NumData, int NumNode, int* cnt, int* NumNodeGroup, int type);
		paillier_ciphertext_t*** Parallel_GSRO_inMultithread(int cnt, paillier_ciphertext_t*** cand, boundary q, paillier_ciphertext_t** alpha, paillier_ciphertext_t* cipher_rand);



		//사용x
		paillier_ciphertext_t* faster_SRO(paillier_ciphertext_t*** qLL, paillier_ciphertext_t*** qRR, paillier_ciphertext_t*** nodeLL, paillier_ciphertext_t*** nodeRR); 
		paillier_plaintext_t*  faster_SRO2(paillier_ciphertext_t** ciper_L, paillier_plaintext_t* alpha);
		int** faster_sRange(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num);
		paillier_ciphertext_t* SPE(paillier_ciphertext_t*** qBound, paillier_ciphertext_t*** nodeLL, paillier_ciphertext_t*** nodeRR);
	

		////////////////////knn///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		int * SkNN_B(paillier_ciphertext_t** ciper, int Q, int k, int row_number);
		int ** SkNN_B(paillier_ciphertext_t*** ciper, paillier_ciphertext_t** qeury, int k, int row_number);
		int** SkNN_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, boundary* node, int k, int NumData, int NumNode); // Proposed skNN with secure Index (ICDE 2014 + secure index)
		int** SkNN_G(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode);
		int** SkNN_PB(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode);
		int** SkNN_PGI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode);
		int** SkNN_PAI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode);

		int * SkNNm_Bob(paillier_ciphertext_t** ciper_result_array, int rand, int k);	
		int ** SkNNm_Bob2(paillier_ciphertext_t*** ciper_result, int rand, int k, int col_num);
		paillier_ciphertext_t*** sNodeRetrievalforkNN(paillier_ciphertext_t*** data, paillier_ciphertext_t*** ciper_qLL_bit, paillier_ciphertext_t*** ciper_qRR_bit, boundary* node,paillier_ciphertext_t** alpha,  int NumData, int NumNode, int* cnt, int* NumNodeGroup);
		paillier_ciphertext_t* SCMP(paillier_ciphertext_t** u, paillier_ciphertext_t** v);
		paillier_ciphertext_t* SCMP_M(paillier_ciphertext_t** u, paillier_ciphertext_t** v);
		paillier_ciphertext_t* SCMP_for_SBD(paillier_ciphertext_t** u, paillier_ciphertext_t** v);
		paillier_ciphertext_t* SCMP_Clustering(paillier_ciphertext_t** u, paillier_ciphertext_t** v);

		paillier_ciphertext_t* SC(paillier_ciphertext_t* u, paillier_ciphertext_t* v);
		paillier_plaintext_t* SETC(paillier_ciphertext_t*** former, paillier_ciphertext_t** formerCnt, paillier_ciphertext_t*** New, paillier_ciphertext_t** NewCnt);


		//SSED
		paillier_ciphertext_t* DP_SSED(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2, int col_num);
		paillier_ciphertext_t* unDP_SSED(paillier_ciphertext_t* ciper1, int col_num);
		int ** SSED_test(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, boundary* node, int k, int NumData, int NumNode);
		paillier_ciphertext_t* OP_SSED(paillier_ciphertext_t** data, paillier_ciphertext_t** Cluster_Center, paillier_ciphertext_t** Center_cnt, int tmp_k, int k);

		//knn_m

		paillier_ciphertext_t** SkNNm_sub(paillier_ciphertext_t** ciper_n, int n);		


		//knn_i

		//knn_g
		paillier_ciphertext_t*** GSRO_sNodeRetrievalforkNN(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha,  int NumData, int NumNode, int* cnt, int* NumNodeGroup);
		int SMSn(paillier_ciphertext_t** ciper, int cnt);
		paillier_ciphertext_t* SMSn_sub(paillier_ciphertext_t* ciper, paillier_ciphertext_t* rand ,int cnt);		
		
		//2015.10.06 작업 시작 ~ 
		paillier_ciphertext_t* GSCMP(paillier_ciphertext_t* u, paillier_ciphertext_t* v);
		paillier_ciphertext_t* GSCMP_sub(paillier_ciphertext_t* ru, paillier_ciphertext_t* rv, paillier_ciphertext_t* ciper_Rand);

		//knn_B
		int ** SkNN_b(paillier_ciphertext_t*** ciper, paillier_ciphertext_t** query, int k, int row_number);
		int * SkNNb_C2(paillier_ciphertext_t** ciper_dist, int k, int row_number);
		int ** SkNNb_Bob(paillier_ciphertext_t*** ciper_result, int rand, int k);
		int ** GSRO_SkNNb(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode);
		//knn_DP
		int** DP_SkNN_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, boundary* node, int k, int NumData, int NumNode);

		////////////////////top-k//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//top-k
		int ** STopk_B(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, int NumData);
		int** STopk_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode);
		int** STopk_G(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode);
		int** STopk_PB(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode);
		int** STopk_PGI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode);
		int** STopk_PAI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode);



		// topk parallel
		void Parallel_GSRO_Topk(paillier_ciphertext_t** cipher_qLL, paillier_ciphertext_t** cipher_qRR, paillier_ciphertext_t** alpha, boundary* node, int NumNode, int* cnt, int* NumNodeGroup, bool type);
		paillier_ciphertext_t*** Parallel_GSRO_sNodeRetrievalforTopk(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup);
		void ComputeScoreinMultithread(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, paillier_ciphertext_t** SCORE, int * cnt, bool type);
		void MAXnMultithread2(int cnt, paillier_ciphertext_t** DIST_MINUS_MIN, paillier_ciphertext_t** DIST, paillier_ciphertext_t* MIN, paillier_ciphertext_t* C_RAND);
		void MAXnMultithread3(int cnt, int s, paillier_ciphertext_t **V, paillier_ciphertext_t ***V2, paillier_ciphertext_t **SCORE, paillier_ciphertext_t ***cand, paillier_ciphertext_t *MIN, paillier_ciphertext_t ***Result);


		int** STopk(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode);
		paillier_ciphertext_t* computeScore(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2);
		paillier_ciphertext_t* computeScore2(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2, paillier_ciphertext_t** coeff, paillier_ciphertext_t* hint);
		paillier_ciphertext_t** Smax_n(paillier_ciphertext_t*** ciper, int number);
		paillier_ciphertext_t** Smax_basic1(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2);
		paillier_ciphertext_t* Smax_basic2(paillier_ciphertext_t** ciper_R, paillier_ciphertext_t** ciper_L, paillier_ciphertext_t** ciper_M, paillier_ciphertext_t* alpha);
		paillier_ciphertext_t** Topk_sub(paillier_ciphertext_t** ciper_n, int n);
		int** STopk_Confirm(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode);
		paillier_ciphertext_t*** sNodeRetrievalforTopk(paillier_ciphertext_t*** data, paillier_ciphertext_t*** ciper_qLL_bit, paillier_ciphertext_t*** ciper_qRR_bit, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup);
		paillier_ciphertext_t*** GSRO_sNodeRetrievalforTopk(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha,  int NumData, int NumNode, int* cnt, int* NumNodeGroup);

		//top-k_g
		int MAXn(paillier_ciphertext_t** ciper, int cnt);

		//Classification
		int** Classification_M(paillier_ciphertext_t*** ciper, paillier_ciphertext_t** qeury, paillier_ciphertext_t** Entire_set, int k, int row_number, int Entire_num);
		int** Classification_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, paillier_ciphertext_t** Entire_set, boundary* node, int k, int NumData, int NumNode, int Entire_num);
		int** Classification_G(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, paillier_ciphertext_t** Entire_set, boundary* node, int k, int NumData, int NumNode, int Entire_num);


		paillier_ciphertext_t** SF_P1(paillier_ciphertext_t** Entire_set, paillier_ciphertext_t** K_set, int w, int k);
		paillier_ciphertext_t*** SF_P2(paillier_ciphertext_t*** S, int w, int k);
		paillier_ciphertext_t* SCMC(paillier_ciphertext_t** Entire_set, paillier_ciphertext_t** K_set, int w, int k);
		paillier_ciphertext_t*** sNodeRetrievalforClassification(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup);


		//Classification parallel
		void SSEDMultiThread(paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t **q, int cnt);
		void SMSnMultithread2(int cnt, paillier_ciphertext_t** DIST_MINUS_MIN, paillier_ciphertext_t** DIST, paillier_ciphertext_t* MIN, paillier_ciphertext_t* C_RAND);
		void SMSnMultithread3(int cnt, int s, paillier_ciphertext_t **V, paillier_ciphertext_t ***V2, paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, paillier_ciphertext_t ***Result);
		paillier_ciphertext_t*** sNodeRetrievalforClassificationMultiThread(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha,	int NumData, int NumNode, int* cnt, int* NumNodeGroup);

		//Clustering 2017 02 01
		paillier_ciphertext_t*** Clustering_m(paillier_ciphertext_t*** ciper, int NumData, int k, int b);
		paillier_ciphertext_t*** Clustering_Grid(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k);
		paillier_ciphertext_t*** PreClustering(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k);
		paillier_ciphertext_t*** extract_sample(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int* Sample_Cnt);
		paillier_ciphertext_t*** Clustering_m_Grid(paillier_ciphertext_t*** ciper, int NumData, int k, int b, paillier_ciphertext_t*** former_Center);
		paillier_ciphertext_t* boundary_dist(boundary node, paillier_ciphertext_t** former_Center, int NumNode);
		paillier_ciphertext_t*** Clustering_Grid_preprocessing(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k);
		paillier_ciphertext_t** G_Smin_bool(paillier_ciphertext_t** data, int number);
};



#endif