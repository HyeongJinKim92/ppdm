#ifndef _PARALLEL_PROTOCOL_H_
#define	 _PARALLEL_PROTOCOL_H_

#include <iostream>
#include <fstream>
#include <gmp.h>
#include <string>
#include <sstream>
#include <thread>
#include <chrono>
#include <vector>
#include <iostream>
#include "paillier.h"
#include "protocol.h"

//RANGE_PB
void SRO_SBD_inThread(std::vector<int> input, paillier_ciphertext_t *** ciper_qLL_bit, paillier_ciphertext_t *** ciper_qRR_bit, paillier_ciphertext_t *** data, paillier_ciphertext_t ** alpha, protocol proto);


//RANGE_PGI
void DP_GSRO_Multithread(paillier_ciphertext_t **qLL, paillier_ciphertext_t **qRR, std::vector<int> input, boundary *node, paillier_ciphertext_t **alpha, protocol proto, int idx);
void sNodeRetrievalforDatainThread(std::vector<int> inputJ, std::vector<int *> inputNodeGroup, int &remained,	std::vector<int> inputCnt, boundary *node, paillier_ciphertext_t ***data, paillier_ciphertext_t **alpha, paillier_ciphertext_t ***cand, protocol proto, int idx);
void DP_GSRO_dataquery_Multithread(paillier_ciphertext_t *** cand, boundary q, std::vector<int> input, paillier_ciphertext_t **alpha, protocol proto, int idx, paillier_ciphertext_t* cipher_rand);

//TOPK PB
void CS_SBD_inThread(std::vector<int>input, paillier_ciphertext_t ** query, paillier_ciphertext_t *** data, paillier_ciphertext_t ** cipher_distance, paillier_ciphertext_t *** cipher_SBD_distance, protocol proto);
void Smax_n_InThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t ***copy_cipher, paillier_ciphertext_t ***min, protocol proto);
void recalculate_dist_inThread(std::vector<int>input, paillier_ciphertext_t*** cipher_SBD_dist, paillier_ciphertext_t** cipher_distance, protocol proto);
void DuringProcessedInTOPK_PB_inThread( std::vector<int>input, paillier_ciphertext_t ** cipher_distance, paillier_ciphertext_t ** cipher_mid, paillier_ciphertext_t * cipher_min, paillier_ciphertext_t *cipher_rand, protocol proto);
void ExtractTOPKInTOPK_PB_inThread( std::vector<int>input, paillier_ciphertext_t*** data, paillier_ciphertext_t** cipher_V, paillier_ciphertext_t ** cipher_thread_result, protocol proto);
void UPDATE_SBD_SCORE_InTOPK_PB_inThread(std::vector<int>input, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_V, protocol proto);

//TOPK PGI
void DP_GSRO_MultithreadforTopk(paillier_ciphertext_t **qLL, paillier_ciphertext_t **qRR, std::vector<int> input, boundary *node, paillier_ciphertext_t **alpha, protocol proto, bool verify);
void ComputeScore_Multithread(paillier_ciphertext_t *** cand, paillier_ciphertext_t ** q, paillier_ciphertext_t ** SCORE, std::vector<int> input, protocol proto, bool idx);
void MAXnInthread2(std::vector<int> inputCnt, paillier_ciphertext_t **DIST_MINUS_MIN, paillier_ciphertext_t **DIST,	paillier_ciphertext_t *MIN, paillier_ciphertext_t *C_RAND, protocol proto, int idx);
void MAXnInThread3_1(int index, int s, std::vector<int> inputCnt, paillier_ciphertext_t **V, paillier_ciphertext_t **DIST,	paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, protocol proto);
void MAXnInThread3_2_v2(int index, int offset, int record_size, int dim, paillier_ciphertext_t** sub_result, paillier_ciphertext_t **V, paillier_ciphertext_t ***V2, paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, protocol proto);
void GSCMPforTopk_inThread(std::vector<int> input, boundary* node, paillier_ciphertext_t** psi, paillier_ciphertext_t* kth_score, paillier_ciphertext_t** q, paillier_ciphertext_t** alpha, protocol proto);

//kNN GB
void recalculate_DISTforkNN_inThread(std::vector<int> input, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_distance, protocol proto);
void UPDATE_SBD_SCORE_InKNN_PB_inThread(std::vector<int>input, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_V, protocol proto);


//CLASSIFICATION_PGI
void SSEDInThread(std::vector<int> inputCnt, paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t **q, int cnt, protocol proto);
void SMSnInthread2(std::vector<int> inputCnt, paillier_ciphertext_t **DIST_MINUS_MIN, paillier_ciphertext_t **DIST,	paillier_ciphertext_t *MIN, paillier_ciphertext_t *C_RAND, protocol proto, int idx);
void SMSnInThread3_1(int index, int s, std::vector<int> inputCnt, paillier_ciphertext_t **V, paillier_ciphertext_t **DIST,	paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, protocol proto);
void SMSnInThread3_2_v2(int index, int offset, int record_size, int dim, paillier_ciphertext_t** sub_result, paillier_ciphertext_t **V, paillier_ciphertext_t ***V2, paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, protocol proto);
void sNodeRetrievalforClassificationinThread(std::vector<int> inputJ, std::vector<int *> inputNodeGroup, int &remained,	std::vector<int> inputCnt, boundary *node, paillier_ciphertext_t ***data, paillier_ciphertext_t **alpha,	paillier_ciphertext_t ***cand, protocol proto, int idx);
void NodeExpansion_Multithread(paillier_ciphertext_t **q, paillier_ciphertext_t *K_DIST, std::vector<int> input, boundary *node, paillier_ciphertext_t **alpha, protocol proto, int idx);



//CLASSIFICATION_PB
void SSED_SBD_inThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t **query, paillier_ciphertext_t ***data, int dim, paillier_ciphertext_t **cipher_distance, paillier_ciphertext_t ***cipher_SBD_distance, protocol proto);
void SSED_inThread(std::vector<int> inputIdx, paillier_ciphertext_t **query, paillier_ciphertext_t ***data, paillier_ciphertext_t **cipher_distance, protocol proto);
void SBD_inThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t **query, paillier_ciphertext_t ***data, paillier_ciphertext_t **cipher_distance, paillier_ciphertext_t ***cipher_SBD_distance, protocol proto);
void Smin_n_InThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t ***copy_cipher, paillier_ciphertext_t ***min, protocol proto);
void caldistanceInThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t ***cipher_SBD_distance, paillier_ciphertext_t **cipher_distance, protocol proto);
void classipmInThread1(int index, std::vector<int> inputIdx, paillier_ciphertext_t **cipher_distance, paillier_ciphertext_t *cipher_min, paillier_ciphertext_t *cipher_rand, paillier_ciphertext_t **cipher_mid, protocol proto);
void classipmInThread2(int index, std::vector<int> inputIdx, paillier_ciphertext_t **cipher_V, paillier_ciphertext_t **cipher_V2, paillier_ciphertext_t ***data, paillier_ciphertext_t **cipher_rabel_result, int s, protocol proto);
void SBOR_inThread(int index, std::vector<int> inputIdx, paillier_ciphertext_t **cipher_V, paillier_ciphertext_t ***cipher_SBD_distance, protocol proto);


//CLUSTERING
void Comp_Cluster_inThread(paillier_ciphertext_t*** cipher, paillier_ciphertext_t*** former_Center, std::vector<int> inputIdx, paillier_ciphertext_t*** NewSumCluster, paillier_ciphertext_t** NewSumCntCluster, protocol proto, int i);


#endif