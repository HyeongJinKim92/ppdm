#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include "protocol.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "util/config.h"
#include "circuit/circuit.h"
/*
#define bool int
#define TRUE 1
#define FALSE 0
*/


CConfig* pConfig = new CConfig();
CCircuit* pCircuit = NULL;

using namespace std;

std::mutex protocol::mtx;

void protocol::set_Param(parsed_query * user_query)
{
	query = user_query->query;
	if(query == RANGE)					rquery_t = user_query->rquery_t;
	else if(query == TOPK)				tquery_t = user_query->tquery_t;
	else if(query == KNN)				kquery_t = user_query->kquery_t;
	else if(query == CLASSIFICATION)	cquery_t = user_query->cquery_t;
	else if(query == KMEANS)			mquery_t = user_query->mquery_t;
	else if(query == ASSOCATION)		{}
	else if(query == TEST)			{}
	
	NumData = user_query->data_n; // data 
	thread_num = user_query->thread_n; // thread
	dim = user_query->dim_n; // dim
	size = user_query->bit_s; // bitsize
	k = user_query->k; // k
	FanOut = (user_query->data_n / pow(2, user_query->tree_level-1))+1; // FanOut
	printf("tree_level : %d\n", user_query->tree_level);
	printf("datanum : %d, dim : %d, bitsize : %d, requiredK : %d, FanOut : %d\n", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->k, FanOut);
	totalNumOfRetrievedNodes = 0;
	if (size >= 15)
	{
		DP = 32;
	}
	else
	{
		DP = 16;
	}

	//for smsn
	smsn_cnt =0;
	smsn_flag= false;
	smsn_tmp_val = 0;
	smsn_tmp_rand = 0;
	smsn_result_idx = 0;
	Print = FALSE;
	total_time = 0.0;
	node_SBD_time = 0.0;
	node_SRO_time = 0.0;
	node_expansion_time = 0.0;
	data_extract_first_time = 0.0;
	data_extract_second_time = 0.0;
	data_SSED_SBD_time = 0.0;
	sMINn_first_time = 0.0;
	sMINn_second_time = 0.0;
	data_SBOR_time = 0.0;
	data_SPE_time = 0.0;


	makeGate();
}


void	protocol::set_Param(int datanum, int dimension, int bitsize, int requiredK, int tree_level)
{
	NumData = datanum;
	dim = dimension;
	size = bitsize;
	k = requiredK;
	FanOut = (datanum / pow(2, tree_level-1))+1;
	printf("tree_level : %d\n", tree_level);
	printf("datanum : %d, dim : %d, bitsize : %d, requiredK : %d, FanOut : %d\n", datanum, dimension, bitsize, requiredK, FanOut);
	totalNumOfRetrievedNodes = 0;
	if (size >= 15)
	{
		DP = 32;
	}
	else
	{
		DP = 16;
	}

	//for smsn
	smsn_cnt =0;
	smsn_flag= false;
	smsn_tmp_val = 0;
	smsn_tmp_rand = 0;
	smsn_result_idx = 0;
	Print = FALSE;
	total_time = 0.0;
	node_SBD_time = 0.0;
	node_SRO_time = 0.0;
	node_expansion_time = 0.0;
	data_extract_first_time = 0.0;
	data_extract_second_time = 0.0;
	data_SSED_SBD_time = 0.0;
	sMINn_first_time = 0.0;
	sMINn_second_time = 0.0;
	data_SBOR_time = 0.0;
	data_SPE_time = 0.0;
}

void protocol::protocol_free(){
	free(pubkey);
	free(prvkey);
	
	paillier_freeplaintext(plain_zero);		paillier_freeplaintext(plain_one);		paillier_freeplaintext(plain_two);		
	paillier_freeplaintext(plain_minus);		paillier_freeplaintext(plain_minus_two);
	paillier_freeplaintext(plain_MAX);		
	
	paillier_freeciphertext(ciper_zero); 	paillier_freeciphertext(ciper_one); paillier_freeciphertext(ciper_minus);	paillier_freeciphertext(l);	// 얘가 뭐하는 애지!!??
	
	free(temp);	// 2차원 포인터 해제로 변경해야 함!!
	paillier_freeciphertext(result);		
	paillier_freeciphertext(ciper_MAX);
}

void protocol::protocol_setkey(paillier_pubkey_t* pubkey, paillier_prvkey_t* prvkey, int modulus)
{
	modul = modulus;
	int i;
	this->pubkey = pubkey;
	this->prvkey = prvkey;

	plain_zero = paillier_plaintext_from_ui(0);	// plain_zero를 0으로 셋팅
	plain_one = paillier_plaintext_from_ui(1);		// plain_one을 1로 셋팅
	plain_two = paillier_plaintext_from_ui(2);		// plain_two을 2로 셋팅

	// plain_minus을 -1로 셋팅
	plain_minus = paillier_plaintext_from_ui(0);
	mpz_sub_ui(plain_minus->m, plain_minus->m, 1);

	plain_minus_two = paillier_plaintext_from_ui(2);	
	mpz_powm(plain_minus_two->m, plain_minus_two->m, plain_minus->m, pubkey->n_squared);
	l = paillier_enc(0, pubkey, plain_minus_two, paillier_get_rand_devurandom);

	plain_MAX = paillier_plaintext_from_ui(pow(2, size/2));


	ciper_zero = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);	// ciper_zero를 암호화 0으로 셋팅
	ciper_one = paillier_enc(0, pubkey, plain_one, paillier_get_rand_devurandom);	// ciper_one을 암호화 1로 셋팅
	ciper_minus = paillier_enc(0, pubkey, plain_minus, paillier_get_rand_devurandom);		// ciper_minus를 암호화 -1로 셋팅
	ciper_MAX = paillier_enc(0, pubkey, plain_MAX, paillier_get_rand_devurandom);


	rand1 = paillier_plaintext_from_ui(5);
	rand2 = paillier_plaintext_from_ui(4);

	ciper_rand1 = paillier_enc(0, pubkey, rand1, paillier_get_rand_devurandom);
	ciper_rand2 = paillier_enc(0, pubkey, rand2, paillier_get_rand_devurandom);



	temp  = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*3);
	for( i = 0 ; i < 3 ; i++ ){
		temp[i] = paillier_create_enc_zero();
	}

	result = paillier_create_enc_zero();
}

paillier_ciphertext_t* protocol::SBN(paillier_ciphertext_t* ciper) {
	paillier_ciphertext_t* result = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(result->c);

	paillier_exp(pubkey, result, ciper, plain_minus);	// multiply -1
	paillier_mul(pubkey, result, result, ciper_one);	// plus 1

	return result;
}

// KHJ add
paillier_ciphertext_t** protocol::SBD(paillier_ciphertext_t* ciper1){
	int i = 0;

	paillier_ciphertext_t* t = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(t->c);

	paillier_ciphertext_t** ciper_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_array_reverse = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);

	for( i = 0 ; i < size ; i++){
		ciper_array[i] = SBD_underBob(ciper1, i);
		paillier_subtract(pubkey, t, ciper1, ciper_array[i]);
		ciper1 = SM_p1(t, l);
		ciper_array_reverse[size - 1 - i] = ciper_array[i];
	}

	paillier_freeciphertext(t);

	return ciper_array_reverse;
}

// added by KHI. 150712
// int extra : 비트 변환 후 가장 마지막 행에 u=0, v=1인 행을 추가해주기 위한 정보
paillier_ciphertext_t** protocol::SBD_for_SRO(paillier_ciphertext_t* ciper1, int extra){
	int i = 0;

	paillier_ciphertext_t* t = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(t->c);


	paillier_ciphertext_t** ciper_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_array_reverse = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	for( i = 0 ; i < size ; i++){
		ciper_array[i] = SBD_underBob(ciper1, i);
		paillier_subtract(pubkey, t, ciper1, ciper_array[i]);
		ciper1 = SM_p1(t, l);
		ciper_array_reverse[size - 1 - i] = ciper_array[i];
	}

	if(extra == 1)
		ciper_array_reverse[size] = ciper_one;
	else
		ciper_array_reverse[size] = ciper_zero;

	paillier_freeciphertext(t);

	return ciper_array_reverse;
}

paillier_ciphertext_t*  protocol::SBD_underBob(paillier_ciphertext_t* ciper1, int round){
	int r = rand()%2;

	paillier_plaintext_t* rand = paillier_plaintext_from_ui(r);
	paillier_ciphertext_t* ciper_rand = paillier_enc(0, pubkey, rand , paillier_get_rand_devurandom);	
	paillier_mul(pubkey, ciper_rand, ciper1, ciper_rand);
	paillier_ciphertext_t* sub_result = SBD_underAlice(ciper_rand);

	if( r % 2 == 0 ){
	
	}else{
		paillier_plaintext_t* one = paillier_plaintext_from_ui(1);	
		paillier_ciphertext_t* ciper_one = paillier_enc(0, pubkey, one, paillier_get_rand_devurandom);
		paillier_subtract(pubkey, sub_result, ciper_one, sub_result); 
	}
	paillier_freeciphertext(ciper_rand);
	paillier_freeplaintext(rand);

	return sub_result;
}

paillier_ciphertext_t*  protocol::SBD_underAlice(paillier_ciphertext_t* ciper){
	paillier_plaintext_t* plain = paillier_dec(0, pubkey, prvkey, ciper);

	mpz_mod(plain->m, plain->m, plain_two->m);

	if( mpz_cmp(plain_zero->m, plain->m) == 0){
		paillier_ciphertext_t* ca = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);

		return ca;
	}else{
		paillier_ciphertext_t* ca = paillier_enc(0, pubkey, plain_one, paillier_get_rand_devurandom);

		return ca;
	}
}

paillier_ciphertext_t* protocol::SM_p1(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2){
	int i = 0;
 /*
	// generate noise (random numbers)
	paillier_plaintext_t* rand1=paillier_plaintext_from_ui(5);
	paillier_plaintext_t* rand2=paillier_plaintext_from_ui(4);

	// encrypt noise 
	paillier_ciphertext_t* ciper_rand1 = paillier_enc(0, pub, rand1, paillier_get_rand_devurandom);
	paillier_ciphertext_t* ciper_rand2 = paillier_enc(0, pub, rand2, paillier_get_rand_devurandom);
*/
	
	paillier_ciphertext_t** temp  = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*3);
	for( i = 0 ; i < 3 ; i++ ){
		temp[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(temp[i]->c);
	}

	// compute (c1 + rd2) and (c2 + rd1)
	paillier_mul(pub, temp[0], ciper1, ciper_rand1);
	paillier_mul(pub, temp[1], ciper2, ciper_rand2);

	// compute (c1 + rd2) x (c2 + rd1)
	result = SM_p2(temp[0], temp[1]);

	// delete c1 x rnd2 
	paillier_exp(pub, temp[2], ciper1, rand2);
	paillier_subtract(pub, result, result, temp[2]);

	// delete c2 x rnd1
	paillier_exp(pub, temp[2], ciper2, rand1);
	paillier_subtract(pub, result, result, temp[2]);

	// delete rnd1 x rnd2
	paillier_exp(pub, temp[2], ciper_rand1, rand2);
	paillier_subtract(pub, result, result, temp[2]);
/*
	paillier_freeplaintext(rand1);		paillier_freeplaintext(rand2);
	paillier_freeciphertext(ciper_rand1);	paillier_freeciphertext(ciper_rand2);
*/
	return result;
}
paillier_ciphertext_t* protocol::SM_p1(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2, int idx){
	int i = 0;
 	
	paillier_ciphertext_t** temp  = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*3);
	for( i = 0 ; i < 3 ; i++ ){
		temp[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(temp[i]->c);
	}

	// compute (c1 + rd2) and (c2 + rd1)
	paillier_mul(pub, temp[0], ciper1, ciper_rand1);
	paillier_mul(pub, temp[1], ciper2, ciper_rand2);

	// compute (c1 + rd2) x (c2 + rd1)
	result = SM_p2(temp[0], temp[1]);

	// delete c1 x rnd2 
	paillier_exp(pub, temp[2], ciper1, rand2);
	paillier_subtract(pub, result, result, temp[2]);

	// delete c2 x rnd1
	paillier_exp(pub, temp[2], ciper2, rand1);
	paillier_subtract(pub, result, result, temp[2]);

	// delete rnd1 x rnd2
	paillier_exp(pub, temp[2], ciper_rand1, rand2);
	paillier_subtract(pub, result, result, temp[2]);
/*
	paillier_freeplaintext(rand1);		paillier_freeplaintext(rand2);
	paillier_freeciphertext(ciper_rand1);	paillier_freeciphertext(ciper_rand2);
*/
	return result;
}

paillier_ciphertext_t* protocol::SM_p2(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2){
	paillier_plaintext_t* plain1 = paillier_dec(0, pubkey, prvkey, ciper1);
	paillier_plaintext_t* plain2 = paillier_dec(0, pubkey, prvkey, ciper2);
	paillier_plaintext_t* result = (paillier_plaintext_t*) malloc(sizeof(paillier_plaintext_t));
	mpz_init(result->m);

	plain_mul(pubkey, result, plain1, plain2);	

	//gmp_printf("The SM_p2 is : %Zd\n",temp);   //test
	//paillier_freeplaintext(plain1);
	//paillier_freeplaintext(plain2);
	
	return paillier_enc(0, pubkey, result, paillier_get_rand_devurandom);
}


paillier_ciphertext_t* protocol::SSED(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2)
{
	// init variables
	paillier_ciphertext_t* dist = paillier_create_enc_zero();
	paillier_ciphertext_t* tmp_dist = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(tmp_dist->c);
	
	paillier_subtract(pubkey, tmp_dist, ciper1, ciper2);
	dist = SM_p1(tmp_dist, tmp_dist);

	paillier_print("dist : ", dist);
	paillier_freeciphertext(tmp_dist);
	return dist;
}

paillier_ciphertext_t* protocol::SSEDm(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2, int col_num)
{
	// init variables
	int i = 0;
	paillier_ciphertext_t* dist = paillier_create_enc_zero();
	paillier_ciphertext_t* tmp_dist = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(tmp_dist->c);

	for( i = 0 ; i < col_num ; i++ )
	{
		paillier_subtract(pubkey, tmp_dist, ciper1[i], ciper2[i]);
		//paillier_print("paillier_subtract tmp_dist : ", tmp_dist);
		tmp_dist = SM_p1(tmp_dist, tmp_dist);	
		//paillier_print("tmp_dist : ", tmp_dist);
		paillier_mul(pubkey, dist, dist, tmp_dist);
		//paillier_print("dist : ", dist);		
	}

	//printf("SSEDm\n");

	//paillier_print("dist : ", dist);

	paillier_freeciphertext(tmp_dist);

	return dist;
}


paillier_ciphertext_t* protocol::SBOR(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2)
{
	paillier_ciphertext_t* result = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(result->c);
	paillier_ciphertext_t* temp1 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp1->c);
	paillier_ciphertext_t* temp2 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp2->c);

	// compute (c1 + c2)
	paillier_mul(pubkey, temp1, ciper1, ciper2);

	// compute  (c1*c2)
	temp2 = SM_p1(ciper1, ciper2);

	// compute (c1 + c2 - c1*c2)
	paillier_subtract(pubkey, result, temp1, temp2);

	//paillier_print("result : ", result);

	paillier_freeciphertext(temp1);		paillier_freeciphertext(temp2);

	return result;
}

paillier_ciphertext_t* protocol::SBXOR(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2)
{
	paillier_ciphertext_t* result = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(result->c);
	paillier_ciphertext_t* temp1 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp1->c);
	paillier_ciphertext_t* temp2 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp2->c);

	// compute (c1 + c2)
	paillier_mul(pubkey, temp1, ciper1, ciper2);

	// compute  (2*c1*c2)
	temp2 = SM_p1(ciper1, ciper2);
	paillier_exp(pubkey, temp2, temp2, plain_two);

	// compute (c1 + c2 -2*c1*c2)
	paillier_subtract(pubkey, result, temp1, temp2);

	//paillier_print("result : ", result);

	paillier_freeciphertext(temp1);		paillier_freeciphertext(temp2);

	return result;
}
paillier_ciphertext_t* protocol::SBXOR(paillier_ciphertext_t* ciper1, paillier_ciphertext_t* ciper2, int idx)
{
	paillier_ciphertext_t* result = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(result->c);
	paillier_ciphertext_t* temp1 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp1->c);
	paillier_ciphertext_t* temp2 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp2->c);

	// compute (c1 + c2)
	paillier_mul(pubkey, temp1, ciper1, ciper2);

	// compute  (2*c1*c2)
	temp2 = SM_p1(ciper1, ciper2, idx);
	paillier_exp(pubkey, temp2, temp2, plain_two);

	// compute (c1 + c2 -2*c1*c2)
	paillier_subtract(pubkey, result, temp1, temp2);

	//paillier_print("result : ", result);

	paillier_freeciphertext(temp1);		paillier_freeciphertext(temp2);

	return result;
}


paillier_ciphertext_t** protocol::Smin_basic1(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2)
{
	int i;
	//BOOL func = FALSE;
	bool func = false;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* ciper_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);

	paillier_ciphertext_t** ciper_W = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_R	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_H	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_L	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_M	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_lambda = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_min	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_O = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_G	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** temp3	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);

	for(i=0; i<size; i++){
		ciper_W[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_R[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_H[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_L[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_M[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_lambda[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_min[i]	= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_O[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_G[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp[i]			= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp2[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp3[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));

		mpz_init(ciper_W[i]->c);
		mpz_init(ciper_R[i]->c);
		mpz_init(ciper_H[i]->c);
		mpz_init(ciper_L[i]->c);
		mpz_init(ciper_M[i]->c);
		mpz_init(ciper_lambda[i]->c);
		mpz_init(ciper_min[i]->c);
		mpz_init(ciper_O[i]->c);
		mpz_init(ciper_G[i]->c);
		mpz_init(temp[i]->c);
		mpz_init(temp2[i]->c);
		mpz_init(temp3[i]->c);
	}

	for(i=0; i<size; i++){
		if(func){	// true :  F : u>v	
			paillier_subtract(pubkey, ciper_W[i], ciper1[i], SM_p1(ciper1[i], ciper2[i]));	// W
			paillier_subtract(pubkey,temp[i],ciper2[i],ciper1[i]);	
			paillier_mul(pubkey,ciper_R[i],temp[i],ciper_Rand_value);	// Gamma
		}else{
			paillier_subtract(pubkey, ciper_W[i], ciper2[i], SM_p1(ciper1[i], ciper2[i]));
			paillier_subtract(pubkey, temp[i], ciper1[i], ciper2[i]);
			paillier_mul(pubkey, ciper_R[i], temp[i], ciper_Rand_value);
		}
		ciper_G[i]=SBXOR(ciper1[i],ciper2[i]);

		if(i==0){
			paillier_exp(pubkey,ciper_H[i], ciper_zero, Rand_value);
			paillier_mul(pubkey,ciper_H[i], ciper_H[i], ciper_G[i]);		
		}else{
			paillier_exp(pubkey,temp[i],ciper_H[i-1],Rand_value);
			paillier_mul(pubkey,ciper_H[i],temp[i],ciper_G[i]);
		}

		paillier_mul(pubkey,ciper_O[i],ciper_H[i],ciper_minus);	// PI
		paillier_exp(pubkey,ciper_O[i],ciper_O[i],Rand_value);
		// paillier_exp(pubkey,temp3[i],ciper_W[i],Rand_value);   // our IDEA
		// paillier_mul(pubkey,ciper_L[i],temp2[i], temp3[i]);	 // our IDEA
		
		paillier_mul(pubkey, ciper_L[i], ciper_O[i], ciper_W[i]);
	}

	alpha = Smin_basic2(ciper_R, ciper_L, ciper_M, alpha);

	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));
	
	for(i=0;i<size;i++){
		paillier_exp(pubkey,tmp,alpha,Rand_value);
		paillier_subtract(pubkey, ciper_lambda[i], ciper_M[i], tmp);
		if(func){
			paillier_mul(pubkey,ciper_min[i],ciper1[i],ciper_lambda[i]);			
		}else{
			paillier_mul(pubkey,ciper_min[i],ciper2[i],ciper_lambda[i]);			
		}
	}

	paillier_freeplaintext(Rand_value);		 
	paillier_freeciphertext(ciper_Rand_value);	paillier_freeciphertext(tmp);
	
	paillier_freeciphertext(alpha);

	for(i=0; i<size; i++){
		paillier_freeciphertext(ciper_W[i]);	paillier_freeciphertext(ciper_R[i]); paillier_freeciphertext(ciper_H[i]);	paillier_freeciphertext(ciper_L[i]);	paillier_freeciphertext(ciper_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	paillier_freeciphertext(temp3[i]);
		paillier_freeciphertext(ciper_lambda[i]);
	}
	
	return ciper_min;
}

paillier_ciphertext_t* protocol::Smin_for_alpha(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2)
{
	int i;
	bool func = false;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* ciper_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);

	paillier_ciphertext_t** ciper_W = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_G	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_H	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_L	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_M	 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_PI = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);

	for(i=0; i<size; i++){
		ciper_W[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_G[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_H[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_L[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_M[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_PI[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));

		temp[i]			= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp2[i]		= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));

		mpz_init(ciper_W[i]->c);
		mpz_init(ciper_G[i]->c);
		mpz_init(ciper_H[i]->c);
		mpz_init(ciper_L[i]->c);
		mpz_init(ciper_M[i]->c);
		mpz_init(ciper_PI[i]->c);

		mpz_init(temp[i]->c);
		mpz_init(temp2[i]->c);
	}

	for(i=0; i<size; i++){
		if(func){	// true :  F : u>v	
			paillier_subtract(pubkey, ciper_W[i], ciper1[i], SM_p1(ciper1[i], ciper2[i]));	// W
		}else{
			paillier_subtract(pubkey, ciper_W[i], ciper2[i], SM_p1(ciper1[i], ciper2[i]));
		}

		// compute G
		ciper_G[i]=SBXOR(ciper1[i],ciper2[i]);

		// compute H
		if(i==0){
			paillier_exp(pubkey,ciper_H[i], ciper_zero, Rand_value);
			paillier_mul(pubkey,ciper_H[i], ciper_H[i], ciper_G[i]);		
		}else{
			paillier_exp(pubkey,temp[i],ciper_H[i-1],Rand_value);
			paillier_mul(pubkey,ciper_H[i],temp[i],ciper_G[i]);
		}

		paillier_mul(pubkey,ciper_PI[i],ciper_H[i],ciper_minus);	// PI
		paillier_exp(pubkey,ciper_PI[i],ciper_PI[i],Rand_value);
		// paillier_exp(pubkey,temp3[i],ciper_W[i],Rand_value);   // our IDEA
		// paillier_mul(pubkey,ciper_L[i],temp2[i], temp3[i]);	 // our IDEA
		
		paillier_mul(pubkey, ciper_L[i], ciper_PI[i], ciper_W[i]);
	}

	alpha = SRO2(ciper_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}

	//gmp_printf("alpha : %Zd(%d line)\n", paillier_dec(0, pubkey, prvkey, alpha), __LINE__);
	
	paillier_freeplaintext(Rand_value);		 
	paillier_freeciphertext(ciper_Rand_value);	paillier_freeciphertext(tmp);
	
	for(i=0; i<size; i++){
		paillier_freeciphertext(ciper_W[i]);	paillier_freeciphertext(ciper_H[i]);	paillier_freeciphertext(ciper_L[i]);	paillier_freeciphertext(ciper_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}
	
	return alpha;
}

paillier_ciphertext_t* protocol::Smin_basic2(paillier_ciphertext_t** ciper_R, paillier_ciphertext_t** ciper_L, paillier_ciphertext_t** ciper_M, paillier_ciphertext_t* alpha)
{
	int i;
	int flag = 0;

	// compute alpha
	for(i=0; i<size; i++)	{
		if (strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, ciper_L[i])), paillier_plaintext_to_str(plain_one)) == 0) {
			alpha = paillier_enc(0, pubkey, plain_one, paillier_get_rand_devurandom);
			//paillier_print("alpha : ", alpha);
			
			flag = 1;
			
			break;
		}
		else {
			alpha = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
			//paillier_print("alpha : ", alpha);
		}
	}

	if(flag == 1) {
		for(i=0; i<size; i++){
			paillier_exp(pubkey, ciper_M[i], ciper_R[i], plain_one);
		}
	}
	else {
		for(i=0; i<size; i++){
			paillier_exp(pubkey, ciper_M[i], ciper_R[i], plain_zero);
		}
	}

	//paillier_print("alpha : ", alpha);

	return alpha;
}

paillier_ciphertext_t** protocol::Smin_n(paillier_ciphertext_t*** ciper, int number){
	int i, j;
	int iter;
	int num = number;
	int k;
	paillier_ciphertext_t*** copy_ciper = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*num);
	for(i=0; i<num ; i++){
		copy_ciper[i] = ciper[i];
	}
	int e, q;
	

/*
	for(j=0;j<10;j++){
		for(i=0; i<10; i++){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, copy_ciper[j][i]));
		}
		printf("\n");
	}
		printf("\n");
*/
	paillier_ciphertext_t** sbd = SBD(paillier_create_enc_zero());
/*
	for(i=0;i<10;i++){
		gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, sbd[i]));
	}
	printf("\n");

	printf("\n");
	printf("\n");
*/
	double n = log10(num)/log10(2) - (int)(log10(num)/log10(2));
	
	
	if(n>0){
		n = (int)(log10(num)/log10(2))+1;
	}else{
		n = (int)(log10(num)/log10(2));
	}
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= num / 2; j++) {
			if (i == 1) {
				e = 2 * j - 2;
				q = 2 * j - 1;
				/*
				for(k=0; k<size; k++){
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, copy_ciper[e][k]));
				}
				printf("\n");
				for(k=0; k<size; k++){
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, copy_ciper[q][k]));
				}
				printf("\n");
				*/
				copy_ciper[e] = Smin_basic1(copy_ciper[e],copy_ciper[q]);
				// copy_ciper[e] = SMSn(copy_ciper[e],copy_ciper[q]);
				// SMSn(..) { 안에서 작은 것을 찾아서, return 해주는 것을 paillier cipertext **로 반환.. }
				copy_ciper[q] = sbd;
				//printf("%d %d\n",e,q);
			} else {
				e = (int)pow(2, i) * (j - 1);
				q = (int)pow(2, i) * j - (int)pow(2, i - 1);
				/*
				for(k=0; k<size; k++){
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, copy_ciper[e][k]));
				}
				printf("\n");
				for(k=0; k<size; k++){
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, copy_ciper[q][k]));
				}
				printf("\n");
				*/
				copy_ciper[e] = Smin_basic1(copy_ciper[e],copy_ciper[q]);
				// copy_ciper[e] = SMSn(copy_ciper[e],copy_ciper[q]);
				// SMSn(..) { 안에서 작은 것을 찾아서, return 해주는 것을 paillier cipertext **로 반환.. }
				copy_ciper[q] = sbd;
				//printf("%d %d\n",e,q);
			}
		}
		if (num % 2 == 0)
			num = num / 2;
		else
			num = num / 2 + 1;

		//printf("%dth iteration finished (out of %d iterations) -> %d %d\n",i, (int)n, e,q);
		//printf("%dth round finished (out of %d round)\n",i, (int)n);
	}
	
	free(sbd);
	return copy_ciper[0];
}

paillier_ciphertext_t** protocol::Smin_bool_n(paillier_ciphertext_t*** ciper, paillier_ciphertext_t** data, int number){
	paillier_ciphertext_t** index_arr;
	paillier_ciphertext_t** middle_result = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*number);

	for(int i = 0 ; i < number; i++){
		middle_result[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(middle_result[i]->c);
		/*
		for(int j = 0 ; j < size ; j++){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper[i][j]));
		}
		printf("\t");
		*/
	}
	//printf("\n");
	paillier_ciphertext_t** min = Smin_n(ciper, number);
/*
	printf("min_bit_arr : ");
	for(int i = 0 ; i < size ; i++){
		gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, min[i]));
	}
	printf("\n");
*/	
	paillier_ciphertext_t* tmp1 = paillier_create_enc(0);
	paillier_ciphertext_t* tmp2 = paillier_create_enc(0);
	int t = 0;
	

	for( int m = size ; m > 0 ; m-- ){
		t = (int)pow(2, m-1);
		tmp1 = paillier_create_enc(t);
		tmp1 = SM_p1(tmp1, min[size-m]);
		paillier_mul(pubkey, tmp2, tmp1, tmp2);
		paillier_freeciphertext(tmp1);
	}
	
	/*	print original data
	for( int i = 0 ; i < number; i++){
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, data[i]));
	}
	*/

	//printf("\n");
	//gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, tmp2));

	//printf("\n");	
	for(int i = 0 ; i < number ; i++){
		paillier_subtract(pubkey, middle_result[i], data[i], tmp2);
		//gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, middle_result[i]));
	}
	//printf("\n");
	
	index_arr = Smin_bool_n_sub(middle_result, number);

		
/*
	for(int i = 0 ; i < number ; i++){
		gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, index_arr[i]));
	}
	printf("\n");
*/
	return index_arr;
}



paillier_ciphertext_t** protocol::Smin_bool_n_sub(paillier_ciphertext_t** ciper, int number){
	paillier_ciphertext_t** index_arr = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*number);
	for( int i = 0 ; i < number ; i++){
		index_arr[i] = paillier_create_enc(0);
	}
	for( int i = 0 ; i < number ; i++){
		if (strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, ciper[i])), paillier_plaintext_to_str(plain_zero)) == 0){
			index_arr[i] = paillier_enc(0, pubkey, plain_one, paillier_get_rand_devurandom);
			break;
		}
	}
/*
	for(int i = 0 ; i < number ; i++){
		gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, index_arr[i]));
	}
	printf("\n");
*/
	return index_arr;
}
// added by KHI. 150711
paillier_ciphertext_t* protocol::SRO(paillier_ciphertext_t*** qLL, paillier_ciphertext_t*** qRR, paillier_ciphertext_t*** nodeLL, paillier_ciphertext_t*** nodeRR)
{
	//printf("\n===== Now SRO starts ===== \n");
	
	int i, j = 0;
	
	bool func = false;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* ciper_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** ciper_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_G	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);
	paillier_ciphertext_t* final_alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(final_alpha->c);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	for(i=0; i<size+1; i++){
		ciper_W[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_G[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_H[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_L[i]	= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_M[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		ciper_PI[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp[i]	= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp2[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));

		mpz_init(ciper_W[i]->c);
		mpz_init(ciper_G[i]->c);
		mpz_init(ciper_H[i]->c);
		mpz_init(ciper_L[i]->c);
		mpz_init(ciper_M[i]->c);
		mpz_init(ciper_PI[i]->c);
		mpz_init(temp[i]->c);
		mpz_init(temp2[i]->c);
	}

	// q.LL 및 node.RR에 대해서 수행
	for(i=0; i<dim; i++) {
		//printf("\n%dth dimension\n", i);
		
		for(j=0; j<size+1; j++) {
			// compute W
			if(func){	// true :  F : u>v	
				paillier_subtract(pubkey, ciper_W[j], qLL[i][j], SM_p1(qLL[i][j], nodeRR[i][j]));	
			}else{
				paillier_subtract(pubkey, ciper_W[j], nodeRR[i][j], SM_p1(qLL[i][j], nodeRR[i][j]));
			}

			// compute G
			ciper_G[j]=SBXOR(qLL[i][j],nodeRR[i][j]);

			// compute H
			if(j==0){
				paillier_exp(pubkey,ciper_H[j], ciper_zero, Rand_value);
				paillier_mul(pubkey,ciper_H[j], ciper_H[j], ciper_G[j]);
			}else{
				paillier_exp(pubkey,temp[j],ciper_H[j-1],Rand_value);
				paillier_mul(pubkey,ciper_H[j],temp[j],ciper_G[j]);
			}
		}

		for(j=0; j<size+1; j++){
			// compute PI
			paillier_mul(pubkey, ciper_PI[j], ciper_H[j], ciper_minus);	// PI

			// compute L with enhanced privacy (new IDEA)
			mpz_init(temp[i]->c);
			paillier_exp(pubkey,temp[j],ciper_PI[j],Rand_value);	// randomize PI
			paillier_exp(pubkey,temp2[j],ciper_W[j],Rand_value); // randomize W
			paillier_mul(pubkey,ciper_L[j],temp[j], temp2[j]);
		}
	
		alpha = SRO2(ciper_L, alpha);

		if(func){
			alpha = SBN(alpha);
		}
	
		// AND operation (using SM Protocol) 
		if(i == 0)
			final_alpha = SM_p1(ciper_one, alpha);	// 최초에는 1과 AND 연산을 해야함. ?箕坪繭捉?0이 나오면 0으로 바뀌는 논리
		else
			final_alpha = SM_p1(final_alpha, alpha);
	}


	// q.RR 및 node.LL에 대해서 수행
	for(i=0; i<dim; i++) {
		//printf("\n%dth dimension\n", i);
		
		for(j=0; j<size+1; j++) {
			// compute W
			if(func){	// true :  F : u>v	
				paillier_subtract(pubkey, ciper_W[j], nodeLL[i][j], SM_p1(nodeLL[i][j], qRR[i][j]));	
			}else{
				paillier_subtract(pubkey, ciper_W[j], qRR[i][j], SM_p1(nodeLL[i][j], qRR[i][j]));
			}

			// compute G
			ciper_G[j]=SBXOR(nodeLL[i][j], qRR[i][j]);

			// compute H
			if(j==0){
				paillier_exp(pubkey,ciper_H[j], ciper_zero, Rand_value);
				paillier_mul(pubkey,ciper_H[j], ciper_H[j], ciper_G[j]);
			}else{
				paillier_exp(pubkey,temp[j],ciper_H[j-1],Rand_value);
				paillier_mul(pubkey,ciper_H[j],temp[j],ciper_G[j]);
			}
		}

		for(j=0; j<size+1; j++){
			// compute PI
			paillier_mul(pubkey, ciper_PI[j], ciper_H[j], ciper_minus);	// PI

			// compute L with enhanced privacy (new IDEA)
			mpz_init(temp[i]->c);
			paillier_exp(pubkey,temp[j],ciper_PI[j],Rand_value);	// randomize PI
			paillier_exp(pubkey,temp2[j],ciper_W[j],Rand_value); // randomize W
			paillier_mul(pubkey,ciper_L[j],temp[j], temp2[j]);
		}
		
		alpha = SRO2(ciper_L, alpha);

		if(func){
			alpha = SBN(alpha);
		}
		
		final_alpha = SM_p1(final_alpha, alpha);

	}

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(ciper_Rand_value);		paillier_freeciphertext(tmp);

	paillier_freeciphertext(alpha);	
	
	for(i=0; i<size+1; i++){
		paillier_freeciphertext(ciper_W[i]);	paillier_freeciphertext(ciper_H[i]);	paillier_freeciphertext(ciper_L[i]);	paillier_freeciphertext(ciper_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}
	return final_alpha;	
}

// added by KHI. 150711
paillier_ciphertext_t*  protocol::SRO2(paillier_ciphertext_t** ciper_L, paillier_ciphertext_t* alpha)
{
	int i;

	// compute alpha
	for(i=0; i<size; i++)	{
		if (strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, ciper_L[i])), paillier_plaintext_to_str(plain_zero)) == 0) {
			alpha = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
			//paillier_print("alpha 0 : ", alpha);
			break;
		}
		else {
			alpha = paillier_enc(0, pubkey, plain_one, paillier_get_rand_devurandom);
			//paillier_print("alpha 1 : ", alpha);
		}
	}

	//paillier_print("alpha : ", alpha);
	
	return alpha;
}

void protocol::makeGate()
{
	// input
	int BitSize = size+2; 
	int NumInputs = 4;
	int InputRange = 0;
	
	vector<USHORT> Ports;
	Ports.resize(1, 8888);
	vector<string> Addrs;
	Addrs.resize(1, "localhost");
	
	int NumParties = 1;
	int PID = 0;
	string tempstr;
	tempstr = "8936097950764538541647693880373941060412422053330581106416547022143872696986309618392492140173286901378451169263228782167065246660381840799235731486130087";
	istringstream pis(tempstr);
	ZZ p;
	pis >> p;
	tempstr = "7660915846360464914746169975675514711978996378800371841846530440188167304050642048045944447898172503094027848348928324282770185052123039167705188057761352";
	istringstream gis(tempstr);
	ZZ g;
	gis >> g;
	tempstr = "2323967942723790652936980"; // seed for server
	string CircCreateName = "smsn";
	vector<int> CircCreateParam(2);
	CircCreateParam[0] = NumInputs;
	CircCreateParam[1] = BitSize;
	// { NumInputs, BitSize }; // input, reps
	pConfig->SetAddrPID(Addrs);
	pConfig->SetPortPID(Ports);
	pConfig->SetNumParties(NumParties);
	pConfig->SetNumInputs(NumInputs);
	pConfig->SetPID(PID);
	pConfig->SetP(p);
	pConfig->SetG(g);
	pConfig->SetSeed(tempstr);
	pConfig->SetCircCreateName(CircCreateName);
	pConfig->SetCircCreateParam(CircCreateParam);

	pCircuit = CREATE_CIRCUIT(pConfig->GetNumParties(), pConfig->GetCircCreateName(),
		pConfig->GetCircCreateParams());
}

int protocol::G_CMP(int mr1, int m1r1, int mr2, int m2r2) // -r1, m1+r1, -r2, m2+r2
{
	double tc1 = clock();

	vector<int> InputDatas(4);
	// input setting
	InputDatas[0] = mr1;
	InputDatas[1] = m1r1;
	InputDatas[2] = mr2;
	InputDatas[3] = m2r2;
	///cout << "Input Set" << endl;

	// circtool get
	GATE* gates = pCircuit->Gates();
	gates[0].val = 0;
	gates[1].val = 1;
	//cout << "Gates Create" << endl;

	for (unsigned i = 0; i<4; i++) // 4 input
	{
		int bits = pCircuit->GetNumVBits(i);
		int start = pCircuit->GetInputStart(i);
		int end = pCircuit->GetInputEnd(i);
		//cout << i << "th Start Gate Num : " << start << ", End Gate Num : " << end << endl;
		
		int j = start;
		for (int j = start; j <= end; j++)
		{
			for (int k = 0; k<bits; k++)
			{
				int mask = (1 << k);
				gates[j].val = (char)!!(InputDatas[i] & mask);
				j++;
			}
		}
	}

	pCircuit->Evaluate();
	pCircuit->Save("SMSnCirc.txt", TRUE);

	int o_start = pCircuit->GetOutputStart(0);

	return gates[o_start].val; // if a > b, then return 1 else return 0
	// if 1
	// then return first *
	// if 0
	// then return second *
}





paillier_ciphertext_t*  protocol::SSR(paillier_ciphertext_t* ciper){
	int rand = 6;
	paillier_plaintext_t* temp = paillier_plaintext_from_ui(rand);
	rand = rand * rand;
	printf("rand : %d\n", rand);

	paillier_plaintext_t* tmp = paillier_plaintext_from_ui(rand);
	paillier_ciphertext_t* ciper_rand = paillier_create_enc(rand);

	paillier_exp(pub, ciper, ciper, tmp);
	paillier_print("mul : ", ciper);
	ciper = SSR_sub(ciper);
	paillier_print("root : ", ciper);
	paillier_divide(pub, ciper, ciper, temp);
	paillier_print("divide : ", ciper);

	free(temp);
	free(tmp);
	free(ciper_rand);
	
	return ciper;
}

paillier_ciphertext_t*  protocol::SSR_sub(paillier_ciphertext_t* ciper){
	paillier_plaintext_t* tmp = paillier_dec(0, pub, prv, ciper);
	int temp = mpz_get_si(tmp->m);
	mpz_sqrt(tmp->m, tmp->m);
	ciper = paillier_enc(0, pub, tmp, paillier_get_rand_devurandom);
	return ciper;
}

void protocol::protocol_1_free(paillier_ciphertext_t* cipher, bool print){
	if(print)
		printf("----protocol_1_free----\n");
	free(cipher);
	cipher = NULL;
}
void protocol::protocol_2_free(paillier_ciphertext_t** cipher, int x, bool print){
	if(print)
		printf("----protocol_2_free----\n");
	int i = 0;
	for( i = 0 ; i < x ; i++ ){
		free(cipher[i]);
		cipher[i] = NULL;
	}
	cipher = NULL;
}
void protocol::protocol_3_free(paillier_ciphertext_t*** cipher, int x, int y, bool print){
	if(print)
		printf("----protocol_3_free----\n");
	int i = 0, j = 0;
	for( i = 0 ; i < x ; i++ ){
		for( j = 0 ; j < y ; j++ ){
			free(cipher[i][j]);
			cipher[i][j] = NULL;
		}
		free(cipher[i]);
		cipher[i] = NULL;
	}
	free(cipher);
	cipher = NULL;
}
void protocol::protocol_4_free(paillier_ciphertext_t**** cipher, int x, int y, int z, bool print){
	if(print)
		printf("----protocol_4_free----\n");
	int i = 0, j = 0, m = 0;
	for( i = 0 ; i < x ; i++ ){
		for( j = 0 ; j < y ; j++ ){
			for( m = 0 ; m < z ; m++ ){
				free(cipher[i][j][m]);
				cipher[i][j][m] = NULL;
			}
			free(cipher[i][j]);
			cipher[i][j] = NULL;
		}
		free(cipher[i]);
		cipher[i] = NULL;
	}
	free(cipher);
	cipher = NULL;
}

void protocol::protocol_1_printf(paillier_ciphertext_t* cipher, char* str, bool print){
	if(print){
		printf("%s : ", str);
		gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, cipher));
	}
}
void protocol::protocol_2_printf(paillier_ciphertext_t** cipher, int x, char* str, bool print){
	if(print){
		int i = 0;
		printf("%s : ", str);
		for( i = 0 ; i < x ; i++ ){
			gmp_printf(" %Zd \t", paillier_dec(0, pubkey, prvkey, cipher[i]));
		}
		printf("\n");
	}
}
void protocol::protocol_3_printf(paillier_ciphertext_t*** cipher, int x, int y, char* str, bool print){
	if(print){
		int i = 0, j = 0;
		printf("===%s===\n", str);
		for( i = 0 ; i < x ; i++ ){
			printf("%d : \t", i);
			for( j = 0 ; j < y ; j++ ){
				gmp_printf(" %Zd \t", paillier_dec(0, pubkey, prvkey, cipher[i][j]));
			}
			printf("\n");
		}
		printf("\n");
	}
}

paillier_ciphertext_t* protocol::SC(paillier_ciphertext_t* u, paillier_ciphertext_t* v){
	paillier_plaintext_t* a = paillier_dec(0, pubkey, prvkey, u);
	paillier_plaintext_t* b = paillier_dec(0, pubkey, prvkey, v);
	paillier_ciphertext_t* r;


	long w; 
	long x;
	long y;
	paillier_plaintext_t* p_y = paillier_plaintext_from_ui(0);

	paillier_plaintext_t* shift = paillier_plaintext_from_ui(2);
	paillier_plaintext_t* tmp = paillier_plaintext_from_ui(0);
	mpz_div(tmp->m,pubkey->n,shift->m);

	mpz_mod(a->m, a->m, pubkey->n);
	//gmp_printf("L : %Zd\n", a->m);
	long wi = mpz_cmp(a->m , tmp->m);
	
	mpz_mod(b->m, b->m, pubkey->n);
	//gmp_printf("R : %Zd\n", b->m);
	long xi = mpz_cmp(b->m , tmp->m);

	mpz_sub(p_y->m, b->m, a->m);
	mpz_mod(p_y->m, p_y->m, pubkey->n);
	long yi = mpz_cmp(p_y->m , tmp->m);
	
	if(wi <= 0){
		w = 1;
	}else{
		w = 0;
	}

	if(xi <= 0){
		x = 1;
	}else{
		x = 0;	
	}
	
	if(yi <= 0){
		y = 1;
	}else{
		y = 0;	
	}

	//gmp_printf("%d %d %d\n", w, x, y);

	if( w == 0 && x == 1 ){
		r = paillier_create_enc(0);
	}else if( w == 1 && x == 0 ){
		r = paillier_create_enc(1);
	}else if( w == 0 && x == 0 && y == 0){
		r = paillier_create_enc(0);
	}else if( w == 0 && x == 0 && y == 1){
		r = paillier_create_enc(1);
	}else if( w == 1 && x == 1 && y == 0){
		r = paillier_create_enc(0);
	}else if( w == 1 && x == 1 && y == 1){
		r = paillier_create_enc(1);
	}

	//gmp_printf("SC!!! %Zd :: %Zd == %Zd\n", paillier_dec(0, pubkey, prvkey, u), paillier_dec(0, pubkey, prvkey, v), paillier_dec(0, pubkey, prvkey, r));
	
	return r;
}

paillier_ciphertext_t** protocol::Smin_bool(paillier_ciphertext_t** data, int number){
	//printf("Start Smin_bool ");
	paillier_ciphertext_t** index_arr = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	
	paillier_plaintext_t** a = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*k);
	paillier_plaintext_t** b = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*k);
	
	for(int i = 0 ; i < k ; i++){
		a[i] = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
		a[i] = paillier_dec(0, pubkey, prvkey, data[i]);
	}

	int min_idx = 0;
	paillier_ciphertext_t* res;
	for(int i = 1 ; i < k ; i++){
		res = SC(data[min_idx], data[i]);
		if (strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, res)), paillier_plaintext_to_str(plain_zero)) == 0) {
			min_idx = i;
		}
	}

	for(int i = 0 ; i < k ; i++){
		if(i == min_idx){
			index_arr[i] = paillier_create_enc(1);
		}else{
			index_arr[i] = paillier_create_enc(0);
		}
	}
	
	for(int i = 0 ; i < k ; i++){
		free(a[i]);
	}
	free(a);free(b);

	//printf("End Smin_bool\n");

	return index_arr;
}

paillier_ciphertext_t** protocol::G_Smin_bool(paillier_ciphertext_t** data, int number){
	//printf("Start Smin_bool ");
	paillier_ciphertext_t** index_arr = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	
	paillier_plaintext_t** a = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*k);
	paillier_plaintext_t** b = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*k);
	
	for(int i = 0 ; i < k ; i++){
		a[i] = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
		a[i] = paillier_dec(0, pubkey, prvkey, data[i]);
	}

	int min_idx = 0;
	paillier_ciphertext_t* res;
	for(int i = 1 ; i < k ; i++){
		res = GSCMP(data[min_idx], data[i]);
		if (strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, res)), paillier_plaintext_to_str(plain_zero)) == 0) {
			min_idx = i;
		}
	}

	for(int i = 0 ; i < k ; i++){
		if(i == min_idx){
			index_arr[i] = paillier_create_enc(1);
		}else{
			index_arr[i] = paillier_create_enc(0);
		}
	}
	
	for(int i = 0 ; i < k ; i++){
		free(a[i]);
	}
	free(a);free(b);

	//printf("End Smin_bool\n");

	return index_arr;
}

paillier_ciphertext_t* protocol::AS_CMP(paillier_ciphertext_t* u, paillier_ciphertext_t* v)
{
	return AS_CMP_sub(u, v);
}


paillier_ciphertext_t* protocol::AS_CMP_sub(paillier_ciphertext_t* u, paillier_ciphertext_t* v)
{
	paillier_ciphertext_t* AS_CMP_RESULT;	

	paillier_plaintext_t* plain_u = paillier_dec(0, pubkey, prvkey, u);
	paillier_plaintext_t* plain_v = paillier_dec(0, pubkey, prvkey, v);
		
	if ( mpz_cmp(plain_u->m, plain_v->m) <=0 )
	{
		AS_CMP_RESULT = paillier_create_enc(1);
	}
	else 
	{
		AS_CMP_RESULT = paillier_create_enc(0);
	}
	return AS_CMP_RESULT;
}
