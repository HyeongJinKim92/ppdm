#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include "ppdm.h"

using namespace std;


char shellquery[128]; 		
protocol proto;
setting setting;


void ppdm_main(char* argv[])
{
	parsed_query * user_query = parsing_query(argv);
	cout << user_query->query << endl;
	proto.set_Param(user_query);
	
	range_main(user_query);

}

int range_main(parsed_query * user_query)
{
	int modulus = user_query->key_s;
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;


	paillier_keygen(user_query->key_s, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub,prv, user_query->key_s);

	int NumNode=0;

	boundary q;
	
		
	// query setting start
	q.LL =  (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*user_query->dim_n);
	q.RR =  (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*user_query->dim_n);

	// query setting end

	for(int i=0; i<user_query->dim_n; i++){
		q.LL[i] = paillier_create_enc(user_query->qLb[i]);
		paillier_print("LL : ", q.LL[i]);
		q.RR[i] = paillier_create_enc(user_query->qUb[i]);
		paillier_print("RR : ", q.RR[i]);
	}

	char kd_filename[128]; 		
	boundary* node = 0;
	
	if(	user_query->rquery_t != RANGE_B )
	{
		sprintf(kd_filename, "input/kd_range/KD_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
		node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);	
	}
		

	char inputFilename[128]; 
	sprintf(inputFilename, "input/range/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** ciper = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);
	
	printf("\n===== sRange start =====\n");
	int result_num=0;
	int** result = 0;
	
	if(user_query->rquery_t == RANGE_B)			result = proto.sRange_M(ciper, q, node, user_query->data_n, NumNode, &result_num); //rangeB
	else if(user_query->rquery_t == RANGE_I) 	result = proto.sRange_I(ciper, q, node, user_query->data_n, NumNode, &result_num);  //rangeI
	else if(user_query->rquery_t == RANGE_GI) 	result = proto.sRange_G(ciper, q, node, user_query->data_n, NumNode, &result_num); //rangeGI
//	else if(user_query->rquery_t == RANGE_PB) 	result = proto.sRange_PM(ciper, q, node, user_query->data_n, NumNode, &result_num); //rangePGI
//	else if(user_query->rquery_t == RANGE_PGI) 	result = proto.sRange_PGI(ciper, q, node, user_query->data_n, NumNode, &result_num); //rangePGI
//	else if(user_query->rquery_t == RANGE_PAI) 	result = proto.sRange_PAI(ciper, q, node, user_query->data_n, NumNode, &result_num); //rangePAI

	if(user_query->rquery_t == RANGE_B)			setting.TimeResult_write_int(shellquery, result, result_num, "RANGE/RANGE_B", proto); //rangeB
	else if(user_query->rquery_t == RANGE_I) 	setting.TimeResult_write_int(shellquery, result, result_num, "RANGE/RANGE_I", proto);  //rangeI
	else if(user_query->rquery_t == RANGE_GI) 	setting.TimeResult_write_int(shellquery, result, result_num, "RANGE/RANGE_GI", proto); //rangeGI
	else if(user_query->rquery_t == RANGE_PB) 	setting.TimeResult_write_int(shellquery, result, result_num, "RANGE/RANGE_PB", proto); //rangePGI
	else if(user_query->rquery_t == RANGE_PGI) 	setting.TimeResult_write_int(shellquery, result, result_num, "RANGE/RANGE_PGI", proto); //rangePAI
	else if(user_query->rquery_t == RANGE_PAI) 	setting.TimeResult_write_int(shellquery, result, result_num, "RANGE/RANGE_PAI", proto); //rangePAI

	for( i = 0 ; i < result_num ; i++ ){
		printf("%d result : ", (i+1) );
		for ( j = 0 ; j < user_query->dim_n ; j++ ){
			printf("%d \t", result[i][j]);
		}
		printf("\n");
	}
	
	//free
	for(i=0; i<2; i++) {
		for(j=0; j<user_query->dim_n; j++) {
			paillier_freeciphertext(node[i].LL[j]);	
			paillier_freeciphertext(node[i].RR[j]);
		}
	}
	
	for(j=0; j<user_query->dim_n; j++) {
			paillier_freeciphertext(q.LL[j]);		
			paillier_freeciphertext(q.RR[j]);
	}

	free(node);
	proto.protocol_free();

	return 0;
}
int topk_main(parsed_query* user_query)
{
	return 0;
}
int knn_main(parsed_query* user_query)
{
	return 0;
}
int classification_main(parsed_query* user_query)
{
	return 0;
}
int kmeans_main(parsed_query* user_query)
{
	return 0;
}

parsed_query* parsing_query(char* argv[])
{
	parsed_query* query = new parsed_query;
	
	query->query = (QUERY)atoi(argv[1]);
	cout << query->query << endl;

	if(query->query == RANGE)				query->rquery_t = (RANGE_TYPE)atoi(argv[2]);
	else if(query->query == TOPK)			query->tquery_t = (TOPK_TYPE)atoi(argv[2]);
	else if(query->query == KNN)				query->kquery_t = (KNN_TYPE)atoi(argv[2]);
	else if(query->query == CLASSIFICATION)	query->cquery_t = (CLASSIFICATION_TYPE)atoi(argv[2]);
	else if(query->query == KMEANS)			query->mquery_t = (KMEANS_TYPE)atoi(argv[2]);
	else if(query->query == ASSOCATION)		{}
	else if(query->query == TEST)			{}


	query->tree_level = atoi(argv[3]);
	query->thread_n = atoi(argv[4]);
	query->data_n = atoi(argv[5]);
	query->dim_n = atoi(argv[6]);
	query->bit_s = atoi(argv[7]);
	query->k = atoi(argv[8]);
	query->key_s = atoi(argv[9]);

	if(query->query == RANGE)
	{
		query->qLb = new int[query->dim_n];
		query->qUb = new int[query->dim_n];
		for ( int i = 0 ; i < query->dim_n ; i++ )
		{
			query->qLb[i] = atoi(argv[10+i]);
			query->qUb[i] = atoi(argv[10+(query->dim_n)+i]);
		}
	}	
	else if(query->query == TOPK)
	{
		query->score = new int[query->dim_n+1];
		for (int i = 0; i < query->dim_n+1 ; i++)
		{
			query->score[i] = atoi(argv[10+i]);
		}
	}	
	else if(query->query == KNN)				
	{
		query->q = new int[query->dim_n];
		for (int i = 0; i < query->dim_n ; i++)
		{
			query->q[i] = atoi(argv[10+i]);
		}
	}
	else if(query->query == CLASSIFICATION)	
	{
		query->q = new int[query->dim_n];
		for (int i = 0; i < query->dim_n ; i++)
		{
			query->q[i] = atoi(argv[10+i]);
		}
	}	
	else if(query->query == KMEANS)			{}
	else if(query->query == ASSOCATION)		{}
	else if(query->query == TEST)			{}
	
	return query;
}

void print_query()
{
	
}
