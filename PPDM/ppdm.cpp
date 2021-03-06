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
	proto.set_Param(user_query);
	if(user_query->query == RANGE)					range_main(user_query);
	else if(user_query->query == TOPK)				topk_main(user_query);
	else if(user_query->query == KNN)				knn_main(user_query);
	else if(user_query->query == CLASSIFICATION)	classification_main(user_query);
	else if(user_query->query == KMEANS)			kmeans_main(user_query);
}

int range_main(parsed_query * user_query)
{
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
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);
	
	printf("\n===== sRange start =====\n");
	int result_num=0;
	int** result = 0;
	
	if(user_query->rquery_t == RANGE_B)			result = proto.sRange_M(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangeB
	else if(user_query->rquery_t == RANGE_I) 	result = proto.sRange_I(cipher, q, node, user_query->data_n, NumNode, &result_num);  //rangeI
	else if(user_query->rquery_t == RANGE_GI) 	result = proto.sRange_G(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangeGI
//	else if(user_query->rquery_t == RANGE_PB) 	result = proto.sRange_PM(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangePGI
//	else if(user_query->rquery_t == RANGE_PGI) 	result = proto.sRange_PGI(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangePGI
//	else if(user_query->rquery_t == RANGE_PAI) 	result = proto.sRange_PAI(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangePAI

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
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;

	paillier_keygen(user_query->key_s, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub, prv, user_query->key_s);

	paillier_ciphertext_t** cipher_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(user_query->dim_n)+1); 

	paillier_ciphertext_t** max_val = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*user_query->dim_n); 

	paillier_ciphertext_t* cipher_zero = paillier_create_enc(0);

	int MAX_VAL = 0;
	MAX_VAL = sqrt(pow(2, user_query->bit_s)/user_query->dim_n);
	MAX_VAL = 512;
	for( i = 0 ; i < user_query->dim_n ; i++ ){
		if(user_query->score[i] < 0) {
			cipher_query[i] = paillier_create_enc(user_query->score[i]*(-1));	// 양수로 바꿔서 일단 암호화 한 후, 
			paillier_subtract(pub, cipher_query[i] , cipher_zero , cipher_query[i]);	 // 0에서 빼서 음수로 전환
		}
		else {
			cipher_query[i] = paillier_create_enc(user_query->score[i]);
		}
		max_val[i] = paillier_create_enc(MAX_VAL);
		gmp_printf("max val : %Zd \n", paillier_dec(0, pub, prv, max_val[i]));
	}
	

	int tmp=abs(user_query->score[0]);
	for(int i=1 ; i < user_query->dim_n ; i++)
	{
		if( tmp < abs(user_query->score[i]) )
			tmp=abs(user_query->score[i]);
	}

	cipher_query[user_query->dim_n] = paillier_create_enc(tmp); 

	//gmp_printf("hint %Zd \n", paillier_dec(0, pub, prv, cipher_query[dim]));

	cout<<"user_query->key_s : " <<user_query->key_s<<endl;

	printf("query : ");

	for(i = 0 ; i < user_query->dim_n + 1 ; i++ ){				
		if(i == user_query->dim_n)	// hint
			printf("hint : ");
			gmp_printf("%Zd \t", paillier_dec(0, pub, prv, cipher_query[i]));
		cout<<endl;
	}
	printf("\n");

	int NumNode = 0;
	char kd_filename[128]; 
	boundary* node = 0;


	if(	user_query->tquery_t != TOPK_B )
	{
		sprintf(kd_filename, "input/kd_topk/KD_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
		node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);	
	}

	char inputFilename[128]; 
	sprintf(inputFilename, "input/topk/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);



	printf("\n===== Top-k start =====\n");

	int result_num=0;
	int** topk = 0;

	
	if(user_query->tquery_t == TOPK_B)			topk = proto.STopk_M(cipher, cipher_query, user_query->data_n); //TOPKB
	else if(user_query->tquery_t == TOPK_I) 	topk =  proto.STopk_I(cipher, cipher_query, node, max_val, user_query->data_n, NumNode);  //TOPKI
	else if(user_query->tquery_t == TOPK_GI) 	topk = proto.STopk_G(cipher, cipher_query, node, max_val, user_query->data_n, NumNode); //TOPKGI
//	else if(user_query->tquery_t == TOPK_PB) 	result = proto.sRange_PM(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); //TOPKPGI
//	else if(user_query->tquery_t == TOPK_PGI) 	result = proto.sRange_PGI(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); //TOPKPGI
//	else if(user_query->tquery_t == TOPK_PAI) 	result = proto.sRange_PAI(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); //TOPKPAI

	if(user_query->tquery_t == TOPK_B)			setting.TimeResult_write_int(shellquery, topk, result_num, "TOPK/TOPK_B", proto); //TOPKB
	else if(user_query->tquery_t == TOPK_I) 	setting.TimeResult_write_int(shellquery, topk, result_num, "TOPK/TOPK_I", proto);  //TOPKI
	else if(user_query->tquery_t == TOPK_GI) 	setting.TimeResult_write_int(shellquery, topk, result_num, "TOPK/TOPK_GI", proto); //TOPKGI
	else if(user_query->tquery_t == TOPK_PB) 	setting.TimeResult_write_int(shellquery, topk, result_num, "TOPK/TOPK_PB", proto); //TOPKPGI
	else if(user_query->tquery_t == TOPK_PGI) 	setting.TimeResult_write_int(shellquery, topk, result_num, "TOPK/TOPK_PGI", proto); //TOPKPAI
	else if(user_query->tquery_t == TOPK_PAI) 	setting.TimeResult_write_int(shellquery, topk, result_num, "TOPK/TOPK_PAI", proto); //TOPKPAI

	for( i = 0 ; i < user_query->k ; i++ ){
		printf("%d result : ", (i+1));
		for ( j = 0 ; j < user_query->dim_n ; j++ ){
			printf("%d \t", topk[i][j]);
		}
		printf("\n");
	}

	proto.protocol_free();
	return 0;
}


int knn_main(parsed_query* user_query)
{
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;


	paillier_keygen(user_query->key_s, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub, prv, user_query->key_s);


	paillier_ciphertext_t** cipher_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*user_query->dim_n); 
	for( i = 0 ; i < user_query->dim_n ; i++ ){
		cipher_query[i] = paillier_create_enc(user_query->q[i]);
	}

	cout<<"user_query->key_s : " <<user_query->key_s<<endl;

	printf("\nquery : ");
	for(i = 0 ; i < user_query->dim_n ; i++ ){				
		gmp_printf("%Zd \t", paillier_dec(0, pub, prv, cipher_query[i]));
	}
	printf("\n\n");

	int NumNode = 0;
	char kd_filename[128]; 
	boundary* node = 0;


	if(	user_query->kquery_t != KNN_B )
	{
		sprintf(kd_filename, "input/kd_knn/KD_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
		node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);
	}

	char inputFilename[128]; 
	sprintf(inputFilename, "input/knn/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);

	int result_num=0;
	int** SkNNm = 0;

	if(user_query->kquery_t == KNN_B)			SkNNm = proto.SkNN_M(cipher, cipher_query, user_query->k, user_query->data_n); //KNNB
	else if(user_query->kquery_t == KNN_I) 		SkNNm = proto.SkNN_I(cipher, cipher_query, node, user_query->k, user_query->data_n, NumNode);  //KNNI
	else if(user_query->kquery_t == KNN_GI) 	SkNNm = proto.SkNN_G(cipher, cipher_query, node, user_query->k, user_query->data_n, NumNode); //KNNGI
//	else if(user_query->kquery_t == KNN_PB) 	SkNNm = proto.sRange_PM(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); //KNNPGI
//	else if(user_query->kquery_t == KNN_PGI) 	SkNNm = proto.sRange_PGI(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); //KNNPGI
//	else if(user_query->kquery_t == KNN_PAI) 	SkNNm = proto.sRange_PAI(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); //KNNPAI

	if(user_query->kquery_t == KNN_B)			setting.TimeResult_write_int(shellquery, SkNNm, result_num, "KNN/KNN_B", proto); //KNNB
	else if(user_query->kquery_t == KNN_I) 		setting.TimeResult_write_int(shellquery, SkNNm, result_num, "KNN/KNN_I", proto);  //KNNI
	else if(user_query->kquery_t == KNN_GI) 	setting.TimeResult_write_int(shellquery, SkNNm, result_num, "KNN/KNN_GI", proto); //KNNGI
	else if(user_query->kquery_t == KNN_PB) 	setting.TimeResult_write_int(shellquery, SkNNm, result_num, "KNN/KNN_PB", proto); //KNNPGI
	else if(user_query->kquery_t == KNN_PGI) 	setting.TimeResult_write_int(shellquery, SkNNm, result_num, "KNN/KNN_PGI", proto); //KNNPAI
	else if(user_query->kquery_t == KNN_PAI) 	setting.TimeResult_write_int(shellquery, SkNNm, result_num, "KNN/KNN_PAI", proto); //KNNPAI

	for( i = 0 ; i < user_query->k ; i++ ){
		printf("%d result : ", (i+1) );
		for ( j = 0 ; j < user_query->dim_n ; j++ ){
			printf("%d \t", SkNNm[i][j]);
		}
		printf("\n");
	}

	proto.protocol_free();
	return 0;
}
int classification_main(parsed_query* user_query)
{
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;
	
	paillier_keygen(user_query->key_s, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub, prv, user_query->key_s);

	paillier_ciphertext_t** cipher_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*user_query->dim_n); 
	for( i = 0 ; i < user_query->dim_n ; i++ ){
		cipher_query[i] = paillier_create_enc(user_query->q[i]);
	}

	cout<<"user_query->key_s : " <<user_query->key_s<<endl;

	printf("\nquery : ");
	for(i = 0 ; i < user_query->dim_n ; i++ ){				
		gmp_printf("%Zd \t", paillier_dec(0, pub, prv, cipher_query[i]));
	}
	printf("\n\n");


	int NumNode = 0;
	boundary* node = 0;


	if(	user_query->cquery_t != CLASSIFICATION_B )
	{
		char kd_filename[128]; 
		sprintf(kd_filename, "input/kd_classification/KD_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
		node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);
	}

	char inputFilename[128]; 
	sprintf(inputFilename, "input/classification/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);
	
	////////Label 설정	

	int Entire_num = 18;
	paillier_ciphertext_t** Entire_set = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(Entire_num));
	for( i = 0 ; i < Entire_num ; i++){
		Entire_set[i] = paillier_create_enc(i);
	}
	int result_num=0;
	int** SkNNm = 0;

	if(user_query->cquery_t == CLASSIFICATION_B)			SkNNm = proto.Classification_M(cipher, cipher_query, Entire_set, user_query->k, user_query->data_n, Entire_num);	
	else if(user_query->cquery_t == CLASSIFICATION_I) 		SkNNm = proto.Classification_I(cipher, cipher_query, Entire_set, node, user_query->k, user_query->data_n, NumNode, Entire_num);	
	else if(user_query->cquery_t == CLASSIFICATION_GI) 		SkNNm = proto.Classification_G(cipher, cipher_query, Entire_set, node, user_query->k, user_query->data_n, NumNode, Entire_num);
//	else if(user_query->cquery_t == CLASSIFICATION_PB) 		SkNNm = proto.sRange_PM(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); 
//	else if(user_query->cquery_t == CLASSIFICATION_PGI) 	SkNNm = proto.sRange_PGI(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); 
//	else if(user_query->cquery_t == CLASSIFICATION_PAI) 	SkNNm = proto.sRange_PAI(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); 

	if(user_query->cquery_t == CLASSIFICATION_B)			setting.TimeResult_write_int(shellquery, SkNNm, result_num, "CLASSIFICATION/CLASSIFICATION_B", proto); 
	else if(user_query->cquery_t == CLASSIFICATION_I) 		setting.TimeResult_write_int(shellquery, SkNNm, result_num, "CLASSIFICATION/CLASSIFICATION_I", proto); 
	else if(user_query->cquery_t == CLASSIFICATION_GI) 		setting.TimeResult_write_int(shellquery, SkNNm, result_num, "CLASSIFICATION/CLASSIFICATION_GI", proto);
	else if(user_query->cquery_t == CLASSIFICATION_PB) 		setting.TimeResult_write_int(shellquery, SkNNm, result_num, "CLASSIFICATION/CLASSIFICATION_PB", proto); 
	else if(user_query->cquery_t == CLASSIFICATION_PGI) 	setting.TimeResult_write_int(shellquery, SkNNm, result_num, "CLASSIFICATION/CLASSIFICATION_PGI", proto);
	else if(user_query->cquery_t == CLASSIFICATION_PAI) 	setting.TimeResult_write_int(shellquery, SkNNm, result_num, "CLASSIFICATION/CLASSIFICATION_PAI", proto);

	proto.protocol_free();
	return 0;
}

int kmeans_main(parsed_query* user_query)
{
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;
	
	paillier_keygen(user_query->key_s, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub, prv, user_query->key_s);


	cout<<"user_query->key_s : " <<user_query->key_s<<endl;

	int NumNode = 0;
	boundary* node = 0;

	char kd_filename[128]; 

	if(	user_query->tquery_t != KMEANS_B )
	{
		char kd_filename[128]; 
		sprintf(kd_filename, "input/Grid/Grid_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
		node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);	
	}

	char inputFilename[128]; 
	sprintf(inputFilename, "input/kmeans/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);


	paillier_ciphertext_t*** Cluster;
	int result_num=0;

	if(user_query->mquery_t == KMEANS_B)			Cluster = proto.Clustering_m(cipher, user_query->data_n, user_query->k, 2);
	else if(user_query->mquery_t == KMEANS_I) 		Cluster = proto.Clustering_Grid(cipher, node, NumNode, user_query->data_n, user_query->k);	
	else if(user_query->mquery_t == KMEANS_GI) 		Cluster = proto.Clustering_Grid_preprocessing(cipher, node, NumNode, user_query->data_n, user_query->k);
//	else if(user_query->mquery_t == KMEANS_PB) 		cluseter = proto.sRange_PM(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); 
//	else if(user_query->mquery_t == KMEANS_PGI) 	cluseter = proto.sRange_PGI(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); 
//	else if(user_query->mquery_t == KMEANS_PAI) 	cluseter = proto.sRange_PAI(ciper, q, node, user_query->data_n, user_query->NumNode, &result_num); 


	printf("result\n");
	int** cluseter = (int**)malloc(sizeof(int*)*user_query->k);
	for(int i = 0 ; i < user_query->k ; i++ ){
		cluseter[i] = (int*)malloc(sizeof(int)*user_query->dim_n);
		for(int j = 0 ; j < user_query->dim_n ; j++ ){
			cluseter[i][j] = mpz_get_ui(paillier_dec(0, pub, prv, Cluster[i][j])->m);
			printf("%d ", cluseter[i][j]);
		}
		printf("\n");
	}

	if(user_query->mquery_t == KMEANS_B)			setting.TimeResult_write_int(shellquery, cluseter, result_num, "KMEANS/KMEANS_B", proto); 
	else if(user_query->mquery_t == KMEANS_I) 		setting.TimeResult_write_int(shellquery, cluseter, result_num, "KMEANS/KMEANS_I", proto); 
	else if(user_query->mquery_t == KMEANS_GI) 		setting.TimeResult_write_int(shellquery, cluseter, result_num, "KMEANS/KMEANS_GI", proto);
	else if(user_query->mquery_t == KMEANS_PB) 		setting.TimeResult_write_int(shellquery, cluseter, result_num, "KMEANS/KMEANS_PB", proto); 
	else if(user_query->mquery_t == KMEANS_PGI) 	setting.TimeResult_write_int(shellquery, cluseter, result_num, "KMEANS/KMEANS_PGI", proto);
	else if(user_query->mquery_t == KMEANS_PAI) 	setting.TimeResult_write_int(shellquery, cluseter, result_num, "KMEANS/KMEANS_PAI", proto);

	proto.protocol_free();
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
