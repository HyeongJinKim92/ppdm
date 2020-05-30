#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>


#include "ppdm.h"

using namespace std;


char shellquery[256]; 		
protocol proto;
setting setting;


void ppdm_main(char* argv[])
{
	parsed_query * user_query = parsing_query(argv);
	print_query(user_query);
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

	boundary* node = 0;
	
	char kd_filename[128]; 		
	sprintf(kd_filename, "input/kd_range/KD_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
	node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);	

	char inputFilename[128]; 
	sprintf(inputFilename, "input/range/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);
	
	int result_num=0;
	int** result = 0;
	

	std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    std::chrono::duration<float> duration_sec;


    startTime = std::chrono::system_clock::now();

	if(user_query->rquery_t == RANGE_B)			result = proto.sRange_B(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangeB
	else if(user_query->rquery_t == RANGE_I) 	result = proto.sRange_I(cipher, q, node, user_query->data_n, NumNode, &result_num);  //rangeI
	else if(user_query->rquery_t == RANGE_GI) 	result = proto.sRange_G(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangeGI
	else if(user_query->rquery_t == RANGE_PB) 	result = proto.sRange_PB(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangePB
	else if(user_query->rquery_t == RANGE_PGI) 	result = proto.sRange_PGI(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangePGI
	else if(user_query->rquery_t == RANGE_PAI) 	result = proto.sRange_PAI(cipher, q, node, user_query->data_n, NumNode, &result_num); //rangePAI


	endTime = std::chrono::system_clock::now();
	duration_sec = endTime - startTime;
	proto.total_time = duration_sec.count();
	std::cout << "RANGE TIME : " << proto.total_time << " sec" << std::endl; 


	if(user_query->rquery_t == RANGE_B)			setting.TimeResult_write_map(shellquery, result, result_num, "RANGE/RANGE_B", proto); //rangeB
	else if(user_query->rquery_t == RANGE_I) 	setting.TimeResult_write_map(shellquery, result, result_num, "RANGE/RANGE_I", proto);  //rangeI
	else if(user_query->rquery_t == RANGE_GI) 	setting.TimeResult_write_map(shellquery, result, result_num, "RANGE/RANGE_GI", proto); //rangeGI
	else if(user_query->rquery_t == RANGE_PB) 	setting.TimeResult_write_map(shellquery, result, result_num, "RANGE/RANGE_PB", proto); //rangePB
	else if(user_query->rquery_t == RANGE_PGI) 	setting.TimeResult_write_map(shellquery, result, result_num, "RANGE/RANGE_PGI", proto); //rangePAI
	else if(user_query->rquery_t == RANGE_PAI) 	setting.TimeResult_write_map(shellquery, result, result_num, "RANGE/RANGE_PAI", proto); //rangePAI

	for( int i = 0 ; i < result_num ; i++ ){
		printf("%d result : ", (i+1) );
		for ( int j = 0 ; j < user_query->dim_n ; j++ ){
			printf("%d \t", result[i][j]);
		}
		printf("\n");
	}	

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
	/*
	cout << "G_CMP : " <<proto.G_CMP(10, 0, 10, 0) << endl;
	cout << "G_CMP : " <<proto.G_CMP(10, 3, 10, 3) << endl;
	cout << "G_CMP : " <<proto.G_CMP(10, 3, 10, 4) << endl;
	cout << "G_CMP : " <<proto.G_CMP(10, 4, 10, 3) << endl;
	cout << "G_CMP : " <<proto.G_CMP(10, 99, 100, 3) << endl;

	cout << "G_CMP : " <<proto.G_CMP(1, 1, 1, 0) << endl;
	cout << "G_CMP : " <<proto.G_CMP(1, 0, 1, 1) << endl;

	cout << "G_CMP : " <<proto.G_CMP(1, 0, 2, 0) << endl;
	cout << "G_CMP : " <<proto.G_CMP(1, 1, 2, 0) << endl;
	cout << "G_CMP : " <<proto.G_CMP(1, 0, 2, 1) << endl;


	cout << "G_CMP : " <<proto.G_CMP(2, 0, 1, 0) << endl;
	cout << "G_CMP : " <<proto.G_CMP(2, 1, 1, 0) << endl;
	cout << "G_CMP : " <<proto.G_CMP(2, 0, 1, 1) << endl;
	*/

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
	boundary* node = 0;

	char kd_filename[128]; 
	sprintf(kd_filename, "input/kd_topk/KD_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
	node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);	

	char inputFilename[128]; 
	sprintf(inputFilename, "input/topk/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);

	int result_num=0;
	int** topk = 0;

	
	std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    std::chrono::duration<float> duration_sec;


    startTime = std::chrono::system_clock::now();


	if(user_query->tquery_t == TOPK_B)			topk = proto.STopk_B(cipher, cipher_query, user_query->data_n); //TOPKB
	else if(user_query->tquery_t == TOPK_I) 	topk =  proto.STopk_I(cipher, cipher_query, node, max_val, user_query->data_n, NumNode);  //TOPKI
	else if(user_query->tquery_t == TOPK_GI) 	topk = proto.STopk_G(cipher, cipher_query, node, max_val, user_query->data_n, NumNode); //TOPKGI
	else if(user_query->tquery_t == TOPK_PB) 	topk = proto.STopk_PB(cipher, cipher_query, node, max_val, user_query->data_n, NumNode); //TOPKPGI
	else if(user_query->tquery_t == TOPK_PGI) 	topk = proto.STopk_PGI(cipher, cipher_query, node, max_val, user_query->data_n, NumNode); //TOPKPGI
	else if(user_query->tquery_t == TOPK_PAI) 	topk = proto.STopk_PAI(cipher, cipher_query, node, max_val, user_query->data_n, NumNode); //TOPKPAI


	endTime = std::chrono::system_clock::now();
	duration_sec = endTime - startTime;
	proto.total_time = duration_sec.count();
	std::cout << "topk TIME : " << proto.total_time << " sec" << std::endl; 



	if(user_query->tquery_t == TOPK_B)			setting.TimeResult_write_map(shellquery, topk, user_query->k, "TOPK/TOPK_B", proto); //TOPKB
	else if(user_query->tquery_t == TOPK_I) 	setting.TimeResult_write_map(shellquery, topk, user_query->k, "TOPK/TOPK_I", proto);  //TOPKI
	else if(user_query->tquery_t == TOPK_GI) 	setting.TimeResult_write_map(shellquery, topk, user_query->k, "TOPK/TOPK_GI", proto); //TOPKGI
	else if(user_query->tquery_t == TOPK_PB) 	setting.TimeResult_write_map(shellquery, topk, user_query->k, "TOPK/TOPK_PB", proto); //TOPKPGI
	else if(user_query->tquery_t == TOPK_PGI) 	setting.TimeResult_write_map(shellquery, topk, user_query->k, "TOPK/TOPK_PGI", proto); //TOPKPAI
	else if(user_query->tquery_t == TOPK_PAI) 	setting.TimeResult_write_map(shellquery, topk, user_query->k, "TOPK/TOPK_PAI", proto); //TOPKPAI

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
	boundary* node = 0;


	char kd_filename[128]; 
	sprintf(kd_filename, "input/kd_knn/KD_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
	node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);


	char inputFilename[128]; 
	sprintf(inputFilename, "input/knn/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);

	int result_num=0;
	int** SkNNm = 0;


	std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    std::chrono::duration<float> duration_sec;


    startTime = std::chrono::system_clock::now();

	if(user_query->kquery_t == KNN_B)			SkNNm = proto.SkNN_B(cipher, cipher_query, user_query->k, user_query->data_n); //KNNB
	else if(user_query->kquery_t == KNN_I) 		SkNNm = proto.SkNN_I(cipher, cipher_query, node, user_query->k, user_query->data_n, NumNode);  //KNNI
	else if(user_query->kquery_t == KNN_GI) 	SkNNm = proto.SkNN_G(cipher, cipher_query, node, user_query->k, user_query->data_n, NumNode); //KNNGI
	else if(user_query->kquery_t == KNN_PB) 	SkNNm = proto.SkNN_PB(cipher, cipher_query, user_query->k, user_query->data_n); //KNNPB
	else if(user_query->kquery_t == KNN_PGI) 	SkNNm = proto.SkNN_PGI(cipher, cipher_query, node, user_query->k, user_query->data_n, NumNode); //KNN PGI
	else if(user_query->kquery_t == KNN_PAI) 	SkNNm = proto.SkNN_PAI(cipher, cipher_query, node, user_query->k, user_query->data_n, NumNode); //KNN PAI

	endTime = std::chrono::system_clock::now();
	duration_sec = endTime - startTime;
	proto.total_time = duration_sec.count();
	std::cout << "kNN TIME : " << proto.total_time << " sec" << std::endl; 




	if(user_query->kquery_t == KNN_B)			setting.TimeResult_write_map(shellquery, SkNNm, result_num, "KNN/KNN_B", proto); //KNNB
	else if(user_query->kquery_t == KNN_I) 		setting.TimeResult_write_map(shellquery, SkNNm, result_num, "KNN/KNN_I", proto);  //KNNI
	else if(user_query->kquery_t == KNN_GI) 	setting.TimeResult_write_map(shellquery, SkNNm, result_num, "KNN/KNN_GI", proto); //KNNGI
	else if(user_query->kquery_t == KNN_PB) 	setting.TimeResult_write_map(shellquery, SkNNm, result_num, "KNN/KNN_PB", proto); //KNNPGI
	else if(user_query->kquery_t == KNN_PGI) 	setting.TimeResult_write_map(shellquery, SkNNm, result_num, "KNN/KNN_PGI", proto); //KNNPAI
	else if(user_query->kquery_t == KNN_PAI) 	setting.TimeResult_write_map(shellquery, SkNNm, result_num, "KNN/KNN_PAI", proto); //KNNPAI

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


	char kd_filename[128]; 
	sprintf(kd_filename, "input/kd_classification/KD_d%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
	node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);

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
	int** classification = 0;



	std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    std::chrono::duration<float> duration_sec;


    startTime = std::chrono::system_clock::now();



	if(user_query->cquery_t == CLASSIFICATION_B)			classification = proto.Classification_B(cipher, cipher_query, Entire_set, user_query->k, user_query->data_n, Entire_num);	
	else if(user_query->cquery_t == CLASSIFICATION_I) 		classification = proto.Classification_I(cipher, cipher_query, Entire_set, node, user_query->k, user_query->data_n, NumNode, Entire_num);	
	else if(user_query->cquery_t == CLASSIFICATION_GI) 		classification = proto.Classification_G(cipher, cipher_query, Entire_set, node, user_query->k, user_query->data_n, NumNode, Entire_num);
	else if(user_query->cquery_t == CLASSIFICATION_PB) 		classification = proto.Classification_PB(cipher, cipher_query, Entire_set, node, user_query->k, user_query->data_n, NumNode, Entire_num);
	else if(user_query->cquery_t == CLASSIFICATION_PGI) 	classification = proto.Classification_PGI(cipher, cipher_query, Entire_set, node, user_query->k, user_query->data_n, NumNode, Entire_num);
	else if(user_query->cquery_t == CLASSIFICATION_PAI) 	classification = proto.Classification_PAI(cipher, cipher_query, Entire_set, node, user_query->k, user_query->data_n, NumNode, Entire_num);


	endTime = std::chrono::system_clock::now();
	duration_sec = endTime - startTime;
	proto.total_time = duration_sec.count();
	std::cout << "RANGE TIME : " << proto.total_time << " sec" << std::endl; 



	if(user_query->cquery_t == CLASSIFICATION_B)			setting.TimeResult_write_map(shellquery, classification, result_num, "CLASSIFICATION/CLASSIFICATION_B", proto); 
	else if(user_query->cquery_t == CLASSIFICATION_I) 		setting.TimeResult_write_map(shellquery, classification, result_num, "CLASSIFICATION/CLASSIFICATION_I", proto); 
	else if(user_query->cquery_t == CLASSIFICATION_GI) 		setting.TimeResult_write_map(shellquery, classification, result_num, "CLASSIFICATION/CLASSIFICATION_GI", proto);
	else if(user_query->cquery_t == CLASSIFICATION_PB) 		setting.TimeResult_write_map(shellquery, classification, result_num, "CLASSIFICATION/CLASSIFICATION_PB", proto); 
	else if(user_query->cquery_t == CLASSIFICATION_PGI) 	setting.TimeResult_write_map(shellquery, classification, result_num, "CLASSIFICATION/CLASSIFICATION_PGI", proto);
	else if(user_query->cquery_t == CLASSIFICATION_PAI) 	setting.TimeResult_write_map(shellquery, classification, result_num, "CLASSIFICATION/CLASSIFICATION_PAI", proto);

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
	sprintf(kd_filename, "input/Grid/Grid_%d_m%d_L%d_h%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s, user_query->tree_level);
	node = setting.KdInfo_read(kd_filename, user_query->dim_n, &NumNode);	


	char inputFilename[128]; 
	sprintf(inputFilename, "input/kmeans/d%d_m%d_L%d.txt", user_query->data_n, user_query->dim_n, user_query->bit_s);
	paillier_ciphertext_t*** cipher = setting.InputData_read(inputFilename, user_query->dim_n, user_query->data_n);


	paillier_ciphertext_t*** Cluster;



	std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    std::chrono::duration<float> duration_sec;

    startTime = std::chrono::system_clock::now();

	if(user_query->mquery_t == KMEANS_B)			Cluster = proto.Clustering_B(cipher, user_query->data_n, user_query->k, 2);
	else if(user_query->mquery_t == KMEANS_I) 		Cluster = proto.Clustering_Grid(cipher, node, NumNode, user_query->data_n, user_query->k);	
	else if(user_query->mquery_t == KMEANS_GI) 		Cluster = proto.Clustering_Grid_preprocessing(cipher, node, NumNode, user_query->data_n, user_query->k);
	else if(user_query->mquery_t == KMEANS_PB) 		Cluster = proto.Clustering_PB(cipher, node, NumNode, user_query->data_n, user_query->k);
	else if(user_query->mquery_t == KMEANS_PGI) 	Cluster = proto.Clustering_PGI(cipher, node, NumNode, user_query->data_n, user_query->k);
	else if(user_query->mquery_t == KMEANS_PAI) 	Cluster = proto.Clustering_PAI(cipher, node, NumNode, user_query->data_n, user_query->k);

	endTime = std::chrono::system_clock::now();
	duration_sec = endTime - startTime;
	proto.total_time = duration_sec.count();
	std::cout << "KMEANS TIME : " << proto.total_time << " sec" << std::endl; 



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

	if(user_query->mquery_t == KMEANS_B)			setting.TimeResult_write_map(shellquery, cluseter, user_query->k, "KMEANS/KMEANS_B", proto); 
	else if(user_query->mquery_t == KMEANS_I) 		setting.TimeResult_write_map(shellquery, cluseter, user_query->k, "KMEANS/KMEANS_I", proto); 
	else if(user_query->mquery_t == KMEANS_GI) 		setting.TimeResult_write_map(shellquery, cluseter, user_query->k, "KMEANS/KMEANS_GI", proto);
	else if(user_query->mquery_t == KMEANS_PB) 		setting.TimeResult_write_map(shellquery, cluseter, user_query->k, "KMEANS/KMEANS_PB", proto); 
	else if(user_query->mquery_t == KMEANS_PGI) 	setting.TimeResult_write_map(shellquery, cluseter, user_query->k, "KMEANS/KMEANS_PGI", proto);
	else if(user_query->mquery_t == KMEANS_PAI) 	setting.TimeResult_write_map(shellquery, cluseter, user_query->k, "KMEANS/KMEANS_PAI", proto);

	proto.protocol_free();
	return 0;
}


/*
	parse the user query 
*/
parsed_query* parsing_query(char* argv[])
{
	parsed_query* query = new parsed_query;
	//query setting
	query->query = (QUERY)atoi(argv[1]);
	
	//query method setting
	if(query->query == RANGE)				query->rquery_t = (RANGE_TYPE)atoi(argv[2]);
	else if(query->query == TOPK)			query->tquery_t = (TOPK_TYPE)atoi(argv[2]);
	else if(query->query == KNN)			query->kquery_t = (KNN_TYPE)atoi(argv[2]);
	else if(query->query == CLASSIFICATION)	query->cquery_t = (CLASSIFICATION_TYPE)atoi(argv[2]);
	else if(query->query == KMEANS)			query->mquery_t = (KMEANS_TYPE)atoi(argv[2]);
	else if(query->query == ASSOCATION)		{}
	else if(query->query == TEST)			{}

	//parameter setting
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
	
	//query merge
	if(query->query == RANGE)
	{	
		switch(query->rquery_t)
		{
			case RANGE_B:
				strcat(shellquery, "RANGE_B");
				break;
			case RANGE_I:
				strcat(shellquery, "RANGE_I");
				break;
			case RANGE_GI:
				strcat(shellquery, "RANGE_GI");
				break;
			case RANGE_PB:
				strcat(shellquery, "RANGE_PB");
				break;
			case RANGE_PGI:
				strcat(shellquery, "RANGE_PGI");
				break;
			case RANGE_PAI:
				strcat(shellquery, "RANGE_PAI");
				break;
			default:
				break;
		}
	}
	else if(query->query == TOPK)
	{
		switch(query->tquery_t)
		{
			case TOPK_B:
				strcat(shellquery, "TOPK_B");
				break;
			case TOPK_I:
				strcat(shellquery, "TOPK_I");
				break;
			case TOPK_GI:
				strcat(shellquery, "TOPK_GI");
				break;
			case TOPK_PB:
				strcat(shellquery, "TOPK_PB");
				break;
			case TOPK_PGI:
				strcat(shellquery, "TOPK_PGI");
				break;
			case TOPK_PAI:
				strcat(shellquery, "TOPK_PAI");
				break;
			default:
				break;
		}
	}
	else if(query->query == KNN)
	{
		switch(query->kquery_t)
		{
			case KNN_B:
				strcat(shellquery, "KNN_B");
				break;
			case KNN_I:
				strcat(shellquery, "KNN_I");
				break;
			case KNN_GI:
				strcat(shellquery, "KNN_GI");
				break;
			case KNN_PB:
				strcat(shellquery, "KNN_PB");
				break;
			case KNN_PGI:
				strcat(shellquery, "KNN_PGI");
				break;
			case KNN_PAI:
				strcat(shellquery, "KNN_PAI");
				break;
			default:
				break;
		}
	}
	else if(query->query == CLASSIFICATION)
	{
		switch(query->cquery_t)
		{
			case CLASSIFICATION_B:
				strcat(shellquery, "CLASSIFICATION_B");
				break;
			case CLASSIFICATION_I:
				strcat(shellquery, "CLASSIFICATION_I");
				break;
			case CLASSIFICATION_GI:
				strcat(shellquery, "CLASSIFICATION_GI");
				break;
			case CLASSIFICATION_PB:
				strcat(shellquery, "CLASSIFICATION_PB");
				break;
			case CLASSIFICATION_PGI:
				strcat(shellquery, "CLASSIFICATION_PGI");
				break;
			case CLASSIFICATION_PAI:
				strcat(shellquery, "CLASSIFICATION_PAI");
				break;
			default:
				break;
		}
	}
	else if(query->query == KMEANS)			
	{
		switch(query->mquery_t)
		{
			case KMEANS_B:
				strcat(shellquery, "KMEANS_B");
				break;
			case KMEANS_I:
				strcat(shellquery, "KMEANS_I");
				break;
			case KMEANS_GI:
				strcat(shellquery, "KMEANS_GI");
				break;
			case KMEANS_PB:
				strcat(shellquery, "KMEANS_PB");
				break;
			case KMEANS_PGI:
				strcat(shellquery, "KMEANS_PGI");
				break;
			case KMEANS_PAI:
				strcat(shellquery, "KMEANS_PAI");
				break;
			default:
				break;
		}
	}

	strcat(shellquery, " TREE_LEVEL : ");
	strcat(shellquery, argv[3]);
	strcat(shellquery, " THREAD : ");
	strcat(shellquery, argv[4]);
	strcat(shellquery, " DATA : ");
	strcat(shellquery, argv[5]);
	strcat(shellquery, " DIMENSION : ");
	strcat(shellquery, argv[6]);
	strcat(shellquery, " BITS : ");
	strcat(shellquery, argv[7]);
	strcat(shellquery, " k : ");
	strcat(shellquery, argv[8]);
	strcat(shellquery, " KEY : ");
	strcat(shellquery, argv[9]);


	strcat(shellquery, "  QUERY : ");
	int i = 0;
	while(argv[10+i] != NULL)
	{	
		strcat(shellquery, argv[10+i]);
		strcat(shellquery, " ");
		i++;
	}

	return query;
}


/*
	Output the query
*/
void print_query(parsed_query * query)
{
	if(query->query == RANGE)
	{
		switch(query->rquery_t){
			case RANGE_B:	
				cout << "RANGE_B" <<endl; break;
			case RANGE_I:	
				cout << "RANGE_I" <<endl; break;
			case RANGE_GI:	
				cout << "RANGE_GI" <<endl; break;
			case RANGE_PB:	
				cout << "RANGE_PB" <<endl; break;
			case RANGE_PGI:	
				cout << "RANGE_PGI" <<endl; break;
			case RANGE_PAI:	
				cout << "RANGE_PAI" <<endl; break;
			default: 
				break;
		}
	}
	else if(query->query == TOPK)
	{
		switch(query->tquery_t){
			case TOPK_B:	
				cout << "TOPK_B" <<endl; break;
			case TOPK_I:	
				cout << "TOPK_I" <<endl; break;
			case TOPK_GI:	
				cout << "TOPK_GI" <<endl; break;
			case TOPK_PB:	
				cout << "TOPK_PB" <<endl; break;
			case TOPK_PGI:	
				cout << "TOPK_PGI" <<endl; break;
			case TOPK_PAI:	
				cout << "TOPK_PAI" <<endl; break;
			default: 
				break;
		}
	}
	else if(query->query == KNN)			
	{
		switch(query->kquery_t){
			case KNN_B:	
				cout << "KNN_B" <<endl; break;
			case KNN_I:	
				cout << "KNN_I" <<endl; break;
			case KNN_GI:	
				cout << "KNN_GI" <<endl; break;
			case KNN_PB:	
				cout << "KNN_PB" <<endl; break;
			case KNN_PGI:	
				cout << "KNN_PGI" <<endl; break;
			case KNN_PAI:	
				cout << "KNN_PAI" <<endl; break;
			default: 
				break;
		}
	}
	else if(query->query == CLASSIFICATION)	
	{
		switch(query->cquery_t){
			case CLASSIFICATION_B:	
				cout << "CLASSIFICATION_B" <<endl; break;
			case CLASSIFICATION_I:	
				cout << "CLASSIFICATION_I" <<endl; break;
			case CLASSIFICATION_GI:	
				cout << "CLASSIFICATION_GI" <<endl; break;
			case CLASSIFICATION_PB:	
				cout << "CLASSIFICATION_PB" <<endl; break;
			case CLASSIFICATION_PGI:	
				cout << "CLASSIFICATION_PGI" <<endl; break;
			case CLASSIFICATION_PAI:	
				cout << "CLASSIFICATION_PAI" <<endl; break;
			default: 
				break;
		}
	}	
	else if(query->query == KMEANS)			
	{
		switch(query->mquery_t){
			case KMEANS_B:	
				cout << "KMEANS_B" <<endl; break;
			case KMEANS_I:	
				cout << "KMEANS_I" <<endl; break;
			case KMEANS_GI:	
				cout << "KMEANS_GI" <<endl; break;
			case KMEANS_PB:	
				cout << "KMEANS_PB" <<endl; break;
			case KMEANS_PGI:	
				cout << "KMEANS_PGI" <<endl; break;
			case KMEANS_PAI:	
				cout << "KMEANS_PAI" <<endl; break;
			default: 
				break;
		}
	}
	else if(query->query == ASSOCATION)		{}
	else if(query->query == TEST)			{}

	cout << "TREE : " << query->tree_level << "\t THREAD : " << query->thread_n << "\t DATA : " << query->data_n << "\t DIM : " << query->dim_n << "\t BITS : " << query->bit_s << "\t KEY : " << query->key_s << endl;
}
