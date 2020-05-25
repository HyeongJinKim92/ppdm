#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <sstream>
#include <fstream>

#include "ppdm.h"
/*
protocol proto;
setting setting;

int* Q = 0;
int* Q2 =0;  //for range
char shellquery[128]; 		


int SRange_main(int datanum, int dimension, int requiredK ,int bitsize, int service_type, int tree_level, int modul) {
	printf("Start SRange(0~3)\n");
	int modulus = modul;
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap;


	paillier_keygen(modulus, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub,prv, modulus);

	int dim = dimension;
	int k = requiredK;
	int NumData = datanum;
	int NumNode=0;

	boundary q;
	
		
	// query setting start
	q.LL =  (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	q.RR =  (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);

	// query setting end

	cout<<"modulus : " <<modulus<<endl;

	for(int i=0; i<dim; i++){
		q.LL[i] = paillier_create_enc(Q[i]);
		paillier_print("LL : ", q.LL[i]);
		q.RR[i] = paillier_create_enc(Q2[i]);
		paillier_print("RR : ", q.RR[i]);
	}

	char kd_filename[128]; 		
	boundary* node = 0;
		
	if(service_type == 0 || service_type == 1 || service_type == 2 ) {
		sprintf(kd_filename, "input/KD_d%d_m%d_L%d_h%d.txt",datanum,dimension,bitsize,tree_level);
		//printf("kd_filename : %s\n",kd_filename);
		node = setting.KdInfo_read(kd_filename,dim,&NumNode);	
		//printf("NumNode : %d\n",NumNode);
	}
	
	

	int result_num=0;
	char inputFilename[128]; 
	sprintf(inputFilename, "input/d%d_m%d_L%d.txt",datanum,dimension,bitsize);
	//printf("inputfilename : %s\n",inputFilename);

	paillier_ciphertext_t*** ciper = setting.InputData_read(inputFilename,dim,NumData);
	
	printf("\n===== sRange start =====\n");
	startTime = clock();
	int** result = 0;
	
	if(service_type == 0){
		printf("Service_type 0 sRange_M");
		result =proto.sRange_M(ciper, q, node, NumData, NumNode, &result_num);
	}
	else if(service_type == 1) {
		printf("Service_type 1 sRange_I");
		result =proto.sRange_I(ciper, q, node, NumData, NumNode,&result_num);  //rangeI
	}
	else if(service_type == 2) {
		printf("Service_type 2 sRange_G");
		result = proto.sRange_G(ciper, q, node, NumData, NumNode,&result_num); //rangeG
	}
	
	endTime = clock();
	proto.total_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("wall time : %f\n", proto.total_time);

	if(service_type==0){
		setting.TimeResult_write_int_forRange(gap,result,k,result_num,NumData,"sRange_M",bitsize,dim,tree_level,modulus,proto);
	}
	else if(service_type == 1) {
		setting.TimeResult_write_int_forRange(gap,result,k,result_num,NumData,"sRange_I",bitsize,dim,tree_level,modulus,proto);
	}
	else if(service_type == 2) {
		setting.TimeResult_write_int_forRange(gap,result,k,result_num,NumData,"sRange_G",bitsize,dim,tree_level,modulus,proto);
	}

	for( i = 0 ; i < result_num ; i++ ){
		printf("%d result : ", (i+1) );
		for ( j = 0 ; j < dim ; j++ ){
			printf("%d \t", result[i][j]);
		}
		printf("\n");
	}
	
	//free
	for(i=0; i<2; i++) {
		for(j=0; j<dim; j++) {
			paillier_freeciphertext(node[i].LL[j]);	
			paillier_freeciphertext(node[i].RR[j]);
		}
	}
	
	for(j=0; j<dim; j++) {
			paillier_freeciphertext(q.LL[j]);		
			paillier_freeciphertext(q.RR[j]);
	}

	free(node);
	proto.protocol_free();

	return 0;
}

int skNN_main(int datanum,int dimension,int requiredK,int bitsize, int service_type, int tree_level, int modul) {
	printf("Start skNN(10~14)\n");
	//int modulus = 1024;
	int modulus = modul;
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap;

	paillier_keygen(modulus, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub, prv, modulus);

	int dim = dimension, k = requiredK;
	int num=0;
	int NumData = datanum;

	paillier_ciphertext_t** ciper_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim); 
	for( i = 0 ; i < dim ; i++ ){
		ciper_query[i] = paillier_create_enc(Q[i]);
	}

	cout<<"modulus : " <<modulus<<endl;

	printf("\nquery : ");
	for(i = 0 ; i < dim ; i++ ){				
		gmp_printf("%Zd \t", paillier_dec(0, pub, prv, ciper_query[i]));
	}
	printf("\n\n");

	int NumNode = 0;
	char kd_filename[128]; 
	boundary* node = 0;

	if(service_type == 10 || service_type == 11 || service_type == 12 || service_type == 13 || service_type == 14) {
		sprintf(kd_filename, "input/KD_d%d_m%d_L%d_h%d.txt",datanum,dim,bitsize,tree_level);
		//printf("kd_filename : %s\n",kd_filename);
		node =setting.KdInfo_read(kd_filename,dim,&NumNode);
		//printf("NumNode : %d\n",NumNode);
	}

	int result_num=0;
	char inputFilename[128]; 
	sprintf(inputFilename, "input/d%d_m%d_L%d.txt",datanum,dim,bitsize);
	//printf("inputfilename : %s\n",inputFilename);

	paillier_ciphertext_t*** ciper = setting.InputData_read(inputFilename,dim,NumData);

	printf("\n===== SKNN start =====\n");
	startTime = clock();
	int** SkNNm = 0;

	if(service_type == 10) {
		printf("Service_type 10 SkNN_m");
		SkNNm = proto.SkNN_M(ciper, ciper_query, k, NumData);
	}
	else if(service_type == 11) {
		printf("Service_type 11 SkNN_I");
		SkNNm = proto.SkNN_I(ciper, ciper_query, node, k, NumData, NumNode);
	}
	else if(service_type == 12){ 
		//skNNG
		printf("Service_type 12 SkNN_G");
		SkNNm = proto.SkNN_G(ciper, ciper_query, node, k, NumData, NumNode); // G로 변경해야함!
		//SkNNm = proto.GSRO_SkNNb(ciper, ciper_query, node, k, NumData, NumNode); // G로 변경해야함!
	}	
	else if(service_type == 13) {
		//skNNb
		printf("Service_type 13 SkNN_b");
		SkNNm = proto.SkNN_b(ciper, ciper_query, k, NumData); 
		//SkNNm = proto.SSED_test(ciper, ciper_query, node, k, NumData, NumNode);
	}
	else if(service_type == 14) {
		//SkNNGI_DP_SSED
		printf("Service_type 14 DP_SkNNsi");
		SkNNm = proto.DP_SkNN_I(ciper, ciper_query, node, k, NumData, NumNode);

	}
	endTime = clock();

	proto.total_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("\ntime : %f\n", proto.total_time);

	if(service_type == 10) {
		setting.TimeResult_write_int(shellquery, gap,SkNNm,k,NumData,"SkNN_m",bitsize,dim,tree_level,modulus,proto);
	}
	else if(service_type == 11) {
		setting.TimeResult_write_int(shellquery, gap,SkNNm,k,NumData,"SkNN_I",bitsize,dim,tree_level,modulus,proto);
	}
	else if(service_type == 12) {
		setting.TimeResult_write_int(shellquery, gap,SkNNm,k,NumData,"SkNN_G",bitsize,dim,tree_level,modulus,proto);   
	}
	else if(service_type == 13) {
		setting.TimeResult_write_int(shellquery, gap,SkNNm,k,NumData,"SkNN_b",bitsize,dim,tree_level,modulus,proto);   
	}
	else if(service_type == 14) {
		setting.TimeResult_write_int(shellquery, gap,SkNNm,k,NumData,"DP_SkNN_I",bitsize,dim,tree_level,modulus,proto);		
	}

	for( i = 0 ; i < k ; i++ ){
		printf("%d result : ", (i+1) );
		for ( j = 0 ; j < dim ; j++ ){
			printf("%d \t", SkNNm[i][j]);
		}
		printf("\n");
	}

	proto.protocol_free();
	return 0;
}

//int sTopk_main(int datanum,int dimension, int requiredK) {
int sTopk_main(int datanum,int dimension,int requiredK,int bitsize, int service_type, int tree_level, int modul){
	printf("Start sTopk(20~23)\n");
	int modulus = modul;
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap;

	paillier_keygen(modulus, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub, prv, modulus);

	int dim = dimension, k = requiredK;
	int num=0;	
	int NumData = datanum;
	
	paillier_ciphertext_t** ciper_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim+1); 

	paillier_ciphertext_t** max_val = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim); 

	paillier_ciphertext_t* ciper_zero = paillier_create_enc(0);

	int MAX_VAL = 0;
	MAX_VAL = sqrt(pow(2, bitsize)/dim);
	MAX_VAL = 512;
	for( i = 0 ; i < dim ; i++ ){
		if(Q[i] <0) {
			ciper_query[i] = paillier_create_enc(Q[i]*(-1));	// 양수로 바꿔서 일단 암호화 한 후, 
			paillier_subtract(pub, ciper_query[i] , ciper_zero , ciper_query[i]);	 // 0에서 빼서 음수로 전환
		}
		else {
			ciper_query[i] = paillier_create_enc(Q[i]);
		}
		max_val[i] = paillier_create_enc(MAX_VAL);
		gmp_printf("max val : %Zd \n", paillier_dec(0, pub, prv, max_val[i]));
	}
	

	int tmp=abs(Q[0]);
	for(int i=1;i<dim;i++){
		if(tmp<abs(Q[i]))
			tmp=abs(Q[i]);
	}

	ciper_query[dim] = paillier_create_enc(tmp); 

	//gmp_printf("hint %Zd \n", paillier_dec(0, pub, prv, ciper_query[dim]));

	cout<<"modulus : " <<modulus<<endl;

	printf("query : ");

	for(i = 0 ; i < dim+1 ; i++ ){				
		if(i == dim)	// hint
			printf("hint : ");
			gmp_printf("%Zd \t", paillier_dec(0, pub, prv, ciper_query[i]));
		cout<<endl;
	}
	printf("\n");

	int NumNode = 0;
	char kd_filename[128]; 
	boundary* node = 0;

	if(service_type == 20 || service_type == 21 || service_type == 22 || service_type == 23) {
		sprintf(kd_filename, "input/topk/KD_d%d_m%d_L%d_h%d.txt",datanum,dim,bitsize,tree_level);
		//printf("kd_filename : %s\n",kd_filename);
		node = setting.KdInfo_read(kd_filename,dim,&NumNode);
		//printf("NumNode : %d\n",NumNode);
	}

	int result_num=0;
	char inputFilename[128]; 
	sprintf(inputFilename, "input/topk/d%d_m%d_L%d.txt",datanum,dim,bitsize);
	printf("inputfilename : %s\n",inputFilename);

	paillier_ciphertext_t*** ciper = setting.InputData_read(inputFilename,dim,NumData);
	printf("\n===== Top-k start =====\n");
	startTime = clock();
	int** topk = 0;

	if(service_type == 20) {
		printf("Service_type 20 STopk_m\n");
		topk = proto.STopk_M(ciper, ciper_query, NumData);
	}else if(service_type == 21) {
		printf("Service_type 21 STopk_I\n");
		topk =  proto.STopk_I(ciper, ciper_query, node, max_val, NumData, NumNode);
	}else if(service_type == 22) {
		printf("Service_type 22 STopk_G\n");
		topk = proto.STopk_G(ciper, ciper_query, node, max_val, NumData, NumNode);
	}

	
	endTime = clock();

	proto.total_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("\ntime : %f\n", proto.total_time);
	
	if(service_type == 20) {
		setting.TimeResult_write_int(shellquery, gap,topk,k,NumData,"STopk_m",bitsize,dim,tree_level,modulus,proto);
	}
	else if(service_type == 21) {
		setting.TimeResult_write_int(shellquery, gap,topk,k,NumData,"STopk_I",bitsize,dim,tree_level,modulus,proto);
	}
	else if(service_type == 22) {
		setting.TimeResult_write_int(shellquery, gap,topk,k,NumData,"STopk_G",bitsize,dim,tree_level,modulus,proto);
	}

	for( i = 0 ; i < k ; i++ ){
		printf("%d result : ", (i+1));
		for ( j = 0 ; j < dim ; j++ ){
			printf("%d \t", topk[i][j]);
		}
		printf("\n");
	}

	//proto.protocol_free();
	return 0;
}
int Classification_main(int datanum, int dimension, int requiredK, int bitsize, int service_type, int tree_level, int modul) {
	printf("Start Classification(30~34)\n");
	// Classification 전체 분류 체계를 전달해줄 인자값이 필요하다.
	//int modulus = 1024;
	int modulus = modul;
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;

	
	time_t startTime = 0;
	time_t endTime = 0;
	float gap;
	
	paillier_keygen(modulus, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub, prv, modulus);

	int dim = dimension , k = requiredK;
	int num = 0;
	int NumData = datanum;

	paillier_ciphertext_t** ciper_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim); 
	for( i = 0 ; i < dim ; i++ ){
		ciper_query[i] = paillier_create_enc(Q[i]);
	}

	cout<<"modulus : " <<modulus<<endl;

	printf("\nquery : ");
	for(i = 0 ; i < dim ; i++ ){				
		gmp_printf("%Zd \t", paillier_dec(0, pub, prv, ciper_query[i]));
	}
	printf("\n\n");

	int NumNode = 0;
	char kd_filename[128]; 
	boundary* node = 0;

	if(service_type == 30 || service_type == 31 || service_type == 32 || service_type == 33 || service_type == 34) {
		sprintf(kd_filename, "input/Classification/KD_d%d_m%d_L%d_h%d.txt",datanum,dim,bitsize,tree_level);
		//printf("kd_filename : %s\n",kd_filename);
		node =setting.KdInfo_read(kd_filename,dim,&NumNode);
		//printf("NumNode : %d\n",NumNode);
	}

	int result_num=0;
	char inputFilename[128]; 
	sprintf(inputFilename, "input/Classification/d%d_m%d_L%d.txt",datanum,dim,bitsize);
	//printf("inputfilename : %s\n",inputFilename);

	paillier_ciphertext_t*** ciper = setting.InputData_read(inputFilename,dim+1,NumData);

	printf("\n===== Classification start =====\n");
	startTime = clock();
	int** SkNNm = 0;
////////Label 설정	
	int Entire_num = 18;
	paillier_ciphertext_t** Entire_set = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(Entire_num));
	for( i = 0 ; i < Entire_num ; i++){
		Entire_set[i] = paillier_create_enc(i);
	}


	if(service_type == 30) {
		printf("Service_type 30 Classification_m\n");
		SkNNm = proto.Classification_M(ciper, ciper_query, Entire_set, k, NumData, Entire_num);	
	}
	else if(service_type == 31) {
		paillier_ciphertext_t* R = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		printf("Service_type 31 Classification_I\n");
		SkNNm = proto.Classification_I(ciper, ciper_query, Entire_set, node, k, NumData, NumNode, Entire_num);	
		//gmp_printf("\nR : %Zd \n", paillier_dec(0, pub, prv, R));
	}
	else if(service_type == 32){ 
		printf("Service_type 32 Classification_G\n");
		SkNNm = proto.Classification_G(ciper, ciper_query, Entire_set, node, k, NumData, NumNode, Entire_num);
	}	
	else if(service_type == 33) {
		printf("Service_type 18\n");
	}
	else if(service_type == 34) {
		printf("Service_type 19\n");
	}
	endTime = clock();

	proto.total_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("\ntime : %f\n", proto.total_time);

	if(service_type == 30) {
		setting.TimeResult_write_int(shellquery, gap, SkNNm, k, NumData, "Classification_m", bitsize, dim, tree_level, modulus, proto);
	}else if(service_type == 31){
		setting.TimeResult_write_int(shellquery, gap, SkNNm, k, NumData, "Classification_I", bitsize, dim+1, tree_level, modulus, proto);
	}else if(service_type == 32){
		setting.TimeResult_write_int(shellquery, gap, SkNNm, k, NumData, "Classification_G", bitsize, dim+1, tree_level, modulus, proto);
	}

	proto.protocol_free();
	return 0;
}
int Clustering_main(int datanum, int dimension, int requiredK, int bitsize, int service_type, int tree_level, int modul) {
	printf("Start Clustering(40~42)\n");
	//int modulus = 1024;
	int modulus = modul;
	paillier_pubkey_t* pub;
	paillier_prvkey_t* prv;
	int i=0, j=0;

	
	time_t startTime = 0;
	time_t endTime = 0;
	float gap;
	
	paillier_keygen(modulus, &pub,&prv, paillier_get_rand_devurandom);

	paillier_setkey(pub, prv);
	proto.protocol_setkey(pub, prv, modulus);

	int dim = dimension , k = requiredK;
	int num = 0;
	int NumData = datanum;

	paillier_ciphertext_t** ciper_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim); 
	for( i = 0 ; i < dim ; i++ ){
		ciper_query[i] = paillier_create_enc(Q[i]);
	}

	cout<<"modulus : " <<modulus<<endl;

	printf("\nquery : ");
	for(i = 0 ; i < dim ; i++ ){				
		gmp_printf("%Zd \t", paillier_dec(0, pub, prv, ciper_query[i]));
	}
	printf("\n\n");

	int NumNode = 0;
	char kd_filename[128]; 
	boundary* node = 0;

	if(service_type == 41 || service_type == 42) {
		sprintf(kd_filename, "input/Grid/Grid_%d_m%d_L%d_h%d.txt",datanum,dim,bitsize,tree_level);
		//printf("kd_filename : %s\n",kd_filename);
		node =setting.KdInfo_read(kd_filename,dim,&NumNode);
		//printf("NumNode : %d\n",NumNode);
	}
	/*
	printf("Grid Index: %d\n", NumNode);
	for(int i = 0 ; i < NumNode ; i++){
		printf("node_id : %d ", node[i].node_id);
		for(int j = 0 ; j < dim ; j++){
			gmp_printf("%Zd ", paillier_dec(0, pub, prv, node[i].LL[j]));
		}
		for(int j = 0 ; j < dim ; j++){
			gmp_printf("%Zd ", paillier_dec(0, pub, prv, node[i].RR[j]));
		}
		printf("%d ", node[i].NumData);
		for(int j = 0 ; j < node[i].NumData ; j++){
			printf("%d ", node[i].indata_id[j]);
		}
		printf("\n");
	}
	*/
/*
	int result_num=0;
	char inputFilename[128]; 
	sprintf(inputFilename, "input/d%d_m%d_L%d.txt",datanum,dim,bitsize);
	//printf("inputfilename : %s\n",inputFilename);

	paillier_ciphertext_t*** ciper = setting.InputData_read(inputFilename,dim,NumData);
	paillier_ciphertext_t*** Cluster;
	printf("\n===== Clustering start =====\n");
	startTime = clock();
	int** cluseter = 0;
////////Label 설정	

	if(service_type == 40) {
		printf("Service_type 40 Clustering_m\n");
		Cluster = proto.Clustering_m(ciper, datanum, k, 2);
	}else if(service_type == 41) {
		printf("Service_type 41 Clustering_Grid\n");
		Cluster = proto.Clustering_Grid(ciper, node, NumNode, datanum, k);
	}else if(service_type == 42) {
		printf("Service_type 42 Clustering_Grid_preprocessing\n");
		Cluster = proto.Clustering_Grid_preprocessing(ciper, node, NumNode, datanum, k);
	}
	endTime = clock();

	proto.total_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("\ntime : %f\n", proto.total_time);
	

	printf("result\n");
	int** result = (int**)malloc(sizeof(int*)*k);
	for(int i = 0 ; i < k ; i++ ){
		result[i] = (int*)malloc(sizeof(int)*dim);
		for(int j = 0 ; j < dim ; j++ ){
			result[i][j] = mpz_get_ui(paillier_dec(0, pub, prv, Cluster[i][j])->m);
			printf("%d ", result[i][j]);
		}
		printf("\n");
	}

	if(service_type == 40) {
		setting.TimeResult_write_int(shellquery, gap, result, k, NumData, "Clustering_m", bitsize, dim, tree_level, modulus, proto);
	}else if(service_type == 41){
		setting.TimeResult_write_int(shellquery, gap, result, k, NumData, "Clustering_Grid", bitsize, dim, tree_level, modulus, proto);
	}else if(service_type == 42){
		setting.TimeResult_write_int(shellquery, gap, result, k, NumData, "Clustering_Grid_preprocessing", bitsize, dim, tree_level, modulus, proto);
	}

	proto.protocol_free();
	return 0;
}*/
int main(int argc, char* argv[]) {

	ppdm_main(argv);
/*
	int service_type = 0;	
	int datanum = 0;
	int dimension = 0;
	int bitsize = 0;
	int requiredK = 0;
	int tree_level = 0;
	int modul = 0;
	//test_main();			
	
	
	if (argc < 6 ) {
		fputs("(usage) : service_type datanum m size(l) k tree_level query_data\n", stderr);
	    fputs("service_type : \n", stderr);
		exit(1);
	}

	service_type = atoi(argv[1]); // range, knn, topk
	datanum = atoi(argv[2]);
	dimension = atoi(argv[3]);
	bitsize = atoi(argv[4]);
	requiredK = atoi(argv[5]);
	tree_level = atoi(argv[6]); // k-d tree
	modul = atoi(argv[7]);

	sprintf(shellquery, "./test %d %d %d %d %d %d %d ", service_type, datanum, dimension, bitsize, requiredK, tree_level, modul);

	cout << service_type << endl;

	Q=(int*)malloc(sizeof(int)*dimension);

	Q2=(int*)malloc(sizeof(int)*dimension);

	
	if(service_type == 0 || service_type == 1 || service_type == 2 || service_type == 3){
		// range query setting 
		for(int i = 0 ; i < dimension ; i++ ){
			Q[i] = atoi(argv[i+8]);
			Q2[i] = atoi(argv[i+8+dimension]);	
		}
		for(int i = 0 ; i < dimension ; i++ ){
			strcat(shellquery, argv[i+8]);
			strcat(shellquery, " ");
		}
		for(int i = 0 ; i < dimension ; i++ ){
			strcat(shellquery, argv[i+8+dimension]);
			strcat(shellquery, " ");
		}	

	}else{
		for(int i=0; i<dimension; i++) {
			Q[i] = atoi(argv[i+8]);
			strcat(shellquery, argv[i+8]);
			strcat(shellquery, " ");
		}
	}

	printf("\n\n");

	proto.set_Param(datanum, dimension, bitsize, requiredK, tree_level);
	
	if(service_type == 0 || service_type == 1 || service_type == 2 || service_type == 3 ) { // 0: Range_basic, 1: sRangeI, 2: sRangeG, 3: Rangeb
		SRange_main(datanum, dimension, requiredK, bitsize, service_type, tree_level, modul);
	}
	else if(service_type == 10 || service_type == 11 || service_type == 12 ||service_type == 13 || service_type == 14) {	 // 2 : skNNm, 3 : skNNI 7 : skNNI+DP
		skNN_main(datanum, dimension, requiredK, bitsize, service_type, tree_level, modul);
	}
	else if(service_type == 20 || service_type== 21 || service_type == 22 || service_type == 23) {
		//sTopk_main(datanum, dimension, requiredK);
		sTopk_main(datanum, dimension, requiredK, bitsize,  service_type,  tree_level, modul);
	}
	else if(service_type == 30 || service_type == 31 || service_type == 32 || service_type == 33 || service_type == 34) {
		Classification_main(datanum, dimension, requiredK, bitsize, service_type, tree_level, modul);
	}
	else if(service_type == 40 || service_type == 41 || service_type == 42){
		Clustering_main(datanum, dimension, requiredK, bitsize, service_type, tree_level, modul);
	}
*/
	return 0;
}


