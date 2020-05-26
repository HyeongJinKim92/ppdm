#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <iostream>
#include "protocol.h"

using namespace std;

int ** protocol::GSRO_SkNNb(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode){
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 10;
	bool verify_flag = false;

	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;
	
	int * idx = (int *)malloc(sizeof(int)*k);
	paillier_ciphertext_t* ciper_rand = paillier_create_enc(rand);
	paillier_ciphertext_t*** cand;
	paillier_ciphertext_t* cipher_Max = paillier_create_enc(10000);
	paillier_ciphertext_t* temp_coord1 ;
	paillier_ciphertext_t* temp_coord2 ;
	paillier_ciphertext_t* temp_coord3 ;
	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*3);
	paillier_ciphertext_t* temp_alpha = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp_alpha->c);
	paillier_ciphertext_t*** temp_cand ;
	paillier_ciphertext_t** ciper_nodedist = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumNode);

	paillier_ciphertext_t** shortestPoint = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	for(i=0; i<dim; i++){
		shortestPoint[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(shortestPoint[i]->c);
	}
	paillier_ciphertext_t* temp_nodedist = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp_nodedist->c);
	
	paillier_ciphertext_t* kth_dist = 0;
	
	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	for(i=0; i<NumNode; i++){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);
		ciper_nodedist[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(ciper_nodedist[i]->c);	
	}

	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t)*k);
	for( i = 0 ; i < k ; i++){
		ciper_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j++){
			ciper_result[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(ciper_result[i][j]->c);
		}
	}
	
	float progress;
	
	while(1) {
		progress = 0.1;
		startTime = clock();
		
		if(!verify_flag)
			printf("\n=== Node SRO start  ===\n");
		else
			printf("\n=== Node expansion start  ===\n");
		
		for(i=0; i<NumNode; i++) {	
			if(i/(float)NumNode >= progress){
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			if(!verify_flag)	{	//  검증 단계에서는 수행하지 않음
				for( i = 0 ; i < NumNode; i++){
					alpha[i] = DP_GSRO( q, q, node[i].LL, node[i].RR);
				}
			}

			else {	// 검증 단계에서 k번째 거리보다 가까이 존재하는 노드를 찾는 과정
				// 각 노드에서 질의와의 최단 점을 찾음
				for(j=0; j<dim; j++) {	
					psi[0] = GSCMP(q[j], node[i].LL[j]);		// 프사이 계산 (1이면 q의 좌표가 노드의 LL bound보다 크고, 0이면 q의 좌표가 작음)
					psi[1] = GSCMP(q[j], node[i].RR[j]);	// 프사이 계산 (1이면 q의 좌표가 노드의 UR bound보다 작고, 0이면 q의 좌표가 큼)
					psi[2] = SBXOR(psi[0], psi[1]);	// psi[2]가 1이면, 해당 차원에서의 질의 좌표가 노드의 LL bound 와 UR bound 사이에 속함을 의미
					temp_coord1 = SM_p1(psi[2], q[j]);		// psi[2]가 1이면, 질의에서의 수선의 발이 최단거리 이므로, 질의의 좌표를 살려야 함
					temp_coord2 = SM_p1(psi[0], node[i].LL[j]);		// psi[0]이 0이면, LL bound의 좌표를 살림
					temp_coord3 = SM_p1(SBN(psi[0]), node[i].RR[j]);		// psi[0]이 1이면, RR bound의 좌표를 살림
					paillier_mul(pubkey, temp_coord3, temp_coord2, temp_coord3);
					temp_coord3 = SM_p1(SBN(psi[2]), temp_coord3);
					paillier_mul(pubkey, shortestPoint[j], temp_coord1, temp_coord3);
				}
				temp_nodedist = DP_SSED(q, shortestPoint, dim);
				paillier_subtract(pubkey, temp_alpha, ciper_one, alpha[i]);
				paillier_mul(pubkey, ciper_nodedist[i], SM_p1(alpha[i], cipher_Max), SM_p1(temp_alpha, temp_nodedist));
				alpha[i] = GSCMP(ciper_nodedist[i], kth_dist);	 // 노드까지의 최단 거리를 계산해서 k번째 결과의 거리와 비교. 노드까지의 거리가 더 작으면 1 셋팅
			}
		}
		
		if(!verify_flag)	{
			endTime = clock();
			node_SRO_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("SRO time : %f\n", node_SRO_time);
		}
		else {
			endTime = clock();
			node_expansion_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node expansion time : %f\n", node_expansion_time);
		}

		cnt = 0;

		if(!verify_flag)	{	// 검증 단계가 아닐 시에는, cand에 저장
			startTime = clock();
			
			cand = GSRO_sNodeRetrievalforkNN(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
		
			endTime = clock();
			data_extract_first_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_first_time);

		}
		else {
			// 검증 단계일 시에는, 이전 결과와 cand를 합침
			startTime = clock();
			cout << "second GSRO_sNodeRetrievalforkNN start" <<endl;
			temp_cand = GSRO_sNodeRetrievalforkNN(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = clock();
			data_extract_second_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_second_time);
			cout << "0" <<endl;
			
			if(cnt == 0) {	// 검증을 위해 추가 탐색이 필요한 노드가 없는 경우를 처리함
				cout<<"break"<<endl;
				break;
			}
			
			
			cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ ){
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = paillier_create_enc_zero();
				}
			}
			// 이전 k개의 결과를 저장
			for(i=0; i<k; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = ciper_result[i][j];
				}
			}
			// 새로 찾은 후보 결과를 저장
			for(i=0; i<cnt; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i+k][j] = temp_cand[i][j];
				}
			}
			cnt = cnt + k;
		}
		paillier_ciphertext_t** ciper_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		for( i = 0 ; i < cnt ; i++ ){
			ciper_distance[i] = DP_SSED(q, cand[i], dim);
			gmp_printf("%d  distance : %Zd\n", i , paillier_dec(0, pub, prv, ciper_distance[i]));
		}
		idx = SkNNb_C2(ciper_distance, k, cnt);
		
		for( i = 0 ; i < k ; i++ ){
			printf("idx :  %d\n", idx[i]);
			for( j = 0 ; j < dim ; j++ ){
				ciper_result[i][j] = cand[idx[i]][j];
			}
		}
		
		if(!verify_flag)	{	
			kth_dist = DP_SSED(q, ciper_result[k-1], dim); // k번째 결과까지의 거리
		}
		
		if(!verify_flag)	{
			verify_flag = true;		// 검증을 수행하기 위해 flag를 true로 변경
		}
		else{
			break;	// 검증을 끝냄
		}
	}
	cout<<"break1"<<endl;
	for( i = 0 ; i < k ; i++ ){
		for( j = 0 ; j < dim ; j++ ){
			paillier_mul(pub, ciper_result[i][j], ciper_result[i][j], ciper_rand);
			gmp_printf("%Zd\t", paillier_dec(0, pub, prv, ciper_result[i][j]));
		}
		cout <<endl;
	}

	free(cipher_Max);
	free(idx);
	
	return SkNNb_Bob(ciper_result, rand, k);
}
// added by KHI. 150808
paillier_ciphertext_t* protocol::GSCMP(paillier_ciphertext_t* u, paillier_ciphertext_t* v){
	//cout << "GSCMP Start" << endl;
	//random value setting = 5
	paillier_ciphertext_t* ciper_Rand_value = paillier_create_enc(5);
	//같은 경우를 판별할 수 필요.
	paillier_plaintext_t* shift = paillier_plaintext_from_ui(2);
	
	//random value 와 arg 인 u , v 를 더해서 저장할 ru , rv에 동적할당
	paillier_ciphertext_t* ru = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	paillier_ciphertext_t* rv = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(ru->c);
	mpz_init(rv->c);

	//ru, rv, rand_value 의 가블드 결과값을 참조할 alpha를 생?? 메모리를 할당할 필요는 없다 왜냐하면 GSCMP_sub에서 alpha라는 변수에 메모리를 할당해 주고 그것의 주소값을 전달하기 때문이다.
	paillier_ciphertext_t* alpha;
	
	//u * 2, v * 2 + 1
	
	paillier_exp(pubkey, ru, u, shift);
	
	paillier_exp(pubkey, rv, v, shift);
	paillier_mul(pubkey, rv, rv, ciper_one); 
	
	paillier_mul(pubkey, ru, ru, ciper_Rand_value); 
	paillier_mul(pubkey, rv, rv, ciper_Rand_value); 
	
	alpha = GSCMP_sub( ru, rv, ciper_Rand_value);
	//메모리 해제
	free(ciper_Rand_value);
	free(ru);
	free(rv);
	//gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, alpha));
	
	return alpha;
}
paillier_ciphertext_t* protocol::GSCMP_sub(paillier_ciphertext_t* ru, paillier_ciphertext_t* rv, paillier_ciphertext_t* ciper_Rand){
	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

	paillier_plaintext_t* plain_U = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	paillier_plaintext_t* plain_V = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	paillier_plaintext_t* plain_R = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	mpz_init(plain_U->m);
	mpz_init(plain_V->m);
	mpz_init(plain_R->m);

	plain_U	= paillier_dec(0, pub, prv, ru);
	plain_V	= paillier_dec(0, pub, prv, rv);
	plain_R	= paillier_dec(0, pub, prv, ciper_Rand);
	
	int a = mpz_get_ui(plain_U->m);
	int b = mpz_get_ui(plain_V->m);
	int R = mpz_get_ui(plain_R->m);
	
	makeGate();
	
	if(G_CMP( b , R , a , R)){
		alpha = paillier_create_enc(1);
	}else{
		alpha = paillier_create_enc(0);
	}

	free(plain_U);
	free(plain_V);
	free(plain_R);
	
	return alpha;
}

paillier_ciphertext_t* protocol::SCMP(paillier_ciphertext_t** u, paillier_ciphertext_t** v)
{
	printf("\n===== Now SCMP starts ===== \n");
	
	int i, j = 0;
	
	bool func = false;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* ciper_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** ciper_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_G		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

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
		
	for(i=0; i<size+1; i++) {
		// compute W
		if(func){	// true :  F : u>v	
			paillier_subtract(pubkey, ciper_W[i], u[i], SM_p1(u[i], v[i]));	
		}else{
			paillier_subtract(pubkey, ciper_W[i], v[i], SM_p1(u[i], v[i]));
		}

		// compute G
		ciper_G[i]=SBXOR(u[i], v[i]);

		// compute H
		if(i==0){
			paillier_exp(pubkey,ciper_H[i], ciper_zero, Rand_value);
			paillier_mul(pubkey,ciper_H[i], ciper_H[i], ciper_G[i]);
		}else{
			paillier_exp(pubkey,temp[i],ciper_H[i-1],Rand_value);
			paillier_mul(pubkey,ciper_H[i],temp[i],ciper_G[i]);
		}
	}

	// debugging
	/*
	printf("\nW : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_W[i]));
	}
	printf("\nG : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_G[i]));
	}

	printf("\nH : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_H[i]));
	}
	*/

	for(i=0; i<size+1; i++){
		// compute PI
		paillier_mul(pubkey, ciper_PI[i], ciper_H[i], ciper_minus);	// PI

		// compute L with enhanced privacy (new IDEA)
		mpz_init(temp[i]->c);
		paillier_exp(pubkey,temp[i],ciper_PI[i],Rand_value);	// randomize PI
		paillier_exp(pubkey,temp2[i],ciper_W[i],Rand_value); // randomize W
		paillier_mul(pubkey,ciper_L[i],temp[i], temp2[i]);
	}
	
	// debugging
	/*
	printf("\nPI : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_PI[i]));
	}

	printf("\nL : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_L[i]));
	}
	printf("\n");
	*/

	alpha = SRO2(ciper_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}
	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(ciper_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size+1; i++){
		paillier_freeciphertext(ciper_W[i]);	paillier_freeciphertext(ciper_H[i]);	paillier_freeciphertext(ciper_L[i]);	paillier_freeciphertext(ciper_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return alpha;
	
}

paillier_ciphertext_t* protocol::SCMP_M(paillier_ciphertext_t** u, paillier_ciphertext_t** v)
{
	int i, j = 0;
	
	bool func = false;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* ciper_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** ciper_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_G		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));

	for(i=0; i<size; i++){
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

	for(i=0; i<size; i++) {
		// compute W
		if(func){	// true :  F : u>v	
			paillier_subtract(pubkey, ciper_W[i], u[i], SM_p1(u[i], v[i]));	
		}else{
			paillier_subtract(pubkey, ciper_W[i], v[i], SM_p1(u[i], v[i]));
		}
		// compute G
		ciper_G[i]=SBXOR(u[i], v[i]);
		// compute H
		if(i==0){
			paillier_exp(pubkey,ciper_H[i], ciper_zero, Rand_value);
			paillier_mul(pubkey,ciper_H[i], ciper_H[i], ciper_G[i]);
		}else{
			paillier_exp(pubkey,temp[i],ciper_H[i-1],Rand_value);
			paillier_mul(pubkey,ciper_H[i],temp[i],ciper_G[i]);
		}
	}

	// debugging
	/*
	printf("\nW : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_W[i]));
	}
	printf("\nG : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_G[i]));
	}

	printf("\nH : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_H[i]));
	}
	*/
	for(i=0; i<size; i++){
		// compute PI
		paillier_mul(pubkey, ciper_PI[i], ciper_H[i], ciper_minus);	// PI

		// compute L with enhanced privacy (new IDEA)
		mpz_init(temp[i]->c);
		paillier_exp(pubkey,temp[i],ciper_PI[i],Rand_value);	// randomize PI
		paillier_exp(pubkey,temp2[i],ciper_W[i],Rand_value); // randomize W
		paillier_mul(pubkey,ciper_L[i],temp[i], temp2[i]);
	}
	
	// debugging
	/*
	printf("\nPI : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_PI[i]));
	}

	printf("\nL : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_L[i]));
	}
	printf("\n");
*/
	alpha = SRO2(ciper_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}
	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(ciper_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size; i++){
		paillier_freeciphertext(ciper_W[i]);	paillier_freeciphertext(ciper_H[i]);	paillier_freeciphertext(ciper_L[i]);	paillier_freeciphertext(ciper_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return alpha;
	
}
paillier_ciphertext_t* protocol::SCMP_Clustering(paillier_ciphertext_t** u, paillier_ciphertext_t** v)
{
	int i, j = 0;
	
	bool func = true;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* ciper_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** ciper_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_G		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** ciper_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));

	for(i=0; i<size; i++){
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

	for(i=0; i<size; i++) {
		// compute W
		if(!func){	// true :  F : u>v	
			paillier_subtract(pubkey, ciper_W[i], u[i], SM_p1(u[i], v[i]));	
		}else{
			paillier_subtract(pubkey, ciper_W[i], v[i], SM_p1(u[i], v[i]));
		}
		// compute G
		ciper_G[i]=SBXOR(u[i], v[i]);
		// compute H
		if(i==0){
			paillier_exp(pubkey,ciper_H[i], ciper_zero, Rand_value);
			paillier_mul(pubkey,ciper_H[i], ciper_H[i], ciper_G[i]);
		}else{
			paillier_exp(pubkey,temp[i],ciper_H[i-1],Rand_value);
			paillier_mul(pubkey,ciper_H[i],temp[i],ciper_G[i]);
		}
	}

	// debugging
	/*
	printf("\nW : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_W[i]));
	}
	printf("\nG : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_G[i]));
	}

	printf("\nH : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_H[i]));
	}
	*/
	for(i=0; i<size; i++){
		// compute PI
		paillier_mul(pubkey, ciper_PI[i], ciper_H[i], ciper_minus);	// PI

		// compute L with enhanced privacy (new IDEA)
		mpz_init(temp[i]->c);
		paillier_exp(pubkey,temp[i],ciper_PI[i],Rand_value);	// randomize PI
		paillier_exp(pubkey,temp2[i],ciper_W[i],Rand_value); // randomize W
		paillier_mul(pubkey,ciper_L[i],temp[i], temp2[i]);
	}
	
	// debugging
	/*
	printf("\nPI : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_PI[i]));
	}

	printf("\nL : ");
	for(i=0; i<size; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_L[i]));
	}
	printf("\n");
*/
	alpha = SRO2(ciper_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}
	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(ciper_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size; i++){
		paillier_freeciphertext(ciper_W[i]);	paillier_freeciphertext(ciper_H[i]);	paillier_freeciphertext(ciper_L[i]);	paillier_freeciphertext(ciper_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return alpha;
	
}
paillier_ciphertext_t* protocol::SCMP_for_SBD(paillier_ciphertext_t** u, paillier_ciphertext_t** v)
{
	//printf("\n===== Now SRO starts ===== \n");
	
	int i, j = 0;
	
	bool func = false;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* ciper_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** ciper_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_G		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** ciper_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

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
		
	for(i=0; i<size; i++) {
		// compute W
		if(func){	// true :  F : u>v	
			paillier_subtract(pubkey, ciper_W[i], u[i], SM_p1(u[i], v[i]));	
		}else{
			paillier_subtract(pubkey, ciper_W[i], v[i], SM_p1(u[i], v[i]));
		}

		// compute G
		ciper_G[i]=SBXOR(u[i], v[i]);

		// compute H
		if(i==0){
			paillier_exp(pubkey,ciper_H[i], ciper_zero, Rand_value);
			paillier_mul(pubkey,ciper_H[i], ciper_H[i], ciper_G[i]);
		}else{
			paillier_exp(pubkey,temp[i],ciper_H[i-1],Rand_value);
			paillier_mul(pubkey,ciper_H[i],temp[i],ciper_G[i]);
		}
	}

	// debugging
	/*
	printf("\nW : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_W[i]));
	}
	printf("\nG : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_G[i]));
	}

	printf("\nH : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_H[i]));
	}
	*/

	for(i=0; i<size; i++){
		// compute PI
		paillier_mul(pubkey, ciper_PI[i], ciper_H[i], ciper_minus);	// PI

		// compute L with enhanced privacy (new IDEA)
		mpz_init(temp[i]->c);
		paillier_exp(pubkey,temp[i],ciper_PI[i],Rand_value);	// randomize PI
		paillier_exp(pubkey,temp2[i],ciper_W[i],Rand_value); // randomize W
		paillier_mul(pubkey,ciper_L[i],temp[i], temp2[i]);
	}
	
	// debugging
	/*
	printf("\nPI : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_PI[i]));
	}

	printf("\nL : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_L[i]));
	}
	printf("\n");
	*/

	alpha = SRO2(ciper_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}
	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(ciper_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size+1; i++){
		paillier_freeciphertext(ciper_W[i]);	paillier_freeciphertext(ciper_H[i]);	paillier_freeciphertext(ciper_L[i]);	paillier_freeciphertext(ciper_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return alpha;
	
}


paillier_ciphertext_t* protocol::DP_SSED(paillier_ciphertext_t** ciper1, paillier_ciphertext_t** ciper2, int col_num)
{
	paillier_plaintext_t* shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* value = paillier_plaintext_from_ui(0);
	paillier_plaintext_t** rand = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*col_num);
	paillier_plaintext_t* tmp= paillier_plaintext_from_ui(0);
	paillier_ciphertext_t* result = paillier_create_enc_zero();

	paillier_ciphertext_t* ciper_value;
	//paillier_ciphertext_t* ciper_tmp= paillier_create_enc_zero();	
	paillier_ciphertext_t* ciper_tmp = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));	
	mpz_init(ciper_tmp->c);
	
	paillier_ciphertext_t* sub =  paillier_create_enc_zero();
/*
	paillier_ciphertext_t** result_array=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*col_num);
	for(int i=0;i<col_num;i++){
		result_array[i]	= paillier_create_enc_zero();
	}
*/

	for(int i=0;i<col_num;i++){
		mpz_ui_pow_ui(shift->m,2,16*i);	
		rand[i] = paillier_plaintext_from_ui(i+100);
		//gmp_printf("rand : %Zd\n",rand[i]);
		
		mpz_mul(tmp->m,rand[i]->m,shift->m);

		mpz_add(value->m,value->m, tmp->m);

		//gmp_printf("value : %Zd\n",value);
	}

	
	ciper_value =  paillier_enc(0, pubkey, value, paillier_get_rand_devurandom);	

	for(int i=0;i<col_num;i++){
		mpz_ui_pow_ui(shift->m,2,16*i);			//shift
		//gmp_printf("shift : %Zd\n",shift);

		paillier_subtract(pubkey,sub,ciper1[i],ciper2[i]); //sub
		//gmp_printf("sub : %Zd\n",paillier_dec(0, pubkey, prvkey, sub));

		//gmp_printf("SM_p1 : %Zd\n",paillier_dec(0, pubkey, prvkey,SM_p1(sub,sub)));
		
		paillier_exp(pubkey, ciper_tmp, sub, shift); //multi		
		//gmp_printf("ciper_tmp : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_tmp));

		paillier_mul(pubkey,ciper_value,ciper_value,ciper_tmp); //E_add
		//gmp_printf("value : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_value));
	}
	
		
	//result_array=unDP_SSED(ciper_value,col_num);
	//paillier_ciphertext_t* test = unDP_SSED(ciper_value,col_num);		// 확인 요망!!

	result = unDP_SSED(ciper_value,col_num);

	//gmp_printf("unDP_SSED result : %Zd\n",paillier_dec(0, pubkey, prvkey, test));

	for(int i=0; i<col_num; i++){
		mpz_mul(tmp->m,rand[i]->m,rand[i]->m);
		ciper_tmp=paillier_enc(0,pub,tmp,paillier_get_rand_devurandom);
		paillier_subtract(pub,result,result,ciper_tmp );

		mpz_mul_ui (rand[i]->m,rand[i]->m,2);
		paillier_subtract(pub,ciper_tmp,ciper1[i],ciper2[i]);
		paillier_exp(pubkey, sub, ciper_tmp, rand[i]);
		paillier_subtract(pub,result,result,sub);
	}
	//free(shift);free(value);free(tmp);
	//gmp_printf("result : %Zd\n",paillier_dec(0, pubkey, prvkey, result));

	return result;		
}

paillier_ciphertext_t* protocol::unDP_SSED(paillier_ciphertext_t* ciper1, int col_num){
	paillier_plaintext_t* shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* tmp = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* result = paillier_dec(0, pub, prv, ciper1);

	paillier_plaintext_t* test = paillier_plaintext_from_ui(0);
	paillier_ciphertext_t* enc_test = paillier_create_enc_zero();
	
	
	//gmp_printf("unDP_SSED val : %Zd\n",result);

	for(int i=0; i<col_num; i++){
		mpz_ui_pow_ui(shift->m,2,(col_num-i-1)*16);
		//gmp_printf("unDP_SSED shift : %Zd\n",shift);
		
		//gmp_printf("unDP_SSED   result : %Zd , shift: %Zd\n",result,shift);
		mpz_div(tmp->m,result->m,shift->m);
		
		//gmp_printf("div_tmp : %Zd\n",tmp);


		mpz_pow_ui(tmp->m,tmp->m,2);
		
		//result_array[col_num-i-1]=paillier_enc(0,pub,tmp,paillier_get_rand_devurandom);
		mpz_add(test->m,test->m,tmp->m);
		
		//gmp_printf("unDP_SSED re : %Zd\n",tmp);

		mpz_mod(result->m,result->m,shift->m);
	}
	enc_test=paillier_enc(0,pub,test,paillier_get_rand_devurandom);
	/*
	for(int i=0; i<col_num; i++){
		gmp_printf("unDP_SSED re : %Zd\n",paillier_dec(0, pub, prv, result_array[i]));
	}
	*/
	//free(shift);free(tmp);
	return enc_test;
	//return result_array;
}

int protocol::SMSn(paillier_ciphertext_t** ciper, int cnt){
	int i = 0, s = 0;
	int buff = 0;
	int R1 = 5;
	int quot = 0;
	int remain = 0;
	int iter = 0;
	int bundle = 0;
	int ciper_idx = 0;
	Print = false;
	bool re = false;
	int * R_array_result = (int *)malloc(sizeof(int)*cnt);
	int * V_array_result = (int *)malloc(sizeof(int)*cnt);
	for( i = 0 ; i < cnt ; i++ )
	{
		V_array_result[i] = 0;
		R_array_result[i] = 0;
	}
	if(modul == 1024 && DP != 32)
	{
		buff = 64;
	}else if(modul == 512 && DP != 32)
	{
		buff = 32;
	}else if( DP == 32 )
	{
		buff = 16;
	}
	quot	= (int)cnt/buff;
	remain	= (int)cnt%buff;

	if(Print)
	{
		cout << " quot : "<< quot << " remain : " << remain << "  " << buff <<endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!SMSn Start!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		for( i = 0 ; i < 200 ; i++ ){
			cout << i+1 << " ciper : ";
			gmp_printf(" %Zd\n",paillier_dec(0, pub, prv, ciper[i]));
		}
	}

	paillier_plaintext_t* shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* tmp= paillier_plaintext_from_ui(0);
	paillier_plaintext_t* p1_array = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* R1_array = paillier_plaintext_from_ui(R1);
	paillier_ciphertext_t* c1_array ;
	paillier_ciphertext_t* cR1_array;
	paillier_ciphertext_t* ciper_tmp = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));	
	mpz_init(ciper_tmp->c);

	paillier_ciphertext_t* SMSn_sub_result;
	//cout<< "setting iter, bundle" <<endl;
	if ( cnt < buff )
	{
		iter = 1;
		bundle = remain;
		re = true;
	}else if ( (cnt >= buff) && (remain == 0) )
	{
		iter = quot;
		bundle = buff;
		re = false;
	}else if ( (cnt >= buff) && (remain > 0) )
	{
		iter = quot + 1;
		bundle = buff;
		re = true;
	}
	//cout << "start iter"<<endl;
	for( s = 0 ; s < iter ; s++ )
	{
		//cout << s+1 << " th iter start" <<endl;
		if( s == quot && (re == true) )
		{
			bundle = remain;
		}
		//cout << bundle << " " <<iter<< " " << DP << " " << buff << " "<< remain << endl;
		for(i=0;i<bundle;i++){
			mpz_ui_pow_ui(shift->m,2,DP*i);
			mpz_mul(tmp->m,R1_array->m,shift->m);
			mpz_add(p1_array->m,p1_array->m, tmp->m);
		}
		c1_array = paillier_enc(0, pubkey, p1_array, paillier_get_rand_devurandom);
		cR1_array = paillier_enc(0, pubkey, p1_array, paillier_get_rand_devurandom);
		for( i = 0 ; i < bundle ; i++){
			mpz_ui_pow_ui(shift->m,2,DP*i);
			paillier_exp(pub, ciper_tmp, ciper[ciper_idx++], shift);
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
		}
		SMSn_sub_result = SMSn_sub(c1_array, cR1_array, bundle);
		R1_array = paillier_dec(0, pub, prv, SMSn_sub_result);
		// UNPacking
		for( i = 0 ; i < bundle; i++){
			//cout << "idx "<<(s*buff+bundle)-1-i<<endl;
			mpz_ui_pow_ui(shift->m,2,(bundle-i-1)*DP);
			mpz_div(tmp->m, R1_array->m, shift->m);
			//gmp_printf("tmpR : %Zd\n", tmp);
			R_array_result[(s*buff+bundle)-1-i] = mpz_get_ui(tmp->m);
			mpz_mod(R1_array->m,R1_array->m,shift->m);
			mpz_div(tmp->m, SMSn_value->m, shift->m);
			//gmp_printf("tmpV : %Zd\n", tmp);
			V_array_result[(s*buff+bundle)-1-i] = mpz_get_ui(tmp->m);
			mpz_mod(SMSn_value->m,SMSn_value->m,shift->m);
		}
		mpz_init(ciper_tmp->c);
		shift = paillier_plaintext_from_ui(0);
		tmp= paillier_plaintext_from_ui(0);
		p1_array = paillier_plaintext_from_ui(0);
		R1_array = paillier_plaintext_from_ui(R1);
	}
	//cout << " quot : "<< quot << " remain : " << remain << endl;


	for(i=0; i<cnt; i++){
		R_array_result[i] = -R_array_result[i];
	}

	if(Print)
	{
		cout<< "~~~~~~~~~~~~~~~~~~~~~~~Data Packing After~~~~~~~~~~~~~~~~~~~~~~~~\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		cout << " quot : "<< quot << " remain : " << remain << endl;
		for( i = 0 ; i < 200 ; i++ )
		{
			cout << i+1 << " ciper : "<< R_array_result[i] + V_array_result[i] <<endl;
		}
	}
	
	if(!smsn_flag){
		smsn_flag = true;
		makeGate();  // 매번해야하는지 ? 체크
	}
	smsn_result_idx = 0;
	smsn_tmp_val = V_array_result[0];
	smsn_tmp_rand = R_array_result[0];
	
	for( i = 0 ; i < cnt ; i++ )
	{
		//cout << smsn_tmp_val+smsn_tmp_rand << " ";
		if(G_CMP(smsn_tmp_val, smsn_tmp_rand, V_array_result[i], R_array_result[i])){
			//cout << ">";
			smsn_tmp_val    = V_array_result[i];	
			smsn_tmp_rand   = R_array_result[i];
			smsn_result_idx = i;
		}else{
			//cout << "<";
		}
		//cout << V_array_result[i]+R_array_result[i] <<endl;
	}
	
	free(shift);
	free(tmp);
	free(p1_array);
	free(R1_array);
	free(c1_array);
	free(cR1_array);
	free(ciper_tmp);
	
	if(Print)
	{
	
		cout << "Min Value : " << smsn_tmp_val+smsn_tmp_rand <<"   idx : "<<smsn_result_idx;
		gmp_printf("  %Zd\n",paillier_dec(0, pub, prv, ciper[smsn_result_idx]));	
	}
	
	return smsn_result_idx;
}

paillier_ciphertext_t* protocol::SMSn_sub(paillier_ciphertext_t* ciper, paillier_ciphertext_t* rand ,int cnt){
	int i = 0;
	int * R1 = (int *)malloc(sizeof(int)*cnt);
	for( i = 0 ; i < cnt ; i ++){
		R1[i] = 4;
	}

	paillier_plaintext_t* shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* tmp = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* rand_tmp = paillier_plaintext_from_ui(0);

	paillier_plaintext_t** R1_array = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*cnt);

	paillier_plaintext_t* value_array = paillier_dec(0, pub, prv, ciper);  //A에서 전달되는 packing 값
	paillier_plaintext_t* rand_array = paillier_dec(0, pub, prv, rand); 

	for(i=0; i<cnt; i++){
		mpz_ui_pow_ui(shift->m,2,DP*i);
		R1_array[i] = paillier_plaintext_from_ui(R1[i]);
		mpz_mul(tmp->m,R1_array[i]->m,shift->m);
		mpz_add(value_array->m,value_array->m, tmp->m);
		mpz_add(rand_array->m,rand_array->m, tmp->m);
	}

	SMSn_value = value_array;
	rand = paillier_enc(0, pub, rand_array, paillier_get_rand_devurandom);

	return rand;
}

int ** protocol::SkNN_b(paillier_ciphertext_t*** ciper, paillier_ciphertext_t** query, int k, int row_number){
	int i = 0, j = 0;
	int rand = 10;
	int * idx = (int *)malloc(sizeof(int)*k);
	paillier_ciphertext_t* ciper_rand = paillier_create_enc(rand);
	paillier_ciphertext_t** ciper_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*row_number);
	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t)*k);

	paillier_ciphertext_t** ciper_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	
	for( i = 0 ; i < dim ; i++ ){
		ciper_query[i] = query[i];
	}
	for( i = 0 ; i < k ; i++){
		ciper_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	}
	
	
	for( i = 0 ; i < row_number ; i++ ){
		ciper_distance[i] = SSEDm(ciper[i] ,ciper_query, dim);
		//paillier_print("ciper_dist : ", ciper_distance[i]);
	}
	idx = SkNNb_C2(ciper_distance, k, row_number);
	
	for( i = 0 ; i < k ; i++ ){
		//printf("idx :  %d\n", idx[i]);
		for( j = 0 ; j < dim ; j++ ){
			ciper_result[i][j] = ciper[idx[i]][j];
			paillier_mul(pub, ciper_result[i][j], ciper_result[i][j], ciper_rand);
			//paillier_print(" ",ciper_result[i][j]);
		}
		//printf("\n");
	}
	node_SBD_time=0;
	node_SRO_time=0;
	node_expansion_time=0;
	data_extract_first_time=0;
	data_extract_second_time=0;
	data_SSED_SBD_time=0;
	sMINn_first_time=0;
	sMINn_second_time=0;
	data_SBOR_time=0;
	data_SPE_time=0;
	totalNumOfRetrievedNodes=1;
	free(idx);
	free(ciper_distance);
	free(ciper_query);
	
	return SkNNb_Bob(ciper_result, rand, k);
}
int ** protocol::SkNNb_Bob(paillier_ciphertext_t*** ciper_result, int rand, int k){
	int i = 0, j = 0;
	int ** result = (int **)malloc(sizeof(int*)*k);
	paillier_plaintext_t* plain;
	for( i = 0 ; i < k ; i ++){
		result[i] = (int *)malloc(sizeof(int)*dim);
		for( j = 0 ; j < dim ; j ++){
			plain = paillier_dec(0, pub, prv, ciper_result[i][j]);
			result[i][j] = mpz_get_si(plain->m);
			result[i][j] = result[i][j] - rand;
			//printf("%d ", result[i][j]);
			free(ciper_result[i][j]);
		}
		//printf("\n");
	}

	return result;
}
int * protocol::SkNNb_C2(paillier_ciphertext_t** ciper_dist, int k, int row_number){
	int i = 0 , j = 0;
	mpz_t temp;
	paillier_plaintext_t** plain = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*row_number);
	int * result = (int *)malloc(sizeof(int)*k);
	int * idx = (int *)malloc(sizeof(int)*row_number);
	int * toint = (int *)malloc(sizeof(int)*row_number);

	for( i = 0 ; i < k; i++){
		result[i] = 0;
	//	printf("%d\n", result[i]);
	}
	for( i = 0 ; i < row_number; i ++){
		toint[i] = 0;
		idx[i] = i;
		plain[i] = paillier_dec(0, pub, prv, ciper_dist[i]);
		toint[i] = mpz_get_si(plain[i]->m);
		//gmp_printf("dec :   %d : %Zd\n",  toint[i] , plain[i]);
	}

	for( i = 0 ; i < row_number-1 ; i++){
		for( j = i+1 ; j < row_number ; j ++ ){
			if( toint[i] > toint[j] ){
				int tmp = toint[i];
				toint[i] = toint[j];
				toint[j] = tmp;
				tmp = idx[i];
				idx[i] = idx[j];
				idx[j] = tmp;
			}
		}
	}
/*
	for( i = 0 ; i < row_number; i ++){
		gmp_printf("dec : %d : %d\n", idx[i] , toint[i]);
	}
*/
	for( i = 0 ; i < k ; i++ ){
		result[i] = idx[i];
		//printf("result :  %d\n", result[i]);
	}
	free(plain);
	free(idx);
	free(toint);
	return result;
}



int * protocol::SkNN_B(paillier_ciphertext_t** ciper, int Q, int k, int row_number){
	int i=0, s=0, n=0, t=0, j=0;
	int rand = 5;
	n = row_number;
	
	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* ciper_binary;
	paillier_ciphertext_t* ciper_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_dist	= paillier_create_enc_zero();
	
	printf("%d\n", Q);
	printf("%d\n", k);
	printf("%d\n", row_number);
	
	paillier_ciphertext_t* ciper_query = paillier_create_enc(Q);

	paillier_print("query : ",ciper_query);
	for( i = 0 ; i < n ; i++){
		paillier_print("i : ",ciper[i]);
	}

	paillier_ciphertext_t** ciper_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** ciper_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** ciper_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** ciper_V;
	paillier_ciphertext_t** ciper_V2 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	
	paillier_ciphertext_t*  ciper_result = paillier_create_enc_zero();
	paillier_ciphertext_t** ciper_result_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	
	paillier_ciphertext_t*** ciper_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*n);
	for( i = 0 ; i < n ; i++ ){
		ciper_mid[i] = paillier_create_enc_zero();	// 초기화만 해줘도 될듯!!??
		ciper_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		for( j = 0 ; j < size ; j ++ ){
			ciper_SBD_distance[i][j] = paillier_create_enc_zero();	// 초기화만 해줘도 될듯!!??
		}
	}

	for( i = 0 ; i < size; i++ ){
		ciper_Smin[i] = ciper_zero;
	}
	
	paillier_ciphertext_t* ciper_rand = paillier_create_enc(rand);
	
	paillier_print("rand : ",ciper_rand);
	
	for( i = 0 ; i < n ; i++ ){
		ciper_distance[i] = SSED(ciper_query, ciper[i]);
		paillier_print("dist : ", ciper_distance[i]);	
		ciper_SBD_distance[i] = SBD(ciper_distance[i]);
		for( j = 0 ; j < size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_SBD_distance[i][j]));
		}
		printf("\n");
	}
	for( s = 0 ; s < k ; s++ ){
		
		ciper_Smin = Smin_n(ciper_SBD_distance, n);

		for( j = 0 ; j < size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smin[j]));
		}
		printf("\n");	
		
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			ciper_binary = paillier_create_enc(t);
			ciper_binary = SM_p1(ciper_binary, ciper_Smin[size-j]);
			paillier_mul(pubkey, ciper_min, ciper_binary, ciper_min);

			
		}
		paillier_print("min : ",ciper_min);
		
		if( s != 0 ){
			for( i = 0 ; i < n ; i++ ){
				for( j = size ; j > 0 ; j--){
					t = (int)pow(2, j-1);
					ciper_binary = paillier_create_enc(t);
					ciper_binary = SM_p1(ciper_binary, ciper_SBD_distance[i][size-j]);
					paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
				}
				ciper_distance[i] = ciper_dist;
				ciper_dist = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	 // 전역 ciper_zero 사용해도 되지 않을까!!??
			}
		}

	
		for( i = 0 ; i < n ; i++ ){
			paillier_print("distance : ",ciper_distance[i]);
		}
		
		for( i = 0 ; i < n ; i++ ){
			paillier_ciphertext_t* ciper_temp;
			paillier_subtract(pubkey, ciper_temp, ciper_distance[i], ciper_min);
			ciper_mid[i] = SM_p1(ciper_temp, ciper_rand);
			paillier_print("dist - dist : ",ciper_temp);
		}
		ciper_V = SkNNm_sub(ciper_mid, n);
		for( i = 0 ; i < n ; i++ ){
			paillier_print("ciper_V : ",ciper_V[i]);
		}
		

		for( i = 0 ; i < n ; i++ ){
			ciper_V2[i] = SM_p1(ciper_V[i], ciper[i]);
			paillier_mul(pubkey, ciper_result, ciper_V2[i], ciper_result);
		}
		paillier_print("ciper_result : ",ciper_result);
		
		ciper_result_array[s] = ciper_result;
		for( i = 0 ; i < n ; i++ ){
			for( j = 0; j < size ; j++ ){
				ciper_SBD_distance[i][j] = SBOR(ciper_V[i], ciper_SBD_distance[i][j]);
			}
		}
		ciper_result = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);
		ciper_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);
	}

	for( i = 0 ; i < n ; i++){
		for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_SBD_distance[i][j]));
			}
		printf("\n");
	}

	for( i = 0 ; i < k ; i++ ){
		paillier_mul(pubkey, ciper_result_array[i],ciper_result_array[i], ciper_rand);
		paillier_print("ciper_result_array : ",ciper_result_array[i]);
	}
	free(ciper_binary);
	free(ciper_min);
	free(pt);
	free(ciper_dist);
	free(ciper_mid);
	free(ciper_Smin);
	free(ciper_V);
	free(ciper_V2);
	free(ciper_result);
	free(ciper_SBD_distance);

	return SkNNm_Bob(ciper_result_array, rand, k);
}

// SkNN subroutine
paillier_ciphertext_t** protocol::SkNNm_sub(paillier_ciphertext_t** ciper_n, int n){
	int i = 0;
	int cnt = 0;
	paillier_plaintext_t* plain = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	paillier_ciphertext_t** ciper_U = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	
	for( i = 0 ; i < n ; i++ ){
		plain = paillier_dec(0, pubkey, prvkey, ciper_n[i]);
		if(  (mpz_cmp(plain_zero->m, plain->m) == 0) && cnt == 0 ){
			ciper_U[i] = paillier_create_enc(1);
			cnt++;
		}else if ((mpz_cmp(plain_zero->m, plain->m) == 0))
		{
			ciper_U[i] = paillier_create_enc(0);
			cnt++;
		}else if((mpz_cmp(plain_zero->m, plain->m) != 0)){
			ciper_U[i] = paillier_create_enc(0);
		}
	}
	if(cnt == n){
		for( i = 0 ; i < n ; i++ ){
			ciper_U[i] = paillier_create_enc(0);
		}
	}
	
	paillier_freeplaintext(plain);

	return ciper_U;
}

// 일차원 SkNN 결과 return 함수
int * protocol::SkNNm_Bob(paillier_ciphertext_t** ciper_result_array, int rand, int k){
	int * kNN = (int *)malloc(sizeof(int)* k);
	int i = 0;
	paillier_plaintext_t* plain;
	for( i = 0 ; i < k ; i++ ){
		plain = paillier_dec(0, pubkey, prvkey, ciper_result_array[i]);
		kNN[i] = mpz_get_ui(plain->m);
		kNN[i] = kNN[i] - rand;
	}
	free(plain);
	free(ciper_result_array);
	for( i = 0 ; i< k ; i++ ){
		printf("result : %d\n", kNN[i]);
	}
	return kNN;
}

// 다차원 SkNN 결과 return 함수
int ** protocol::SkNNm_Bob2(paillier_ciphertext_t*** ciper_result, int rand, int k, int col_num){
	int i = 0 , j = 0;
	paillier_plaintext_t* plain;

	int ** kNN = (int **)malloc(sizeof(int*)*k);
	for( i = 0 ; i < k ; i++ )
		kNN[i] = (int*)malloc(sizeof(int)*col_num);
	
	for( i = 0 ; i < k ; i++ ){
		for( j = 0 ; j < col_num ; j++ ){
			//paillier_subtract(pubkey, ciper_result[i][j], ciper_result[i][j], ciper_rand);
			plain = paillier_dec(0, pubkey, prvkey, ciper_result[i][j]);
			kNN[i][j] = mpz_get_ui(plain->m);
			kNN[i][j] = kNN[i][j] - rand ; 
		}
	}

	/*
	for( i = 0 ; i < k ; i++ ){
		printf("%d result : ", (i+1) );
		for ( j = 0 ; j < col_num ; j++ ){
			printf("%d \t", kNN[i][j]);
		}
		printf("\n");
	}
	*/
	paillier_freeplaintext(plain);
	for(i=0; i<k; i++) {
		for(j=0; j<dim; j++) {
			paillier_freeciphertext(ciper_result[i][j]);
		}
		free(ciper_result[i]);
	}
	free(ciper_result);

	return kNN;
}

// 다차원 배열 SkNN
int** protocol::SkNN_B(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, int k, int row_number){
	int i=0, s=0, n=0, t=0, j=0;
	int rand = 5;
	n = row_number;
	
	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* ciper_binary;
	paillier_ciphertext_t* ciper_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_dist	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_dist = paillier_create_enc_zero();

	paillier_ciphertext_t** ciper_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t*** ciper_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*n);
	paillier_ciphertext_t** ciper_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** ciper_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** ciper_V=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t*** ciper_V2 = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*n);

	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	for( i = 0 ; i < size; i++ ){
		ciper_Smin[i] = ciper_zero;
	}

	for( i = 0 ; i < n ; i++ ){
		ciper_distance[i] 	= ciper_zero;
		ciper_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		ciper_V[i]	= ciper_zero;
		ciper_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	}

	for( i = 0 ; i < k ; i++ ){
		ciper_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim; j ++ ){
			ciper_result[i][j] = paillier_create_enc_zero();
		}
	}


	printf("\n=== Data SSED & SBD start ===\n");
	startTime = clock();
	for( i = 0 ; i < n ; i++ ){
		ciper_distance[i] = SSEDm(query, data[i], dim);
		ciper_SBD_distance[i] = SBD(ciper_distance[i]);
		
		/*
		paillier_print("dist : ", ciper_distance[i]);	
		for( j = 0 ; j < size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_SBD_distance[i][j]));
		}
		printf("\n");
		*/		
	}
	endTime = clock();
	data_SSED_SBD_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("data SSED & SBD time : %f\n", (float)(endTime-startTime)/(CLOCKS_PER_SEC));

	for( s = 0 ; s < k ; s++ ){
		printf("\n%dth sMINn start \n", s+1);
		startTime = clock();
		ciper_Smin = Smin_n(ciper_SBD_distance, n);	// bit로 표현된 암호화 min 거리 추출
		endTime = clock();
		gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		sMINn_first_time += gap;
		printf("%dth sMINn time : %f\n", s+1, gap);
	
		/*
		for( j = 0 ; j < size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smin[j]));
		}
		printf("\n");	
		*/
		
		// bit로 표현된 암호화 min 거리를 암호화 정수로 변환
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			ciper_binary = paillier_create_enc(t);
			ciper_binary = SM_p1(ciper_binary, ciper_Smin[size-j]);
			paillier_mul(pubkey, ciper_min, ciper_binary, ciper_min);
			paillier_freeciphertext(ciper_binary);
			//paillier_print("min : ", ciper_min);
		}
		paillier_print("min dist : ", ciper_min);
		
		printf("\n== recalculate query<->data distances ===\n");
		// 이전 iteration에서 min값으로 선택된 데이터의 거리가 secure하게 MAX로 변환되었기 때문에, 
		// 질의-데이터 간 거리 계산을 모두 다시 수행
		if( s != 0 ){
			for( i = 0 ; i < n ; i++ ){
				for( j = size ; j > 0 ; j--){
					t = (int)pow(2, j-1);
					ciper_binary = paillier_create_enc(t);
					ciper_binary = SM_p1(ciper_binary, ciper_SBD_distance[i][size-j]);
					paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
				}
				ciper_distance[i] = ciper_dist;
				ciper_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
			}
		}
		
		/*
		for( i = 0 ; i < n ; i++ ){
			paillier_print("distance : ", ciper_distance[i]);
		}
		*/

		// 질의-데이터 거리와 min 거리와의 차를 구함 (min 데이터의 경우에만 0으로 만들기 위함)
		for( i = 0 ; i < n ; i++ ){
			paillier_subtract(pubkey, temp_dist, ciper_distance[i], ciper_min);
			ciper_mid[i] = SM_p1(temp_dist, ciper_rand);
			//paillier_print("dist - dist : ",ciper_mid[i]);
		}

		ciper_V = SkNNm_sub(ciper_mid, n);

		/*
		for( i = 0 ; i < n ; i++ ){
			printf("%d : ", i);
			paillier_print("ciper_V : ", ciper_V[i]);
		}
		*/

		// min 데이터 추출
		for( i = 0 ; i < n ; i++ ){
			for( j = 0 ; j < dim; j++ ){
				ciper_V2[i][j] = SM_p1(ciper_V[i], data[i][j]);
				paillier_mul(pubkey, ciper_result[s][j], ciper_V2[i][j], ciper_result[s][j]);
			}
		}

		printf("ciper_result : ");	
		for( j = 0 ; j < dim; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey,ciper_result[s][j]));
		}
		printf("\n");		
		
		// Data SBOR 수행
		startTime = clock();
		for( i = 0 ; i < n ; i++ ){
			for( j = 0; j < size ; j++ ){
				ciper_SBD_distance[i][j] = SBOR(ciper_V[i], ciper_SBD_distance[i][j]);
			}
		}
		endTime = clock();
		data_SBOR_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);

		ciper_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
	}
	
	// user(Bob)에게 결과 전송을 위해 random 값 삽입
	for(i=0;i<k;i++){
		//printf("%d final result : ", i);
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,ciper_result[i][j],ciper_result[i][j],ciper_rand);
			//gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, ciper_result[i][j]));
		}
		//printf("\n");
	}

	paillier_freeciphertext(temp_dist);	
	
	return SkNNm_Bob2(ciper_result, rand, k, dim);
}



paillier_ciphertext_t*** protocol::sNodeRetrievalforkNN(paillier_ciphertext_t*** data, paillier_ciphertext_t*** ciper_qLL_bit, paillier_ciphertext_t*** ciper_qRR_bit, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now sNodeRetrievalforkNN starts =====\n");
	//printf("NumNode : %d\n", NumNode);
	int i=0, j=0, m=0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	int** node_group;

	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);
	if(*NumNodeGroup == 0)
		return 0;

	printf("set_num : %d\n", *NumNodeGroup);

	for(i=0; i<*NumNodeGroup; i++) {
		printf("%dth Node Group : ", i+1);
		for(j=1; j<=node_group[i][0]; j++) {	// 0번지에 해당 노드 그룹에 몇개의 노드가 있는지가 저장되어 있음
			printf("%d ", node_group[i][j]);
		}
		printf("\n");
	}

	int nodeId = 0;
	int dataId = 0;
	int z = 0;
	int remained = 0;	  // 노드 그룹 내에서 아직 처리할 데이터가 남아있는 노드가 몇개인지 저장함
	
	paillier_ciphertext_t** tmp = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	for( i = 0 ; i < dim; i++ ){
		tmp[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(tmp[i]->c);
	}

	paillier_ciphertext_t*** cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(*NumNodeGroup)*FanOut);
		
	for( i = 0 ; i < *NumNodeGroup*FanOut ; i++ ){
		cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j ++ ){
			cand[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cand[i][j]->c);
		}
	}

	float progress = 0.1;
	
	// 노드 그룹 별로 데이터 추출을 통해, 질의 영역을 포함하는 노드 내 데이터 추출
	for(i=0; i<*NumNodeGroup; i++) {
		remained = node_group[i][0];	 // 해당 노드 그룹 내에서 아직 처리할 데이터가 남아있는 노드가 몇개인지 저장함

		for(j=0; j<FanOut; j++) {
			if(remained == 0)	// 해당 노드 그룹 내에서 더 이상 처리할 데이터가 없다면, 다음 노드로 넘어감
				break;

			if( (i*FanOut + j) / (float)(*NumNodeGroup*FanOut) >= progress)	{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			for(m=1; m<=node_group[i][0]; m++) {	 // 0번지에 해당 노드 그룹에 몇개의 노?弱?있는?側?저장되어 있음
				nodeId = node_group[i][m];	 // 노드 그룹에서 노드 ID를 하나씩 꺼냄
				//printf("selected node ID : %d\n", nodeId);

				if(node[nodeId].NumData >= j+1) 
				{
					dataId = node[nodeId].indata_id[j];	// 해당 노드에 저장된 데이터 ID를 하나씩 꺼냄
					//printf("selected data ID : %d\n", dataId);

					for(z=0; z<dim; z++) {
						if(m == 1) {
							//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, pubkey, prvkey,data[dataId][z]), paillier_dec(0, pubkey, prvkey, alpha[nodeId]));
							//printf("cnt : %d,  dim : %d\n", *cnt, z);
							cand[*cnt][z] = SM_p1(data[dataId][z], alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
						} else {
							//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, pubkey, prvkey,data[dataId][z]), paillier_dec(0, pubkey, prvkey, alpha[nodeId]));
							//printf("cnt : %d,  dim : %d\n", *cnt, z);
							tmp[z] = SM_p1(data[dataId][z], alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, tmp[z]));
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
							//gmp_printf("cnt : %d, (%Zd)\n", *cnt, paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
						}
					}
					//printf("\n");

					if(node[nodeId].NumData == j+1)		// 해당 노드가 마지막 데이터를 처리한다면, remained를 1 감소시킴
						remained--;
				}
				else {		// 해당 노드에는 데이터가 없지만, 동일 노드 그룹 내 다른 노드에는 아직 처리할 데이터가 있는 경우를 핸들링
					for(z=0; z<dim; z++) {
						if(m == 1) {
							//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, pubkey, prvkey, ciper_MAX), paillier_dec(0, pubkey, prvkey, alpha[nodeId]));
							//printf("cnt : %d,  dim : %d\n", *cnt, z);
							cand[*cnt][z] = SM_p1(ciper_MAX, alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
						} else {
							//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, pubkey, prvkey, ciper_MAX), paillier_dec(0, pubkey, prvkey, alpha[nodeId]));
							//printf("cnt : %d,  dim : %d\n", *cnt, z);
							tmp[z] = SM_p1(ciper_MAX, alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, tmp[z]));
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
							//gmp_printf("(%Zd)\n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
						}
					}
					//printf("\n");
				}
			}
			(*cnt)++;	// 노드 그룹의 노드들을 한바퀴 돌고나면, 데이터 하나가 완성됨
		}
	}
	//printf("\n");

	/*
	for(i=0; i<*cnt; i++) {
		gmp_printf("%dth data -> coord : ", i);
		for(j=0; j<dim; j++) {
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
		}
		printf("\n");
	}
	*/

	return cand;
}



// Proposed skNN with secure Index
int** protocol::SkNN_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode) {
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;

	bool verify_flag = false;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_ciphertext_t*** cand;
	paillier_ciphertext_t*** temp_cand ;
	paillier_ciphertext_t* temp_coord1 ;
	paillier_ciphertext_t* temp_coord2 ;
	paillier_ciphertext_t* temp_coord3 ;
	
	paillier_ciphertext_t*** ciper_nodedist_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumNode);
	paillier_ciphertext_t** temp_nodedist_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t*** ciper_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** ciper_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t**** ciper_nodeLL_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t**** ciper_nodeRR_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t** shortestPoint = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	for(i=0; i<dim; i++){
		shortestPoint[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(shortestPoint[i]->c);
	}

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* ciper_binary;
	paillier_ciphertext_t* ciper_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_dist	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_dist = paillier_create_enc_zero();
	paillier_ciphertext_t* kth_dist = 0;
	paillier_ciphertext_t** kth_dist_bit = 0;

	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*3);
	paillier_ciphertext_t* temp_alpha;
	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	for(i=0; i<NumNode; i++){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);

		ciper_nodedist_bit[i] = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*(size+1));
		for( j = 0 ; j < size+1 ; j++ ){ 
			ciper_nodedist_bit[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(ciper_nodedist_bit[i][j]->c);
		}

		ciper_nodeLL_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
		ciper_nodeRR_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	}

	for(i = 0 ; i < size+1 ; i++ ){ 
		temp_nodedist_bit[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(temp_nodedist_bit[i]->c);
	}

	for( i = 0 ; i < k ; i++ ){
		ciper_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim; j ++ ){
			ciper_result[i][j] = paillier_create_enc_zero();
		}
	}
		
	// query 비트 변환 수행
	for(i=0; i<dim; i++) {
		ciper_qLL_bit[i] = SBD_for_SRO(q[i], 0);		// query LL bound 변환
		ciper_qRR_bit[i] = SBD_for_SRO(q[i], 1);		// query RR bound 변환
	}

	/* 
	printf("Query LL bound\n");
	for(i=0; i<dim; i++) {
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qLL_bit[i][j]));
		}
		printf("\n");
	}
	printf("Query RR bound\n");
	for(i=0; i<dim; i++) {	
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qRR_bit[i][j]));
		}
		printf("\n");
	}
	*/

	printf("\n=== Node SBD start ===\n");

	startTime = clock();

	// 각 노드 비트 변환 수행 
	for(i=0; i<NumNode; i++) {	
		for(j=0; j<dim; j++) {
			ciper_nodeLL_bit[i][j] = SBD_for_SRO(node[i].LL[j], 0);				// node LL bound 변환
			ciper_nodeRR_bit[i][j] = SBD_for_SRO(node[i].RR[j], 1);	 			// node RR bound 변환
		}
		/*		
		printf("%dth node LL bound\n", i);
		for(j=0; j<dim; j++) {
			for(m=0; m<size+1; m++) {
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeLL_bit[i][j][m]));
			}
			printf("\n");
		}
		printf("%dth node RR bound\n", i);
		for(j=0; j<dim; j++) {
			for(m=0; m<size+1; m++) {
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeRR_bit[i][j][m]));
			}
			printf("\n");
		}
		*/
	}
	endTime = clock();
	node_SBD_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("node SBD time: %f\n", node_SBD_time);

	float progress;
	
	while(1) {
		progress = 0.1;
		startTime = clock();
		
		if(!verify_flag)
			printf("\n=== Node SRO start  ===\n");
		else
			printf("\n=== Node expansion start  ===\n");
		
		for(i=0; i<NumNode; i++) {	
			if(i/(float)NumNode >= progress)
			{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			if(!verify_flag)	{	//  검증 단계에서는 수행하지 않음
				alpha[i] = SRO(ciper_qLL_bit, ciper_qRR_bit, ciper_nodeLL_bit[i], ciper_nodeRR_bit[i]);		
				//gmp_printf("alpha : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, alpha[i]), __LINE__);
			}
			else {	// 검증 단계에서 k번째 거리보다 가까이 존재하는 노드를 찾는 과정
				// 각 노드에서 질의와의 최단 점을 찾음
				for(j=0; j<dim; j++) {	
					psi[0] = SCMP(ciper_qLL_bit[j], ciper_nodeLL_bit[i][j]);		// 프사이 계산 (1이면 q의 좌표가 노드의 LL bound보다 크고, 0이면 q의 좌표가 작음)
					//gmp_printf("psi0 : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, psi[0]), __LINE__);
					psi[1] = SCMP(ciper_qLL_bit[j], ciper_nodeRR_bit[i][j]);	// 프사이 계산 (1이면 q의 좌표가 노드의 UR bound보다 작고, 0이면 q의 좌표가 큼)
					//gmp_printf("psi1 : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, psi[1]), __LINE__);
					psi[2] = SBXOR(psi[0], psi[1]);	// psi[2]가 1이면, 해당 차원에서의 질의 좌표가 노드의 LL bound 와 UR bound 사이에 속함을 의미
					//gmp_printf("psi2 : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, psi[2]), __LINE__);

					temp_coord1 = SM_p1(psi[2], q[j]);		// psi[2]가 1이면, 질의에서의 수선의 발이 최단거리 이므로, 질의의 좌표를 살려야 함
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord1), __LINE__);
					temp_coord2 = SM_p1(psi[0], node[i].LL[j]);		// psi[0]이 0이면, LL bound의 좌표를 살림
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord2), __LINE__);
					temp_coord3 = SM_p1(SBN(psi[0]), node[i].RR[j]);		// psi[0]이 1이면, RR bound의 좌표를 살림
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord3), __LINE__);
					paillier_mul(pubkey, temp_coord3, temp_coord2, temp_coord3);
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord3), __LINE__);
					temp_coord3 = SM_p1(SBN(psi[2]), temp_coord3);
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord3), __LINE__);

					paillier_mul(pubkey, shortestPoint[j], temp_coord1, temp_coord3);
					//gmp_printf("coord: %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, shortestPoint[j]), __LINE__);
				}

				/*
				for(j=0; j<dim; j++) {	
					gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, shortestPoint[j]), __LINE__);
				}
				printf("\n");
				*/
				
				//gmp_printf("--------------------------------SSEDm %Zd ", paillier_dec(0, pubkey, prvkey, SSEDm(q, shortestPoint, dim)) , __LINE__);
				temp_nodedist_bit = SBD_for_SRO(SSEDm(q, shortestPoint, dim), 0);
				for(j=0; j<size+1; j++) {
					ciper_nodedist_bit[i][j] = SBOR(ciper_nodedist_bit[i][j], temp_nodedist_bit[j]);
				}

				alpha[i] = SCMP(ciper_nodedist_bit[i], kth_dist_bit);	 // 노드까지의 최단 거리를 계산해서 k번째 결과의 거리와 비교. 노드까지의 거리가 더 작으면 1 셋팅

				
				// 비트 변환된 노드 정보를 변경했던 버전 (노드까지의 거리 계산 시, size(L)의 범위가 벗어나서 변경)
				//	alpha[i] = SCMP(SBD_for_SRO(SSEDm(q, shortestPoint, dim), 0), kth_dist_bit);	 // 노드까지의 최단 거리를 계산해서 k번째 결과의 거리와 비교. 노드까지의 거리가 더 작으면 1 셋팅
			}


			if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
				printf("\n%dth node overlaps the query region.\n", i);
			
			if(!verify_flag)	{	//  검증 단계에서는 수행하지 않음
				// 이미 검색이 완료된 노드가 재 탐색되는 것을 방지하기 위해, bound를 MAX로 변환
				for(m=0; m<size+1; m++) {
					ciper_nodedist_bit[i][m] = alpha[i];
				}

				/*
				// 비트 변환된 노드 정보를 변경했던 버전 (노드까지의 거리 계산 시, size(L)의 범위가 벗어나서 변경)
				for(j=0; j<dim; j++) {
					for(m=0; m<size+1; m++) {
						ciper_nodeLL_bit[i][j][m] = SBOR(ciper_nodeLL_bit[i][j][m], alpha[i]);
						ciper_nodeRR_bit[i][j][m] = SBOR(ciper_nodeRR_bit[i][j][m], alpha[i]);
					}
				}
				*/
			}


			/*			
			printf("%dth node LL bound\n", i);
			for(j=0; j<dim; j++) {
				for(m=0; m<size+1; m++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeLL_bit[i][j][m]));
				}
				printf("\n");
			}
			printf("%dth node RR bound\n", i);
			for(j=0; j<dim; j++) {
				for(m=0; m<size+1; m++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeRR_bit[i][j][m]));
				}
				printf("\n");
			}
			*/
		}
		printf("\n");


		if(!verify_flag)	{
			endTime = clock();
			node_SRO_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("SRO time : %f\n", node_SRO_time);
		}
		else {
			endTime = clock();
			node_expansion_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node expansion time : %f\n", node_expansion_time);
		}

		/*
		printf("alpha\n");
		for(i=0; i<NumNode; i++) {
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, alpha[i]));
		}
		printf("\n");
		*/

		cnt = 0;

		if(!verify_flag)	{	// 검증 단계가 아닐 시에는, cand에 저장
			startTime = clock();
			
			cand = sNodeRetrievalforkNN(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
		
			endTime = clock();
			data_extract_first_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_first_time);

		}
		else {		// 검증 단계일 시에는, 이전 결과와 cand를 합침
			startTime = clock();

			temp_cand	= sNodeRetrievalforkNN(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = clock();
			data_extract_second_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_second_time);

			if(cnt == 0) {	// 검증을 위해 추가 탐색이 필요한 노드가 없는 경우를 처리함
				break;
			}

			cand =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ ){
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = paillier_create_enc_zero();
				}
			}

			// 이전 k개의 결과를 저장
			for(i=0; i<k; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = ciper_result[i][j];
				}
			}

			// 새로 찾은 후보 결과를 저장
			for(i=0; i<cnt; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i+k][j] = temp_cand[i][j];
				}
			}

			cnt = cnt + k;
		}

		paillier_ciphertext_t** ciper_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t*** ciper_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);
		paillier_ciphertext_t** ciper_V=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t*** ciper_V2 = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);

		paillier_ciphertext_t** ciper_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t** ciper_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		
		for( i = 0 ; i < size; i++ ){
			ciper_Smin[i] = ciper_zero;
		}

		for( i = 0 ; i < cnt ; i++ ){
			ciper_distance[i] 	= ciper_zero;
			ciper_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
			ciper_V[i]	= ciper_zero;
			ciper_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		}

		printf("\n=== Data SSED & SBD start ===\n");
		startTime = clock();
		for( i = 0 ; i < cnt ; i++ ){
			ciper_distance[i] = SSEDm(q, cand[i], dim);
			ciper_SBD_distance[i] = SBD(ciper_distance[i]);
			
			/*
			paillier_print("dist : ", ciper_distance[i]);	
			for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_SBD_distance[i][j]));
			}
			printf("\n");
			*/
		}
		endTime = clock();
		data_SSED_SBD_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		printf("data SSED & SBD time : %f\n", (float)(endTime-startTime)/(CLOCKS_PER_SEC));

		for( s = 0 ; s < k ; s++ ){
			//printf("\n%dth sMINn start \n", s+1);
			startTime = clock();
			ciper_Smin = Smin_n(ciper_SBD_distance, cnt);	// bit로 표현된 암호화 min 거리 추출
			endTime = clock();
			gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);

			if(!verify_flag)	{
				sMINn_first_time += gap;
				printf("%dth sMINn time : %f\n", s+1, gap);
			}
			else	{
				sMINn_second_time += gap;
				printf("%dth sMINn time : %f\n", s+1, gap);
			}
		
			/*	
			for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smin[j]));
			}
			printf("\n");	
			*/

			// bit로 표현된 암호화 min 거리를 암호화 정수로 변환
			for( j = size ; j > 0 ; j-- ){
				t = (int)pow(2, j-1);
				ciper_binary = paillier_create_enc(t);
				ciper_binary = SM_p1(ciper_binary, ciper_Smin[size-j]);
				paillier_mul(pubkey, ciper_min, ciper_binary, ciper_min);
				paillier_freeciphertext(ciper_binary);
				//paillier_print("min : ", ciper_min);
			}
			paillier_print("min dist : ", ciper_min);

			printf("\n== recalculate query<->data distances ===\n");
			// 이전 iteration에서 min값으로 선택된 데이터의 거리가 secure하게 MAX로 변환되었기 때문에, 
			// 질의-데이터 간 거리 계산을 모두 다시 수행
			if( s != 0 ){
				for( i = 0 ; i < cnt; i++ ){
					for( j = size ; j > 0 ; j--){
						t = (int)pow(2, j-1);
						ciper_binary = paillier_create_enc(t);
						ciper_binary = SM_p1(ciper_binary, ciper_SBD_distance[i][size-j]);
						paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
					}
					ciper_distance[i] = ciper_dist;
					ciper_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}
			/*
			for( i = 0 ; i < cnt ; i++ ){
				gmp_printf("distance : %Zd", ciper_distance[i]);
			}
			*/
			
			// 질의-데이터 거리와 min 거리와의 차를 구함 (min 데이터의 경우에만 0으로 만들기 위함)
			for( i = 0 ; i < cnt ; i++ ){
				paillier_subtract(pubkey, temp_dist, ciper_distance[i], ciper_min);
				ciper_mid[i] = SM_p1(temp_dist, ciper_rand);
				//paillier_print("dist - dist : ",ciper_mid[i]);
			}

			ciper_V = SkNNm_sub(ciper_mid, cnt);
			
			/*
			for( i = 0 ; i < cnt ; i++ ){
				printf("%d : ", i);
				paillier_print("ciper_V : ", ciper_V[i]);
			}
			*/

			// min 데이터 추출
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0 ; j < dim; j++ ){
					//gmp_printf("cand : %Zd(%d line)\n", paillier_dec(0, pubkey, prvkey, cand[i][j]), __LINE__);
					ciper_V2[i][j] = SM_p1(ciper_V[i], cand[i][j]);
					if(i==0) {
						ciper_result[s][j] = ciper_V2[i][j];
					}
					else {
						paillier_mul(pubkey, ciper_result[s][j], ciper_V2[i][j], ciper_result[s][j]);
					}
				}
			}
/*
			printf("ciper_result : ");		
			for( j = 0 ; j < dim; j++ ){
				gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey,ciper_result[s][j]));
			}
			printf("\n");		
*/			
			// Data SBOR 수행
			startTime = clock();
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0; j < size ; j++ ){
					ciper_SBD_distance[i][j] = SBOR(ciper_V[i], ciper_SBD_distance[i][j]);
				}
			}
			endTime = clock();
			data_SBOR_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		
			ciper_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
		}

		if(!verify_flag)	{	
			kth_dist = SSEDm(q, ciper_result[k-1], dim); // k번째 결과까지의 거리
			//gmp_printf("%dth dist : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_dist));

			kth_dist_bit = SBD_for_SRO(kth_dist, 1);

			
			/*
			// 비트 변환된 노드 정보를 변경했던 버전 (노드까지의 거리 계산 시, size(L)의 범위가 벗어나서 변경)
			for( i = 0 ; i < NumNode; i++ ){
				for( j = 0 ; j < dim; j++ ){
					for( m = size ; m > 0 ; m--){
						t = (int)pow(2, m-1);
						ciper_binary = paillier_create_enc(t);
						ciper_binary = SM_p1(ciper_binary, ciper_nodeLL_bit[i][j][size-m]);
						paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
					}
					node[i].LL[j] = ciper_dist;
					ciper_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);

					for( m = size ; m > 0 ; m--){
						t = (int)pow(2, m-1);
						ciper_binary = paillier_create_enc(t);
						ciper_binary = SM_p1(ciper_binary, ciper_nodeRR_bit[i][j][size-m]);
						paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
					}
					node[i].RR[j] = ciper_dist;
					ciper_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}
			*/


			/*
			for(i=0; i<dim; i++) {
				paillier_subtract(pubkey, temp_dist, q[i], kth_dist);	
				ciper_qLL_bit[i] = SBD_for_SRO(temp_dist, 0);		// query LL bound 변환
				gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, temp_dist));

				// LL bound의 경우, 음수 값이 되는 경우를 처리 (0으로 변경)
				temp_alpha = Smin_for_alpha(ciper_qLL_bit[i], SBD(paillier_create_enc_zero()));
				gmp_printf("alpha : %Zd(%d line)\n", paillier_dec(0, pubkey, prvkey, temp_alpha), __LINE__);
				for(j=0; j<size+1; j++) {
					ciper_qLL_bit[i][j] = SM_p1(ciper_qLL_bit[i][j], temp_alpha);	
				}

				paillier_mul(pubkey, temp_dist, q[i], kth_dist);	
				ciper_qRR_bit[i] = SBD_for_SRO(temp_dist, 1);		// query RR bound 변환
				gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, temp_dist));
			}
		*/

			/*
			printf("Query LL bound\n");
			for(i=0; i<dim; i++) {
				for(j=0; j<size+1; j++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qLL_bit[i][j]));
				}
				printf("\n");
			}
			printf("Query RR bound\n");
			for(i=0; i<dim; i++) {	
				for(j=0; j<size+1; j++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qRR_bit[i][j]));
				}
				printf("\n");
			}
			*/
		}

		if(!verify_flag)	{
			verify_flag = true;		// 검증을 수행하기 위해 flag를 true로 변경
		}
		else
			break;	// 검증을 끝냄
	}
		

	// user(Bob)에게 결과 전송을 위해 random 값 삽입
	for(i=0; i<k; i++){
		//printf("%d final result : ", i);
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,ciper_result[i][j],ciper_result[i][j],ciper_rand);
			//gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, ciper_result[i][j]));
		}
		//printf("\n");
	}

	paillier_freeciphertext(temp_dist);	
	
	return SkNNm_Bob2(ciper_result, rand, k, dim);
}

// Proposed skNN with secure Index + SMSn
int** protocol::SkNN_G(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode) {
	printf("\n=== GSRO_SMSn_SkNNsi start ===\n");
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;
	
	char Blink[10];
	
	float progress;
	
	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;
	paillier_ciphertext_t*		MAX				= paillier_create_enc(pow(2, size)-1);
	paillier_ciphertext_t*		C_RAND			= paillier_create_enc(rand);
	paillier_ciphertext_t*		K_DIST = 0;	
	paillier_ciphertext_t***	cand;
	paillier_ciphertext_t***	temp_cand;

	paillier_ciphertext_t*		temp_coord1		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	paillier_ciphertext_t*		temp_coord2		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	paillier_ciphertext_t*		temp_coord3		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));

	paillier_ciphertext_t*		TMP_alpha		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));	
	paillier_ciphertext_t*		TMP_NodeDIST	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	paillier_ciphertext_t**		alpha			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	paillier_ciphertext_t**		NodeDIST		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	paillier_ciphertext_t**		psi				= (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*3);
	paillier_ciphertext_t**		shortestPoint	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t***	ciper_result	= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	
	mpz_init(temp_coord1->c);
	mpz_init(temp_coord2->c);
	mpz_init(temp_coord3->c);
	mpz_init(TMP_alpha->c);
	mpz_init(TMP_NodeDIST->c);
	for( i = 0 ; i < 3 ; i ++)
	{
		psi[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(psi[i]->c);
	}
	for( i = 0 ; i < k ; i ++ )
	{
		ciper_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j++ )
		{
			ciper_result[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(ciper_result[i][j]->c);
		}
	}
	for( i = 0 ; i < dim ; i++ ){
		shortestPoint[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(shortestPoint[i]->c);
	}
	for( i = 0 ; i < NumNode ; i ++ )
	{
		NodeDIST[i]	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		alpha[i]	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(NodeDIST[i]->c);
		mpz_init(alpha[i]->c);
	}


	while(1)
	{
		progress = 0.1;
		startTime = clock();

		if(!verify_flag)
		{	
			cout<< "!!!!!!!!GSRO Start!!!!!!!!"<<endl;
			{
				for( i = 0 ; i < NumNode ; i ++ )
				{
					alpha[i] = DP_GSRO(q, q, node[i].LL, node[i].RR);
					if(Print)
					{
						cout<< i+1 <<" alpha : ";
						gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
					}
				}
			}
		}else{
			cout<< "!!!!!!!!Seconde Node Expansion!!!!!!!!"<<endl;
			for( i  = 0 ; i < NumNode ; i++ )
			{	
				//cout<< "NumNode : "<< i+1 <<endl;
				for(j=0; j<dim; j++) {	
					psi[0]			= GSCMP(q[j], node[i].LL[j]);		// 프사이 계산 (1이면 q의 좌표가 노드의 LL bound보다 크고, 0이면 q의 좌표가 작음)
					psi[1]			= GSCMP(q[j], node[i].RR[j]);	// 프사이 계산 (1이면 q의 좌표가 노드의 UR bound보다 작고, 0이면 q의 좌표가 큼)
					psi[2]			= SBXOR(psi[0], psi[1]);	// psi[2]가 1이면, 해당 차원에서의 질의 좌표가 노드의 LL bound 와 UR bound 사이??속함을 의미
					//gmp_printf("psi : %Zd %Zd %Zd\n", paillier_dec(0, pubkey, prvkey, psi[0]), paillier_dec(0, pubkey, prvkey, psi[1]), paillier_dec(0, pubkey, prvkey, psi[2]));
					temp_coord1		= SM_p1(psi[2], q[j]);		// psi[2]가 1이면, 질의에서의 수선의 발이 최단거리 이므로, 질의의 좌표를 살려야 함
					temp_coord2		= SM_p1(psi[0], node[i].LL[j]);		// psi[0]이 0이면, LL bound의 좌표를 살림
					temp_coord3		= SM_p1(SBN(psi[0]), node[i].RR[j]);		// psi[0]이 1이면, RR bound의 좌표를 살림
					//gmp_printf("temp_coord : %Zd %Zd %Zd\n", paillier_dec(0, pubkey, prvkey, temp_coord1), paillier_dec(0, pubkey, prvkey, temp_coord2), paillier_dec(0, pubkey, prvkey, temp_coord3));
					paillier_mul(pubkey, temp_coord3, temp_coord2, temp_coord3);
					temp_coord3		= SM_p1(SBN(psi[2]), temp_coord3);
					paillier_mul(pubkey, shortestPoint[j], temp_coord1, temp_coord3);
				}
				TMP_NodeDIST		= DP_SSED(q, shortestPoint, dim);
				//gmp_printf("TMP_NodeDIST %Zd\n", paillier_dec(0, pubkey, prvkey, TMP_NodeDIST));
				paillier_subtract(pubkey, TMP_alpha, ciper_one, alpha[i]);
				//gmp_printf("TMP_alpha %Zd\n", paillier_dec(0, pubkey, prvkey, TMP_alpha));
				paillier_mul(pubkey, NodeDIST[i], SM_p1(alpha[i], MAX), SM_p1(TMP_alpha, TMP_NodeDIST));				
				//gmp_printf("NodeDIST %Zd\n", paillier_dec(0, pubkey, prvkey, NodeDIST[i]));
				alpha[i]			= GSCMP(NodeDIST[i], K_DIST);
				//gmp_printf("alpha : %Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
				
				if(Print)
				{
					cout<< i+1 <<"Node Expansion alpha : ";
					gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
				}
			}
		}
		if(Print)
		{
			cout<<"enter any keys\n";
			cin>>Blink;
		}
		if(!verify_flag)	
		{
			endTime = clock();
			node_SRO_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("GSRO time : %f\n", node_SRO_time);
		}
		else 
		{
			endTime = clock();
			node_expansion_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node expansion time : %f\n", node_expansion_time);
		}
		cnt = 0;

		if(!verify_flag)
		{	// 검??단계가 아닐 시에는, cand에 저장
			startTime = clock();
			cand = GSRO_sNodeRetrievalforkNN(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_first_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_first_time);
		}else{
			startTime = clock();
			cout << "second GSRO_sNodeRetrievalforkNN start" <<endl;
			temp_cand = GSRO_sNodeRetrievalforkNN(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = clock();
			data_extract_second_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_second_time);

			if(cnt == 0)
			{	// 검증을 위해 추가 탐색이 필요한 노드가 없는 경우를 처리함
				break;
			}

			cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ )
			{
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
				for( j = 0 ; j < dim; j ++ )
				{
					cand[i][j] = paillier_create_enc_zero();
				}
			}

			// 이전 k개의 결과를 저장
			//cout << "cand : result "<<endl;
			for(i=0; i<k; i++)
			{
				for( j = 0 ; j < dim; j ++ )
				{
					cand[i][j] = ciper_result[i][j];
				}
			}
			// 새로 찾은 후보 결과를 저장
			for(i=0; i<cnt; i++)
			{
				for( j = 0 ; j < dim; j ++ )
				{
					cand[i+k][j] = temp_cand[i][j];
				}
			}
			cnt = cnt + k;
		}

		if(Print)
		{
			cout << "cnt : "<<cnt<<endl;
			cout << "!!!!EXTRACT CAND LIST!!!!"<<endl;
			for( i = 0 ; i < cnt; i ++)
			{
				for( j = 0 ; j < dim ; j ++ )
				{
					gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, cand[i][j]));
				}
				cout<<endl;
			}
			cout<<"enter any keys\n";
			cin>>Blink;
		}

		int idx = 0;
		paillier_ciphertext_t*		MIN				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		paillier_ciphertext_t**		ORIGIN_DIST		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t**		DIST			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t**		V				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t**		DIST_MINUS_MIN	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t***	V2				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);

		mpz_init(MIN->c);
		for( i = 0 ; i < cnt ; i ++ )
		{
			DIST_MINUS_MIN[i]	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			DIST[i]				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			V[i]				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			V2[i]				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
			
			mpz_init(DIST_MINUS_MIN[i]->c);
			mpz_init(DIST[i]->c);
			mpz_init(V[i]->c);
			for( j = 0 ; j < dim ; j++ ){
				V2[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
				mpz_init(V2[i][j]->c);
			}
		}
		for( i = 0 ; i < cnt ; i ++)
		{
			DIST[i]			= DP_SSED(cand[i], q,dim);

			if(Print)
			{
				cout<<"!!!!====SSED====!!!"<<endl;
				cout<< i+1 << " : Q ";
				for( j = 0 ; j < dim ; j++ )
				{
					gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, q[j]));
				}
				cout<< " | Data : ";
				for( j = 0 ; j < dim ; j++ )
				{
					gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
				}
				gmp_printf("|\t\t %Zd \n", paillier_dec(0, pubkey, prvkey, DIST[i]));
			}
		}
		for(s = 0 ; s < k ; s ++ )
		{
			cout << s+1 <<" th KNN Start !!!!!!!!!!!!"<<endl;
			if(Print)
			{
				for( i = 0 ; i < cnt ; i ++)
				{
					cout<< i+1 <<" DIST : ";
					gmp_printf(" %Zd \n", paillier_dec(0, pubkey, prvkey, DIST[i]));
				}
			}
			
			startTime = clock();

			idx = SMSn(DIST, cnt);

			endTime = clock();
			if(!verify_flag) {
				sMINn_first_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			}
			else {
				sMINn_second_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			}
			

			MIN = DIST[idx];
			for( i = 0 ; i < cnt ; i++ )
			{
				paillier_subtract(pubkey, DIST_MINUS_MIN[i], DIST[i], MIN);
				DIST_MINUS_MIN[i] = SM_p1(DIST_MINUS_MIN[i], C_RAND);
				if(Print)
				{
					gmp_printf("DIST-MIN : %Zd \n", paillier_dec(0, pubkey, prvkey, DIST_MINUS_MIN[i]));
				}
				
			}
			V = SkNNm_sub(DIST_MINUS_MIN, cnt);
			for( i = 0 ; i < cnt ; i++ )
			{
				if(Print)
				{
					gmp_printf("V : %Zd \n", paillier_dec(0, pubkey, prvkey, V[i]));

				}
				
				paillier_subtract(pubkey, TMP_alpha, ciper_one, V[i]);
				paillier_mul(pubkey, DIST[i], SM_p1(V[i], MAX), SM_p1(TMP_alpha, DIST[i]));
				for( j = 0 ; j < dim; j++ )
				{
					V2[i][j] = SM_p1(V[i], cand[i][j]);
					if( i == 0 )
					{
						ciper_result[s][j] = V2[i][j];
					}
					else
					{
						paillier_mul(pubkey, ciper_result[s][j], V2[i][j], ciper_result[s][j]);
					}
				}
			}
		}
		if(Print)
		{
			for( i = 0 ; i < k ; i++ )
			{
				cout << "MIDDLE RESULT : ";
				for( j = 0 ; j < dim ; j++ )
				{
					gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_result[i][j]));
				}
				cout <<endl;
			}
		}
		if(!verify_flag)
		{	
			K_DIST = DP_SSED(q, ciper_result[k-1], dim); // k번째 결과까지의 거리
		}

		if(!verify_flag)
			verify_flag = true;
		else
			break;
	}

	for(i=0; i<k; i++){
		NodeDIST[i] = DP_SSED(ciper_result[i], q, dim);
		gmp_printf("%d th DIST : %Zd\n", i+1 , paillier_dec(0, pubkey, prvkey, NodeDIST[i]));
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,ciper_result[i][j],ciper_result[i][j],C_RAND);
		}
	}

	free(cand);
	free(temp_cand);
	free(MAX);
	free(C_RAND);
	free(TMP_alpha);
	free(TMP_NodeDIST);
	free(alpha);
	free(NodeDIST);
	free(psi);
	free(shortestPoint);

	return SkNNm_Bob2(ciper_result, rand, k, dim);
}


// Proposed skNN with secure Index + packing SSED
int** protocol::DP_SkNN_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode) {
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;

	bool verify_flag = false;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_ciphertext_t*** cand;
	paillier_ciphertext_t*** temp_cand ;
	paillier_ciphertext_t* temp_coord1 ;
	paillier_ciphertext_t* temp_coord2 ;
	paillier_ciphertext_t* temp_coord3 ;
	
	paillier_ciphertext_t*** ciper_nodedist_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumNode);
	paillier_ciphertext_t** temp_nodedist_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t*** ciper_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** ciper_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t**** ciper_nodeLL_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t**** ciper_nodeRR_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t** shortestPoint = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	for(i=0; i<dim; i++){
		shortestPoint[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(shortestPoint[i]->c);
	}

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* ciper_binary;
	paillier_ciphertext_t* ciper_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_dist	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_dist = paillier_create_enc_zero();
	paillier_ciphertext_t* kth_dist = 0;
	paillier_ciphertext_t** kth_dist_bit = 0;

	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*3);
	paillier_ciphertext_t* temp_alpha;
	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	for(i=0; i<NumNode; i++){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);

		ciper_nodedist_bit[i] = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*(size+1));
		for( j = 0 ; j < size+1 ; j++ ){ 
			ciper_nodedist_bit[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(ciper_nodedist_bit[i][j]->c);
		}

		ciper_nodeLL_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
		ciper_nodeRR_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	}

	for(i = 0 ; i < size+1 ; i++ ){ 
		temp_nodedist_bit[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(temp_nodedist_bit[i]->c);
	}

	for( i = 0 ; i < k ; i++ ){
		ciper_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim; j ++ ){
			ciper_result[i][j] = paillier_create_enc_zero();
		}
	}
		
	// query 비트 변환 수행
	for(i=0; i<dim; i++) {
		ciper_qLL_bit[i] = SBD_for_SRO(q[i], 0);		// query LL bound 변환
		ciper_qRR_bit[i] = SBD_for_SRO(q[i], 1);		// query RR bound 변환
	}

	/* 
	printf("Query LL bound\n");
	for(i=0; i<dim; i++) {
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qLL_bit[i][j]));
		}
		printf("\n");
	}
	printf("Query RR bound\n");
	for(i=0; i<dim; i++) {	
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qRR_bit[i][j]));
		}
		printf("\n");
	}
	*/

	printf("\n=== Node SBD start ===\n");

	startTime = clock();

	// 각 노드 비트 변환 수행 
	for(i=0; i<NumNode; i++) {	
		for(j=0; j<dim; j++) {
			ciper_nodeLL_bit[i][j] = SBD_for_SRO(node[i].LL[j], 0);				// node LL bound 변환
			ciper_nodeRR_bit[i][j] = SBD_for_SRO(node[i].RR[j], 1);	 			// node RR bound 변환
		}
		/*		
		printf("%dth node LL bound\n", i);
		for(j=0; j<dim; j++) {
			for(m=0; m<size+1; m++) {
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeLL_bit[i][j][m]));
			}
			printf("\n");
		}
		printf("%dth node RR bound\n", i);
		for(j=0; j<dim; j++) {
			for(m=0; m<size+1; m++) {
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeRR_bit[i][j][m]));
			}
			printf("\n");
		}
		*/
	}
	endTime = clock();
	node_SBD_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("node SBD time: %f\n", node_SBD_time);

	float progress;
	
	while(1) {
		progress = 0.1;
		startTime = clock();
		
		if(!verify_flag)
			printf("\n=== Node SRO start  ===\n");
		else
			printf("\n=== Node expansion start  ===\n");
		
		for(i=0; i<NumNode; i++) {	
			if(i/(float)NumNode >= progress)
			{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			if(!verify_flag)	{	//  검증 단계에서는 수행하지 않음
				alpha[i] = SRO(ciper_qLL_bit, ciper_qRR_bit, ciper_nodeLL_bit[i], ciper_nodeRR_bit[i]);		
				//gmp_printf("alpha : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, alpha[i]), __LINE__);
			}
			else {	// 검증 단계에서 k번째 거리보다 가까이 존재하는 노드를 찾는 과정
				// 각 노드에서 질의와의 최단 점을 찾음
				for(j=0; j<dim; j++) {	
					psi[0] = SCMP(ciper_qLL_bit[j], ciper_nodeLL_bit[i][j]);		// 프사이 계산 (1이면 q의 좌표가 노드의 LL bound보다 크고, 0이면 q의 좌표가 작음)
					//gmp_printf("psi0 : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, psi[0]), __LINE__);
					psi[1] = SCMP(ciper_qLL_bit[j], ciper_nodeRR_bit[i][j]);	// 프사이 계산 (1이면 q의 좌표가 노드의 UR bound보다 작고, 0이면 q의 좌표가 큼)
					//gmp_printf("psi1 : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, psi[1]), __LINE__);
					psi[2] = SBXOR(psi[0], psi[1]);	// psi[2]가 1이면, 해당 차원에서의 질의 좌표가 노드의 LL bound 와 UR bound 사이에 속함을 의미
					//gmp_printf("psi2 : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, psi[2]), __LINE__);

					temp_coord1 = SM_p1(psi[2], q[j]);		// psi[2]가 1이면, 질의에서의 수선의 발이 최단거리 이므로, 질의의 좌표를 살려야 함
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord1), __LINE__);
					temp_coord2 = SM_p1(psi[0], node[i].LL[j]);		// psi[0]이 0이면, LL bound의 좌표를 살림
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord2), __LINE__);
					temp_coord3 = SM_p1(SBN(psi[0]), node[i].RR[j]);		// psi[0]이 1이면, RR bound의 좌표를 살림
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord3), __LINE__);
					paillier_mul(pubkey, temp_coord3, temp_coord2, temp_coord3);
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord3), __LINE__);
					temp_coord3 = SM_p1(SBN(psi[2]), temp_coord3);
					//gmp_printf("temp : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, temp_coord3), __LINE__);

					paillier_mul(pubkey, shortestPoint[j], temp_coord1, temp_coord3);
					//gmp_printf("coord: %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, shortestPoint[j]), __LINE__);
				}

				/*
				for(j=0; j<dim; j++) {	
					gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, shortestPoint[j]), __LINE__);
				}
				printf("\n");
				*/
				
				//gmp_printf("--------------------------------SSEDm %Zd ", paillier_dec(0, pubkey, prvkey, SSEDm(q, shortestPoint, dim)) , __LINE__);
				temp_nodedist_bit = SBD_for_SRO(SSEDm(q, shortestPoint, dim), 1);
				for(j=0; j<size+1; j++) {
					ciper_nodedist_bit[i][j] = SBOR(ciper_nodedist_bit[i][j], temp_nodedist_bit[j]);
				}

				alpha[i] = SCMP(ciper_nodedist_bit[i], kth_dist_bit);	 // 노드까지의 최단 거리를 계산해서 k번째 결과의 거리와 비교. 노드까지의 거리가 더 작으면 1 셋팅

				
				// 비트 변환된 노드 정보를 변경했던 버전 (노드까지의 거리 계산 시, size(L)의 범위가 벗어나서 변경)
				//	alpha[i] = SCMP(SBD_for_SRO(SSEDm(q, shortestPoint, dim), 0), kth_dist_bit);	 // 노드까지의 최단 거리를 계산해서 k번째 결과의 거리와 비교. 노드까지의 거리가 더 작으면 1 셋팅
			}


			if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
				printf("\n%dth node overlaps the query region.\n", i);
			
			if(!verify_flag)	{	//  검증 단계에서는 수행하지 않음
				// 이미 검색이 완료된 노드가 재 탐색되는 것을 방지하기 위해, bound를 MAX로 변환
				for(m=0; m<size+1; m++) {
					ciper_nodedist_bit[i][m] = alpha[i];
				}

				/*
				// 비트 변환된 노드 정보를 변경했던 버전 (노드까지의 거리 계산 시, size(L)의 범위가 벗어나서 변경)
				for(j=0; j<dim; j++) {
					for(m=0; m<size+1; m++) {
						ciper_nodeLL_bit[i][j][m] = SBOR(ciper_nodeLL_bit[i][j][m], alpha[i]);
						ciper_nodeRR_bit[i][j][m] = SBOR(ciper_nodeRR_bit[i][j][m], alpha[i]);
					}
				}
				*/
			}


			/*			
			printf("%dth node LL bound\n", i);
			for(j=0; j<dim; j++) {
				for(m=0; m<size+1; m++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeLL_bit[i][j][m]));
				}
				printf("\n");
			}
			printf("%dth node RR bound\n", i);
			for(j=0; j<dim; j++) {
				for(m=0; m<size+1; m++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeRR_bit[i][j][m]));
				}
				printf("\n");
			}
			*/
		}
		printf("\n");


		if(!verify_flag)	{
			endTime = clock();
			node_SRO_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("SRO time : %f\n", node_SRO_time);
		}
		else {
			endTime = clock();
			node_expansion_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node expansion time : %f\n", node_expansion_time);
		}

		/*
		printf("alpha\n");
		for(i=0; i<NumNode; i++) {
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, alpha[i]));
		}
		printf("\n");
		*/

		cnt = 0;

		if(!verify_flag)	{	// 검증 단계가 아닐 시에는, cand에 저장
			startTime = clock();
			
			cand = sNodeRetrievalforkNN(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
		
			endTime = clock();
			data_extract_first_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_first_time);

		}
		else {		// 검증 단계일 시에는, 이전 결과와 cand를 합침
			startTime = clock();

			temp_cand	= sNodeRetrievalforkNN(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = clock();
			data_extract_second_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_second_time);

			if(cnt == 0) {	// 검증을 위해 추가 탐색이 필요한 노드가 없는 경우를 처리함
				break;
			}

			cand =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ ){
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = paillier_create_enc_zero();
				}
			}

			// 이전 k개의 결과를 저장
			for(i=0; i<k; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = ciper_result[i][j];
				}
			}

			// 새로 찾은 후보 결과를 저장
			for(i=0; i<cnt; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i+k][j] = temp_cand[i][j];
				}
			}

			cnt = cnt + k;
		}

		paillier_ciphertext_t** ciper_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t*** ciper_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);
		paillier_ciphertext_t** ciper_V=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t*** ciper_V2 = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);

		paillier_ciphertext_t** ciper_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t** ciper_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		
		for( i = 0 ; i < size; i++ ){
			ciper_Smin[i] = ciper_zero;
		}

		for( i = 0 ; i < cnt ; i++ ){
			ciper_distance[i] 	= ciper_zero;
			ciper_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
			ciper_V[i]	= ciper_zero;
			ciper_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		}

		printf("\n=== Data SSED & SBD start ===\n");
		startTime = clock();
		for( i = 0 ; i < cnt ; i++ ){
			ciper_distance[i] = DP_SSED(q, cand[i], dim);
			ciper_SBD_distance[i] = SBD(ciper_distance[i]);
			
			/*
			paillier_print("dist : ", ciper_distance[i]);	
			for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_SBD_distance[i][j]));
			}
			printf("\n");
			*/
		}
		endTime = clock();
		data_SSED_SBD_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		printf("data SSED & SBD time : %f\n", (float)(endTime-startTime)/(CLOCKS_PER_SEC));

		for( s = 0 ; s < k ; s++ ){
			//printf("\n%dth sMINn start \n", s+1);
			startTime = clock();
			ciper_Smin = Smin_n(ciper_SBD_distance, cnt);	// bit로 표현된 암호화 min 거리 추출
			endTime = clock();
			gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);

			if(!verify_flag)	{
				sMINn_first_time += gap;
				printf("%dth sMINn time : %f\n", s+1, gap);
			}
			else	{
				sMINn_second_time += gap;
				printf("%dth sMINn time : %f\n", s+1, gap);
			}
		
			/*	
			for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smin[j]));
			}
			printf("\n");	
			*/

			// bit로 표현된 암호화 min 거리를 암호화 정수로 변환
			for( j = size ; j > 0 ; j-- ){
				t = (int)pow(2, j-1);
				ciper_binary = paillier_create_enc(t);
				ciper_binary = SM_p1(ciper_binary, ciper_Smin[size-j]);
				paillier_mul(pubkey, ciper_min, ciper_binary, ciper_min);
				paillier_freeciphertext(ciper_binary);
				//paillier_print("min : ", ciper_min);
			}
			//paillier_print("min dist : ", ciper_min);

			printf("\n== recalculate query<->data distances ===\n");
			// 이전 iteration에서 min값으로 선택된 데이터의 거리가 secure하게 MAX로 변환되었기 때문에, 
			// 질의-데이터 간 거리 계산을 모두 다시 수행
			if( s != 0 ){
				for( i = 0 ; i < cnt; i++ ){
					for( j = size ; j > 0 ; j--){
						t = (int)pow(2, j-1);
						ciper_binary = paillier_create_enc(t);
						ciper_binary = SM_p1(ciper_binary, ciper_SBD_distance[i][size-j]);
						paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
					}
					ciper_distance[i] = ciper_dist;
					ciper_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}
			/*
			for( i = 0 ; i < cnt ; i++ ){
				gmp_printf("distance : %Zd", ciper_distance[i]);
			}
			*/
			
			// 질의-데이터 거리와 min 거리와의 차를 구함 (min 데이터의 경우에만 0으로 만들기 위함)
			for( i = 0 ; i < cnt ; i++ ){
				paillier_subtract(pubkey, temp_dist, ciper_distance[i], ciper_min);
				ciper_mid[i] = SM_p1(temp_dist, ciper_rand);
				//paillier_print("dist - dist : ",ciper_mid[i]);
			}

			ciper_V = SkNNm_sub(ciper_mid, cnt);
			
			/*
			for( i = 0 ; i < cnt ; i++ ){
				printf("%d : ", i);
				paillier_print("ciper_V : ", ciper_V[i]);
			}
			*/

			// min 데이터 추출
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0 ; j < dim; j++ ){
					//gmp_printf("cand : %Zd(%d line)\n", paillier_dec(0, pubkey, prvkey, cand[i][j]), __LINE__);
					ciper_V2[i][j] = SM_p1(ciper_V[i], cand[i][j]);
					if(i==0) {
						ciper_result[s][j] = ciper_V2[i][j];
					}
					else {
						paillier_mul(pubkey, ciper_result[s][j], ciper_V2[i][j], ciper_result[s][j]);
					}
				}
			}
/*
			printf("ciper_result : ");		
			for( j = 0 ; j < dim; j++ ){
				gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey,ciper_result[s][j]));
			}
			printf("\n");		
			*/
			// Data SBOR 수행
			startTime = clock();
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0; j < size ; j++ ){
					ciper_SBD_distance[i][j] = SBOR(ciper_V[i], ciper_SBD_distance[i][j]);
				}
			}
			endTime = clock();
			data_SBOR_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		
			ciper_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
		}

		if(!verify_flag)	{	
			kth_dist = DP_SSED(q, ciper_result[k-1], dim); // k번째 결과까지의 거리
			//gmp_printf("%dth dist : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_dist));

			kth_dist_bit = SBD_for_SRO(kth_dist, 1);

			
			/*
			// 비트 변환된 노드 정보를 변경했던 버전 (노드까지의 거리 계산 시, size(L)의 범위가 벗어나서 변경)
			for( i = 0 ; i < NumNode; i++ ){
				for( j = 0 ; j < dim; j++ ){
					for( m = size ; m > 0 ; m--){
						t = (int)pow(2, m-1);
						ciper_binary = paillier_create_enc(t);
						ciper_binary = SM_p1(ciper_binary, ciper_nodeLL_bit[i][j][size-m]);
						paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
					}
					node[i].LL[j] = ciper_dist;
					ciper_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);

					for( m = size ; m > 0 ; m--){
						t = (int)pow(2, m-1);
						ciper_binary = paillier_create_enc(t);
						ciper_binary = SM_p1(ciper_binary, ciper_nodeRR_bit[i][j][size-m]);
						paillier_mul(pubkey, ciper_dist, ciper_binary, ciper_dist);
					}
					node[i].RR[j] = ciper_dist;
					ciper_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}
			*/


			/*
			for(i=0; i<dim; i++) {
				paillier_subtract(pubkey, temp_dist, q[i], kth_dist);	
				ciper_qLL_bit[i] = SBD_for_SRO(temp_dist, 0);		// query LL bound 변환
				gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, temp_dist));

				// LL bound의 경우, 음수 값이 되는 경우를 처리 (0으로 변경)
				temp_alpha = Smin_for_alpha(ciper_qLL_bit[i], SBD(paillier_create_enc_zero()));
				gmp_printf("alpha : %Zd(%d line)\n", paillier_dec(0, pubkey, prvkey, temp_alpha), __LINE__);
				for(j=0; j<size+1; j++) {
					ciper_qLL_bit[i][j] = SM_p1(ciper_qLL_bit[i][j], temp_alpha);	
				}

				paillier_mul(pubkey, temp_dist, q[i], kth_dist);	
				ciper_qRR_bit[i] = SBD_for_SRO(temp_dist, 1);		// query RR bound 변환
				gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, temp_dist));
			}
		*/

			/*
			printf("Query LL bound\n");
			for(i=0; i<dim; i++) {
				for(j=0; j<size+1; j++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qLL_bit[i][j]));
				}
				printf("\n");
			}
			printf("Query RR bound\n");
			for(i=0; i<dim; i++) {	
				for(j=0; j<size+1; j++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qRR_bit[i][j]));
				}
				printf("\n");
			}
			*/
		}

		if(!verify_flag)	{
			verify_flag = true;		// 검증을 수행하기 위해 flag를 true로 변경
		}
		else
			break;	// 검증을 끝냄
	}
		

	// user(Bob)에게 결과 전송을 위해 random 값 삽입
	for(i=0; i<k; i++){
		//printf("%d final result : ", i);
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,ciper_result[i][j],ciper_result[i][j],ciper_rand);
			//gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, ciper_result[i][j]));
		}
		//printf("\n");
	}

	paillier_freeciphertext(temp_dist);	
	
	return SkNNm_Bob2(ciper_result, rand, k, dim);
}
paillier_ciphertext_t*** protocol::GSRO_sNodeRetrievalforkNN(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha,  int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now GSRO_sNodeRetrievalforkNN starts =====\n");
	//printf("NumNode : %d\n", NumNode);
	int i=0, j=0, m=0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	int** node_group;
	
	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);
	printf("111set_num : %d\n", *NumNodeGroup);

	if(*NumNodeGroup == 0)
		return 0;

	for(i=0; i<*NumNodeGroup; i++) {
		printf("%dth Node Group : ", i);
		for(j=1; j<=node_group[i][0]; j++) {	// 0번지에 해당 노드 그룹에 몇개의 노드가 있는지가 저장되어 있음
			printf("%d ", node_group[i][j]);
		}
		printf("\n");
	}
	int nodeId = 0;
	int dataId = 0;
	int z = 0;
	int remained = 0;	  // 노드 그룹 내에서 아직 처리할 데이터가 남아있는 노드가 몇개인??저장함

	paillier_ciphertext_t** tmp = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	for( i = 0 ; i < dim; i++ ){
		tmp[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(tmp[i]->c);
	}

	paillier_ciphertext_t*** cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(*NumNodeGroup)*FanOut);
		
	for( i = 0 ; i < *NumNodeGroup*FanOut ; i++ ){
		cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j ++ ){
			cand[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cand[i][j]->c);
		}
	}
	float progress = 0.1;

	// 노드 그룹 별로 데이터 추출을 통해, 질의 영역을 포함하는 노드 내 데이터 추출
	for(i=0; i<*NumNodeGroup; i++) {
		remained = node_group[i][0];	 // 해당 노드 그룹 내에서 아직 처리할 데이터가 남아있는 노드가 몇개인지 저장함

		for(j=0; j<FanOut; j++) {
			//cout << "extract data s" << endl; 
			if(remained == 0)	// 해당 노드 그룹 내에서 더 이상 처리할 데이터가 없다면, 다음 노드로 넘어감
				break;

			if( (i*FanOut + j) / (float)(*NumNodeGroup*FanOut) >= progress)	{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			for(m=1; m<=node_group[i][0]; m++) {	 // 0번지에 해당 노드 그룹에 몇개의 노드가 있는지가 저장되어 있음
				nodeId = node_group[i][m];	 // 노드 그룹에서 노드 ID를 하나씩 꺼냄
				if (nodeId == -1){
					if(node[nodeId].NumData == j+1)		// 해당 노드가 마지막 데이터를 처리한다면, remained를 1 감소시킴
						remained--;
					continue;
				}
				
				//printf("selected node ID : %d\n", nodeId);

				if(node[nodeId].NumData >= j+1) 
				{
					dataId = node[nodeId].indata_id[j];	// 해당 노드에 저장된 데?謙?ID를 하나씩 꺼냄
					//printf("selected data ID : %d\n", dataId);

					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(data[dataId][z], alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
						} else {
							tmp[z] = SM_p1(data[dataId][z], alpha[nodeId]);  // 해당 데???ID의 실제 데이터에 접근
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}

					if(node[nodeId].NumData == j+1)		// 해당 노드가 마지막 데이터를 처리한다면, remained를 1 감소시킴
						remained--;
				}
				else {		// 해당 노드에는 데이터가 없지만, 동일 노드 그룹 내 다른 노드에는 아직 처리할 데이터가 있는 경우를 핸들링
					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(ciper_MAX, alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
						} else {
							tmp[z] = SM_p1(ciper_MAX, alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}
				}
			}
			(*cnt)++;	// 노드 그룹의 노드들을 한바퀴 돌고나면, 데이터 하나가 완성됨
			//cout << "extract data end" << endl;
		}
	}
	return cand;
}
int ** protocol::SSED_test(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, boundary* node, int k, int NumData, int NumNode){
	int i = 0, j = 0, s = 0;
	int ** original_data = (int **)malloc(sizeof(int*)*NumData);
	int * compute_data = (int *)malloc(sizeof(int)*NumData);
	int * result_data = (int *)malloc(sizeof(int)*NumData);
	int ** final_data = (int **)malloc(sizeof(int*)*k);
	Print = true;
	paillier_plaintext_t* Invert_int = paillier_plaintext_from_ui(0);
	paillier_ciphertext_t * Invert_Ciper = paillier_create_enc_zero();
	paillier_ciphertext_t ** SSED = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);
	/*
	for( i = 0 ; i < NumData ; i ++ )
	{
		SSED[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(SSED[i]->c);
	}
	*/
	int * Q = (int *)malloc(sizeof(int)*dim);
	for( i = 0 ; i < dim ; i ++ ){
		Invert_int = paillier_dec(0, pubkey, prvkey, query[i]);
		Q[i] = mpz_get_ui(Invert_int->m);
	}
	if(Print){
		cout << "Query : ";
		for( i = 0 ; i < dim ; i++ )
			cout << Q[i] << "  ";
		cout<<endl;
	}
	for( i = 0 ; i < NumData ; i ++){
		original_data[i] = (int *)malloc(sizeof(int)*dim);
		final_data[i] = (int *)malloc(sizeof(int)*dim);
		for( j = 0 ; j < dim ; j ++ ){
			Invert_int = paillier_dec(0, pubkey, prvkey, data[i][j]);
			original_data[i][j] = mpz_get_ui(Invert_int->m);
		}
	}
	if(Print){
		cout << "Original Data" << endl;
		for( i = 0 ; i < NumData ; i ++ ){
			cout << i+1<< " : ";
			for( j = 0 ; j < dim ; j++){
				cout << original_data[i][j] << "  ";
			}
			cout <<endl;
		}
	}
	/*
	for( i = 0 ; i < NumData; i ++ ){
		SSED[i] = DP_SSED(data[i],query,dim);
		Invert_int = paillier_dec(0, pubkey, prvkey, SSED[i]);
		compute_data[i] = mpz_get_ui(Invert_int->m);
	}
	*/
	/*
	SMSn(SSED, NumData);
	*/
	for( i = 0 ; i < NumData; i ++ ){
		result_data[i] = 0;
		for( j = 0 ; j < dim ; j++){
			result_data[i] = result_data[i] + (Q[j]-original_data[i][j]) * (Q[j]-original_data[i][j]) ;
		}
	}
	int swap;
	int m = 0;
	for(i=0;i<NumData-1;i++){
		for(j=i;j<NumData;j++){
			/*
			if(compute_data[i] > compute_data[j]){
				swap = compute_data[j];
				compute_data[j] = compute_data[i];
				compute_data[i] = swap;
			}
			*/
			if(result_data[i] > result_data[j]){
				swap = result_data[j];
				result_data[j] = result_data[i];
				result_data[i] = swap;
				for(m = 0 ; m < dim ; m++){
					swap = original_data[i][m];
					original_data[i][m] = original_data[j][m];
					original_data[j][m] = swap;
				}
			}
		}
	}
	
	int bre = 0;
	//cout << "Compute SSED"<<endl;
	for( i = 0 ; i < k ; i ++){
		cout << i+1 << " Q : ";
		for( j = 0 ; j < dim ; j++){
			cout << Q[j] << "  ";
		}cout << "  |Data :  ";
		for( j = 0 ; j < dim ; j++){
			final_data[i][j] = original_data[i][j];
			cout << original_data[i][j] << "  ";
		}
		cout << " ==? " << result_data[i] <<endl;
	}
	
	free(original_data);
	free(compute_data);
	free(Invert_int);
	free(Q);
	return final_data;
}