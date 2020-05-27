#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include<iostream>
#include "protocol.h"

using namespace std;

paillier_ciphertext_t* protocol::GSRO(paillier_ciphertext_t** qLL, paillier_ciphertext_t** qRR, paillier_ciphertext_t** nodeLL, paillier_ciphertext_t** nodeRR){
	//cout << "\nGSRO start" <<endl;
	int i = 0;
	int random = 0;
	int cnt = 0 ;
	int * Flag_array = (int *)malloc(sizeof(int)*(dim*2));
	int * iR1_array = (int *)malloc(sizeof(int)*(dim*2));
	int * iR2_array = (int *)malloc(sizeof(int)*(dim*2));
	paillier_ciphertext_t** F1_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));
	paillier_ciphertext_t** F2_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));

	paillier_ciphertext_t** R1_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));
	paillier_ciphertext_t** R2_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));
	
	for( i = 0 ; i < dim*2 ; i++ )
	{
		F1_array[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		F2_array[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		R1_array[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		R2_array[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(F1_array[i]->c);
		mpz_init(F2_array[i]->c);
		mpz_init(R1_array[i]->c);
		mpz_init(R2_array[i]->c);
		iR1_array[i] = 5;
		iR2_array[i] = 5;
		R1_array[i] = paillier_create_enc(iR1_array[i]);
		R2_array[i] = paillier_create_enc(iR2_array[i]);
		//	gmp_printf("%d : %Zd\t", i, F1_array[i]);		gmp_printf("%Zd\t", F2_array[i]);		gmp_printf("%Zd\t", R1_array[i]);		gmp_printf("%Zd\t", R2_array[i]);		cout << endl;
	}
	srand((unsigned int)time(NULL));
	for( i = 0 ; i < dim ; i ++ )
	{
		random = rand()%2;
		if (random==0)
		{
			Flag_array[i] = 0;
			paillier_mul(pub,F1_array[i],qLL[cnt],R1_array[i]);
			paillier_mul(pub,F2_array[i],nodeRR[cnt++],R2_array[i]);
		}
		else
		{
			Flag_array[i] = 1;
			paillier_mul(pub,F1_array[i],nodeRR[cnt],R2_array[i]);
			paillier_mul(pub,F2_array[i],qLL[cnt++],R1_array[i]);
		}
	}
	for( i = dim, cnt = 0 ; i < dim*2 ; i++ )
	{
		random = rand()%2;
		if( random == 0 )
		{
			Flag_array[i] = 0;
			paillier_mul(pub,F1_array[i],nodeLL[cnt],R1_array[i]);
			paillier_mul(pub,F2_array[i],qRR[cnt++],R2_array[i]);
		}
		else
		{
			Flag_array[i] = 1;
			paillier_mul(pub,F1_array[i],qRR[cnt],R2_array[i]);
			paillier_mul(pub,F2_array[i],nodeLL[cnt++],R1_array[i]);
		}
	}
	F1_array = GSRO_sub(F1_array,F2_array,iR1_array,iR2_array);

	for( i = 0 ; i < dim*2 ; i ++){
		//F1_array[i] = paillier_create_enc(Flag_array[i]);
		if(Flag_array[i] == 1 ){
			F1_array[i] = SBN(F1_array[i]);
		}
		//F1_array[i] = SBXOR(F1_array[i],F2_array[i]);
//	 	gmp_printf("%Zd\t",paillier_dec(0, pub, prv, F1_array[i]));
	}
	for( i = 0 ; i < dim*2 ; i ++){
		F1_array[0] = SM_p1(F1_array[0],F1_array[i]);
	}

	//gmp_printf("\n\nrap : %Zd\n", paillier_dec(0, pub, prv, F1_array[0]));
	return F1_array[0];
}

paillier_ciphertext_t** protocol::GSRO_sub(paillier_ciphertext_t** A, paillier_ciphertext_t** B, int * R1, int * R2){
	int i = 0;
	paillier_plaintext_t* plain;
	paillier_ciphertext_t** result = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));
	int * iresult = (int *)malloc(sizeof(int)*dim*2);
	int * iA = (int *)malloc(sizeof(int)*dim*2);
	int * iB = (int *)malloc(sizeof(int)*dim*2);
	for( i = 0 ; i < dim*2 ; i ++){
		plain = paillier_dec(0, pub, prv, A[i]);
		iA[i] = mpz_get_ui(plain->m);
		plain = paillier_dec(0, pub, prv, B[i]);
		iB[i] = mpz_get_ui(plain->m);
	}
	/*
	cout << endl;
	cout << endl;
	
	cout << " GSRO_sub " <<endl;
	for( i = 0 ; i < dim*2 ; i ++){
		cout << iA[i] << "\t";
	}
	cout << endl;
	for( i = 0 ; i < dim*2 ; i ++){
		cout << iB[i] << "\t";
	}
	cout << endl;
	for( i = 0 ; i < dim*2 ; i ++){
		cout << R1[i] << "\t";
	}
	cout << endl;
	for( i = 0 ; i < dim*2 ; i ++){
		cout << R2[i] << "\t";
	}
	cout << endl;
	*/
	for(i=0; i < dim*2; i++){
		iresult[i] = G_CMP(iB[i],R2[i],iA[i],R1[i]);
	/*	if( (iA[i]-R1[i]) < (iB[i]-R2[i]) ){
			iresult[i] = 1;
		}else{
			iresult[i] = 0;
		}
		cout << iresult[i] << "\t";
	*/	result[i] = paillier_create_enc(iresult[i]);
	}
	return result;
}

paillier_ciphertext_t* protocol::DP_GSRO(paillier_ciphertext_t** qLL, paillier_ciphertext_t** qRR, paillier_ciphertext_t** nodeLL, paillier_ciphertext_t** nodeRR){
	//cout << "DP_GSRO start" <<endl;
	int i = 0;
	int * Flag_array = (int *)malloc(sizeof(int)*(dim*2));
	int * R1 = (int *)malloc(sizeof(int)*(dim*2));
	int * R2 = (int *)malloc(sizeof(int)*(dim*2));
	for( i = 0 ; i < dim*2 ; i++){
		Flag_array[i] = 0;
		R1[i] = 5;
		R2[i] = 5;
	}
	paillier_ciphertext_t** F1_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2)); 
	paillier_ciphertext_t** F2_array; 

	paillier_plaintext_t* shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* p1_array = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* p2_array = paillier_plaintext_from_ui(0);
	paillier_plaintext_t** R1_array = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*(dim*2));
	paillier_plaintext_t** R2_array = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*(dim*2));
	paillier_plaintext_t* tmp= paillier_plaintext_from_ui(0);
	
	paillier_ciphertext_t* c1_array;
	paillier_ciphertext_t* c2_array;
	paillier_ciphertext_t* ciper_tmp = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));	
	mpz_init(ciper_tmp->c);

	paillier_ciphertext_t** temp_qLL = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** temp_qRR = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** temp_nodeLL = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** temp_nodeRR = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);

	for( i = 0 ; i < dim ; i ++){
		temp_qLL[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp_qRR[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp_nodeLL[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp_nodeRR[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(temp_qLL[i]->c);
		mpz_init(temp_qRR[i]->c);
		mpz_init(temp_nodeLL[i]->c);
		mpz_init(temp_nodeRR[i]->c);
	}

	// after shifting and add random value under Array P1
	for(i=0;i<dim*2;i++){
		mpz_ui_pow_ui(shift->m,2,DP*i);	
		R1_array[i] = paillier_plaintext_from_ui(R1[i]);
		mpz_mul(tmp->m,R1_array[i]->m,shift->m);
		mpz_add(p1_array->m,p1_array->m, tmp->m);
	}
	// after shifting and add random value under Array P2
	for(i=0;i<dim*2;i++){
		mpz_ui_pow_ui(shift->m,2,DP*i);	
		R2_array[i] = paillier_plaintext_from_ui(R2[i]);
		mpz_mul(tmp->m,R2_array[i]->m,shift->m);
		mpz_add(p2_array->m,p2_array->m, tmp->m);
	}
	
	// enc Array P1, enc Array P2
	c1_array = paillier_enc(0, pubkey, p1_array, paillier_get_rand_devurandom);
	c2_array = paillier_enc(0, pubkey, p2_array, paillier_get_rand_devurandom);

	int random = 0 ;
	int cnt = 0 ;	

	// 2qLL 2nodeLL
	for( i = 0 ; i < dim ; i++){
		mpz_ui_pow_ui(shift->m,2,1);
		paillier_exp(pubkey, temp_qLL[i], qLL[i], shift);
		paillier_exp(pubkey, temp_nodeLL[i], nodeLL[i], shift);
	}
	// 2qRR+1 2nodeRR+1
	for( i = 0 ; i < dim ; i++){
		paillier_exp(pubkey, temp_nodeRR[i], nodeRR[i], shift); 
		paillier_mul(pubkey, temp_nodeRR[i], temp_nodeRR[i], ciper_one); 
		
		paillier_exp(pubkey, temp_qRR[i], qRR[i], shift);
		paillier_mul(pubkey, temp_qRR[i], temp_qRR[i], ciper_one); 
	}

	// add to Array QueryLL, nodeRR value
	for( i = 0 ; i < dim ; i++ ){
		random = rand()%2;
		mpz_ui_pow_ui(shift->m,2,DP*i);
		if (random == 0){
			Flag_array[i] = 0;
			paillier_exp(pubkey, ciper_tmp, temp_qLL[cnt], shift); //multi
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
			paillier_exp(pubkey, ciper_tmp, temp_nodeRR[cnt++], shift); //multi
			paillier_mul(pubkey, c2_array, c2_array, ciper_tmp); //E_add			
		}else{
			Flag_array[i] = 1;
			paillier_exp(pubkey, ciper_tmp, temp_nodeRR[cnt], shift); //multi
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
			paillier_exp(pubkey, ciper_tmp, temp_qLL[cnt++], shift); //multi
			paillier_mul(pubkey, c2_array, c2_array, ciper_tmp); //E_add
		}
	}
	// add to Array QueryRR, nodeLL value	
	for( i = dim, cnt = 0; i < dim*2 ; i++){
		random = rand()%2;
		mpz_ui_pow_ui(shift->m,2,DP*i);
		if (random == 0){
			Flag_array[i] = 0;
			paillier_exp(pubkey, ciper_tmp, temp_nodeLL[cnt], shift); //multi
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
			paillier_exp(pubkey, ciper_tmp, temp_qRR[cnt++], shift); //multi
			paillier_mul(pubkey, c2_array, c2_array, ciper_tmp); //E_add			
		}else{
			Flag_array[i] = 1;
			paillier_exp(pubkey, ciper_tmp, temp_qRR[cnt], shift); //multi
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
			paillier_exp(pubkey, ciper_tmp, temp_nodeLL[cnt++], shift); //multi
			paillier_mul(pubkey, c2_array, c2_array, ciper_tmp); //E_add			
		}
	}

	F2_array = unDP_GSRO(c1_array, c2_array, R1, R2);
	
	
	for( i = 0 ; i < dim*2 ; i ++){
		F1_array[i] = paillier_create_enc(Flag_array[i]);
		F1_array[i] = SBXOR(F1_array[i],F2_array[(dim*2)-i-1]);
	}
	for( i = 0 ; i < dim*2 ; i ++){
		F1_array[0] = SM_p1(F1_array[0],F1_array[i]);
	}

	//gmp_printf("%Zd\n", paillier_dec(0, pub, prv, F1_array[0]));
	free(Flag_array);
	free(R1);
	free(R2);
	free(shift);
	free(p1_array);
	free(p2_array);
	free(tmp);
	free(ciper_tmp);
	for(i=0;i<dim*2;i++){
		free(R1_array[i]);
		free(R2_array[i]);
		free(F2_array[i]);
	}
	free(R1_array);
	free(R2_array);
	free(F2_array);
	for(i=0;i<dim;i++){
		free(temp_qLL[i]);
		free(temp_qRR[i]);
		free(temp_nodeLL[i]);
		free(temp_nodeRR[i]);
	}
	return F1_array[0];
}
paillier_ciphertext_t** protocol::unDP_GSRO(paillier_ciphertext_t* A, paillier_ciphertext_t* B, int * R1, int * R2){
	//cout << "unDP_GSRO start" <<endl;
	paillier_ciphertext_t** result	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));
	int * iresult = (int *)malloc(sizeof(int)*dim*2);
	int i = 0;
	paillier_plaintext_t* shift		= paillier_plaintext_from_ui(0);
	paillier_plaintext_t* tmp		= paillier_plaintext_from_ui(0);
	paillier_plaintext_t* A_array	= paillier_dec(0, pub, prv, A);
	paillier_plaintext_t* B_array	= paillier_dec(0, pub, prv, B);
	int * A_array_result			= (int *)malloc(sizeof(int)*(dim*2));
	int * B_array_result			= (int *)malloc(sizeof(int)*(dim*2));
	
	for( i = 0 ; i < dim*2 ; i++){
		mpz_ui_pow_ui(shift->m,2,((dim*2)-i-1)*DP);
		mpz_div(tmp->m,A_array->m,shift->m);
		A_array_result[i] = mpz_get_ui(tmp->m);
		mpz_mod(A_array->m,A_array->m,shift->m);
	}
	for( i = 0 ; i < dim*2 ; i++){
		mpz_ui_pow_ui(shift->m,2,((dim*2)-i-1)*DP);
		mpz_div(tmp->m,B_array->m,shift->m);
		B_array_result[i] = mpz_get_ui(tmp->m);
		mpz_mod(B_array->m,B_array->m,shift->m);
	}

	makeGate();
	for(i=0; i < dim*2; i++){
		iresult[i] = G_CMP(B_array_result[i],R2[i],A_array_result[i],R1[i]);
		result[i] = paillier_create_enc(iresult[i]);
	}
	free(shift);
	free(tmp);
	free(A_array);
	free(B_array);
	free(A_array_result);
	free(B_array_result);
	free(iresult);
	return result;
}	


paillier_ciphertext_t*** protocol::GSRO_sNodeRetrievalforRange(paillier_ciphertext_t*** data, paillier_ciphertext_t** ciper_qLL, paillier_ciphertext_t** ciper_qRR, boundary* node, int NumData, int NumNode, int* cnt, int* NumNodeGroup, int type)
{
	printf("\n===== Now GSRO_sNodeRetrievalforRange starts =====\n");
	printf("NumNode : %d\n", NumNode);
	int i=0, j=0, m=0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	int** node_group;

	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	for(i=0; i<NumNode; i++){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);
	}

	// 각 노드 GSRO 호출
	for(i=0; i<NumNode; i++) {	
		startTime = clock();
		alpha[i] = DP_GSRO(ciper_qLL, ciper_qRR, node[i].LL, node[i].RR);		
		endTime = clock();
		node_SRO_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pub, prv, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
			printf("%dth node overlaps the query region.\n", i);
	}

	printf("node SBD time: %f\n", node_SBD_time);
	printf("SRO time : %f\n", node_SRO_time);

	startTime = clock();
	
	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);

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
		
	for( i = 0 ; i < (*NumNodeGroup)*FanOut ; i++ ){
		cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j ++ ){
			cand[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cand[i][j]->c);
		}
	}


	// 노드 그룹 별로 데이터 추출을 통해, 질의 영역을 포함하는 노드 내 데이터 추출
	for(i=0; i<*NumNodeGroup; i++) {
		remained = node_group[i][0];	 // 해당 노드 그룹 내에서 아직 처리할 데이터가 남아있는 노드가 몇개인지 저장함

		for(j=0; j<FanOut; j++) {
			if(remained == 0)	// 해당 노드 그룹 내에서 더 이상 처리할 데이터가 없다면, 다음 노드로 넘어감
				break;

			for(m=1; m<=node_group[i][0]; m++) {	 // 0번지에 해당 노드 그룹에 몇개의 노드가 있는지가 저장되어 있음
				nodeId = node_group[i][m];	 // 노드 그룹에서 노드 ID를 하나씩 꺼냄
				if(node[nodeId].NumData >= j+1) 
				{
					dataId = node[nodeId].indata_id[j];	// 해당 노드에 저장된 데이터 ID를 하나씩 꺼냄
					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(data[dataId][z], alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
						} else {
							tmp[z] = SM_p1(data[dataId][z], alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
							paillier_mul(pub, cand[*cnt][z], cand[*cnt][z], tmp[z]);
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
							paillier_mul(pub, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}
				}
			}
			(*cnt)++;	// 노드 그룹의 노드들을 한바퀴 돌고나면, 데이터 하나가 완성됨
		}
	}

	endTime = clock();
	data_extract_first_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("Node Retrieval time : %f\n", data_extract_first_time);

	

	return cand;
}

int** protocol::sRange_G(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num)
{
	printf("\n===== Now sRangeG starts =====\n");
	//printf("NumNode : %d\n", NumNode);
	int i=0, j=0, m=0;
	int rand = 5;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	int NumNodeGroup = 0;

	paillier_ciphertext_t* ciper_rand = paillier_create_enc(rand);
	paillier_print("rand : ",ciper_rand);


	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수

	paillier_ciphertext_t*** cand ;

	cand = GSRO_sNodeRetrievalforRange(data, q.LL, q.RR, node, NumData, NumNode, &cnt, &NumNodeGroup, 1);

	

	totalNumOfRetrievedNodes += NumNodeGroup;

	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
	
	// cand 비트 변환 수행 및 SPE호출
	for(i=0; i<cnt; i++) {
		startTime = clock();
		alpha[i] = DP_GSRO(cand[i], cand[i], q.LL, q.RR);
		//alpha[i]=SPE(cand_bit,ciper_nodeLL_bit, ciper_nodeRR_bit);		
		endTime = clock();
		data_SPE_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		//printf("data SPE time : %f\n", (float)(endTime-startTime)/(CLOCKS_PER_SEC));
	
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,cand[i][j],cand[i][j],ciper_rand);
		}
	}

	

	paillier_ciphertext_t*** result;
	result = sRange_result(alpha, cand, cnt, NumNodeGroup, result_num);
	

	return FsRange_Bob(result, rand, (*result_num), dim);
}


paillier_ciphertext_t*** protocol::sRange_result(paillier_ciphertext_t** alpha, paillier_ciphertext_t*** cand, int cnt, int NumNodeGroup, int * result_num){
	int i,j;

	cout<<"sRange_result NumNodeGroup : "<<NumNodeGroup<<endl;
	cout<<"sRange_result FanOut : "<<FanOut<<endl;

	paillier_ciphertext_t*** result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumNodeGroup*FanOut);
	for(i=0; i<cnt; i++){
		if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pub, prv, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
		{
			result[(*result_num)]=cand[i]; //result에 삽입
			(*result_num)++;
		}
	}
	return result;
}