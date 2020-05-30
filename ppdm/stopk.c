#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>

#include "protocol.h"
#include <iostream>

using namespace std;

paillier_ciphertext_t*** protocol::sNodeRetrievalforTopk(paillier_ciphertext_t*** data, paillier_ciphertext_t*** ciper_qLL_bit, paillier_ciphertext_t*** ciper_qRR_bit, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now sNodeRetrievalforTopk starts =====\n");
	int i=0, j=0, m=0;

	int** node_group;

	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);
	if(*NumNodeGroup == 0)
		return 0;


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
		remained = node_group[i][0];	 // 해당 노드 그룹 내에서 아직 처리할 데이터가 남?팀獵?노드가 몇개인지 저장함
		for(j=0; j<FanOut; j++) {
			if(remained == 0)	// 해당 노드 그룹 내에서 더 이상 처리할 데이터가 없다면, 다음 노드로 넘어감
				break;

			if( (i*FanOut + j) / (float)(*NumNodeGroup*FanOut) >= progress)	{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

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
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}

					if(node[nodeId].NumData == j+1)		// 해당 노드가 마지막 데이터를 처리한다면, remained를 1 감소시킴
						remained--;
				}
				else {		// 해당 노드에는 데이터가 없지만, 동일 노드 그룹 내 다른 노드에는 아??처리할 데이터가 있는 경우를 핸들링
					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(cipher_zero, alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
						} else {
							tmp[z] = SM_p1(cipher_zero, alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}

				}
			}
			(*cnt)++;	// 노드 그룹의 노드들을 한바퀴 돌고나면, 데이터 하나가 완성됨
		}
	}
	return cand;
}

paillier_ciphertext_t*** protocol::GSRO_sNodeRetrievalforTopk(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha,  int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now GSRO_sNodeRetrievalforkNN starts =====\n");
	//printf("NumNode : %d\n", NumNode);
	int i=0, j=0, m=0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	int** node_group;
	
	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);
	printf("set_num : %d\n", *NumNodeGroup);

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
					dataId = node[nodeId].indata_id[j];	// 해당 노드에 저장된 데이터 ID를 하나씩 꺼냄
					//printf("selected data ID : %d\n", dataId);

					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(data[dataId][z], alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
						} else {
							tmp[z] = SM_p1(data[dataId][z], alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}

					if(node[nodeId].NumData == j+1)		// 해당 노드가 마지막 데이터를 처리한다면, remained를 1 감소시킴
						remained--;
				}
				else {		// 해당 노드에는 데이터가 없지만, 동일 노드 그룹 내 다른 노드에는 아직 처리할 데이터가 있는 경우를 핸들링
					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(cipher_zero, alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
						} else {
							tmp[z] = SM_p1(cipher_zero, alpha[nodeId]);  // 해당 데이터 ID의 실제 데이터에 접근
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


// 151024. added by KHJ.
int** protocol::STopk_Confirm(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode){
	cout<< "!!!TopK Authentication!!!"<<endl;
	int i = 0, j = 0, m = 0, swap, idx;
	int ** original_data = (int **)malloc(sizeof(int*)*NumData);
	int ** result_data = (int **)malloc(sizeof(int*)*NumData);
	int * compute_data = (int *)malloc(sizeof(int)*NumData);
	int * result = (int *)malloc(sizeof(int)*NumData);
	int * query = (int *)malloc(sizeof(int)*dim);
	int hint_int = 0;

	paillier_ciphertext_t** ciper_Smax		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t*** SBD_data		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	for( i = 0 ; i < NumData ; i++ )
	{
		SBD_data[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	}
	paillier_ciphertext_t * hint = paillier_create_enc_zero();
	paillier_ciphertext_t** coeff = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_plaintext_t* Invert_int = paillier_plaintext_from_ui(0);
	paillier_ciphertext_t ** C_result = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);

	hint = q[dim];
	Invert_int = paillier_dec(0, pubkey, prvkey, q[dim]);
	hint_int = mpz_get_ui(Invert_int->m);
	cout<<"Query Set"<<endl;
	for(i=0; i<dim; i++){
		coeff[i] = paillier_create_enc_zero();
		paillier_mul(pubkey, coeff[i], q[i], hint);
		Invert_int = paillier_dec(0, pubkey, prvkey, coeff[i]);
		query[i] = mpz_get_ui(Invert_int->m);
		query[i] = query[i] - hint_int;
		cout << query[i]<<" ";	
	}cout<<endl;
	for(i = 0 ; i < NumData; i++){
		original_data[i] = (int *)malloc(sizeof(int)*dim);
		result_data[i] = (int *)malloc(sizeof(int)*dim);
		original_data[i][j] = 0;
		result_data[i][j] = 0;
	}
	//cout <<"SBD"<<endl;
	for( i = 0 ; i < NumData; i++){
		C_result[i] = computeScore(data[i],q);
		//SBD_data[i] = SBD(C_result[i]);
		//Invert_int = paillier_dec(0, pubkey, prvkey, C_result[i]);
		//compute_data[i] = mpz_get_ui(Invert_int->m);
		for( j = 0 ; j < dim ; j++){
			Invert_int = paillier_dec(0, pubkey, prvkey, data[i][j]);
			original_data[i][j] = mpz_get_ui(Invert_int->m);
			result_data[i][j] = mpz_get_ui(Invert_int->m);
		}
	}
	idx = MAXn(C_result, NumData);

	gmp_printf("MAX Data : %Zd\n",  paillier_dec(0, pubkey, prvkey, C_result[idx]));	
	/*
	cout<<"=================MAX================="<<endl;
	for(j=0;j<size;j++)
	{
			gmp_printf("%Zd",  paillier_dec(0, pubkey, prvkey, ciper_Smax[j]));	
	}
	cout <<endl;
	*/
	for(i=0;i<NumData;i++){
		int sum = 0;
		for(j=0;j<dim;j++){
			sum = sum + query[j]*original_data[i][j];
		}
		result[i] = sum;
	}
	for(i=0;i<NumData-1;i++){
		for(j=i;j<NumData;j++){
			if(result[i] < result[j]){
				swap = result[j];
				result[j] = result[i];
				result[i] = swap;
				for(m = 0 ; m < dim ; m++){
					swap = result_data[i][m];
					result_data[i][m] = result_data[j][m];
					result_data[j][m] = swap;
				}
			}
			if(compute_data[i] < compute_data[j]){
				swap = compute_data[j];
				compute_data[j] = compute_data[i];
				compute_data[i] = swap;
				for(m = 0 ; m < dim ; m++){
					swap = original_data[i][m];
					original_data[i][m] = original_data[j][m];
					original_data[j][m] = swap;
				}
			}
		}
	}
	/*
	for(i = 0 ; i < 1000 ; i++ ){
		cout<< i <<" : ";
		for(j=0;j<dim;j++){
			cout<<original_data[i][j]<<" ";
		}
		cout << " | compute : " << compute_data[i] <<endl;
	}
	*/
	cout <<"============ Topk result ============="<<endl;
	for(i = 0 ; i < k;i++){
		cout<< i <<" : ";
		cout << "compute : " << result[i] <<" | ";
		for(j=0;j<dim;j++){
			cout<<result_data[i][j]<<" ";
		}
		cout<<endl;
	} 
	
	return result_data;
}

paillier_ciphertext_t* protocol::computeScore(paillier_ciphertext_t** cipher1, paillier_ciphertext_t** cipher2)
{
	// init variables
	int i = 0;
	paillier_ciphertext_t* score = paillier_create_enc_zero();
	paillier_ciphertext_t* tmp_score = 0; 

	for( i = 0 ; i < dim; i++ )
	{
		tmp_score = SM_p1(cipher1[i], cipher2[i]);	
		//gmp_printf("tmp_score : %Zd\n", paillier_dec(0, pubkey, prvkey, tmp_score));
		paillier_mul(pubkey, score, score, tmp_score);
		//gmp_printf("score : %Zd\n", paillier_dec(0, pubkey, prvkey, score));		
	}

	paillier_freeciphertext(tmp_score);

	return score;
}

// 150809. added by KHI.
paillier_ciphertext_t* protocol::computeScore2(paillier_ciphertext_t** cipher1, paillier_ciphertext_t** cipher2, paillier_ciphertext_t** coeff, paillier_ciphertext_t* hint)
{
	// init variables
	int i = 0;
	paillier_ciphertext_t* score = paillier_create_enc_zero();
	paillier_ciphertext_t* tmp_score = 0; 

	for( i = 0 ; i < dim; i++ )
	{
		tmp_score = SM_p1(coeff[i], cipher2[i]);	
		//gmp_printf("tmp_score : %Zd\n", paillier_dec(0, pubkey, prvkey, tmp_score));
		paillier_subtract(pubkey, tmp_score, tmp_score, SM_p1(hint, cipher2[i]));
		//gmp_printf("tmp_score : %Zd\n", paillier_dec(0, pubkey, prvkey, tmp_score));
		paillier_mul(pubkey, score, score, tmp_score);
		//gmp_printf("score : %Zd\n", paillier_dec(0, pubkey, prvkey, score));		
	}

	paillier_freeciphertext(tmp_score);

	return score;
}


// 150809. added by KHI.
paillier_ciphertext_t** protocol::Smax_n(paillier_ciphertext_t*** ciper, int number){
	int i, j;
	int iter;
	int num = number;
	int k;
	paillier_ciphertext_t*** copy_ciper = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*num);
	for(i=0; i<num ; i++){
		copy_ciper[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		for(j=0; j<size;j++){
			copy_ciper[i][j] = ciper[i][j];
		}
	}
	int e, q;
	
	paillier_ciphertext_t** sbd = SBD(paillier_create_enc_zero());

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

				copy_ciper[e] = Smax_basic1(copy_ciper[e],copy_ciper[q]);
				copy_ciper[q] = sbd;
				//printf("%d %d\n",e,q);
			} else {
				e = (int)pow(2, i) * (j - 1);
				q = (int)pow(2, i) * j - (int)pow(2, i - 1);

				copy_ciper[e] = Smax_basic1(copy_ciper[e],copy_ciper[q]);
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

// 150809. added by KHI.
paillier_ciphertext_t** protocol::Smax_basic1(paillier_ciphertext_t** cipher1, paillier_ciphertext_t** cipher2)
{
	int i;
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
			paillier_subtract(pubkey, ciper_W[i], cipher1[i], SM_p1(cipher1[i], cipher2[i]));	// W
			paillier_subtract(pubkey,temp[i],cipher2[i],cipher1[i]);	
			paillier_mul(pubkey,ciper_R[i],temp[i],ciper_Rand_value);	// Gamma
		}else{
			paillier_subtract(pubkey, ciper_W[i], cipher2[i], SM_p1(cipher1[i], cipher2[i]));
			paillier_subtract(pubkey, temp[i], cipher1[i], cipher2[i]);
			paillier_mul(pubkey, ciper_R[i], temp[i], ciper_Rand_value);
		}
		ciper_G[i]=SBXOR(cipher1[i],cipher2[i]);

		if(i==0){
			paillier_exp(pubkey,ciper_H[i], cipher_zero, Rand_value);
			paillier_mul(pubkey,ciper_H[i], ciper_H[i], ciper_G[i]);		
		}else{
			paillier_exp(pubkey,temp[i],ciper_H[i-1],Rand_value);
			paillier_mul(pubkey,ciper_H[i],temp[i],ciper_G[i]);
		}

		paillier_mul(pubkey,ciper_O[i],ciper_H[i],cipher_minus);	// PI
		paillier_exp(pubkey,ciper_O[i],ciper_O[i],Rand_value);
		// paillier_exp(pubkey,temp3[i],ciper_W[i],Rand_value);   // our IDEA
		// paillier_mul(pubkey,ciper_L[i],temp2[i], temp3[i]);	 // our IDEA
		
		paillier_mul(pubkey, ciper_L[i], ciper_O[i], ciper_W[i]);
	}

	alpha = Smax_basic2(ciper_R, ciper_L, ciper_M, alpha);

	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));
	
	for(i=0;i<size;i++){
		paillier_exp(pubkey,tmp,alpha,Rand_value);
		paillier_subtract(pubkey, ciper_lambda[i], ciper_M[i], tmp);
		if(func){
			paillier_mul(pubkey,ciper_min[i],cipher1[i],ciper_lambda[i]);			
		}else{
			paillier_mul(pubkey,ciper_min[i],cipher2[i],ciper_lambda[i]);			
		}
	}
/*
	for(i=0; i<size; i++){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_min[i]));
	}
	printf("\n");
*/
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

// 150809. added by KHI.
paillier_ciphertext_t* protocol::Smax_basic2(paillier_ciphertext_t** ciper_R, paillier_ciphertext_t** ciper_L, paillier_ciphertext_t** ciper_M, paillier_ciphertext_t* alpha)
{
	int i;
	int flag = 0;

	// compute alpha
	for(i=0; i<size; i++)	{
		if (strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, ciper_L[i])), paillier_plaintext_to_str(plain_one)) == 0) {
			alpha = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
			//paillier_print("alpha : ", alpha);
			
			flag = 1;
			
			break;
		}
		else {
			alpha = paillier_enc(0, pubkey, plain_one, paillier_get_rand_devurandom);
			//paillier_print("alpha : ", alpha);
		}
	}

	if(flag == 1) {
		for(i=0; i<size; i++){
			paillier_exp(pubkey, ciper_M[i], ciper_R[i], plain_zero);
		}
	}
	else {
		for(i=0; i<size; i++){
			paillier_exp(pubkey, ciper_M[i], ciper_R[i], plain_one);
		}
	}

	//paillier_print("alpha : ", alpha);

	return alpha;
}
int** protocol::STopk_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode) {
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;

	paillier_ciphertext_t * hint = paillier_create_enc_zero();
	paillier_ciphertext_t*** cand;
	paillier_ciphertext_t*** temp_cand ;

	paillier_ciphertext_t** coeff = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t*** ciper_coeff_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t** temp_q = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** ciper_hint_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** ciper_kthscore_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** ciper_tempscore_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t*** ciper_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** ciper_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t**** ciper_nodeLL_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t**** ciper_nodeRR_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);
	paillier_ciphertext_t* cipher_binary;
	paillier_ciphertext_t* cipher_MAX	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_scores	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_score = paillier_create_enc_zero();
	paillier_ciphertext_t* kth_score = 0;

	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	paillier_ciphertext_t* temp_alpha;
	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);	 // 프사이
	for(i=0; i<NumNode; i++){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);
		ciper_nodeLL_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
		ciper_nodeRR_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	}
	for( i = 0 ; i < k ; i++ ){
		ciper_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim; j ++ ){
			ciper_result[i][j] = paillier_create_enc_zero();
		}
	}
	for(i=0; i<dim; i++){
		coeff[i] = paillier_create_enc_zero();
		temp_q[i] = paillier_create_enc_zero();
		psi[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);
	}
	


	startTime = std::chrono::system_clock::now(); // startTime check
	hint = q[dim];
	if(Print){
		gmp_printf("hint : %Zd  /  %Zd\n", paillier_dec(0, pubkey, prvkey, hint),paillier_dec(0, pubkey, prvkey, q[dim]));
		cout <<"q[i] + hint : "<<endl;
	}
	for(i=0; i<dim; i++) {
		paillier_mul(pubkey, coeff[i], q[i], hint);
		if(Print){
			gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, coeff[i]));
			gmp_printf("max val : %Zd \n", paillier_dec(0, pub, prv, max_val[i]));
		}
	}
	cout <<endl;
	ciper_hint_bit = SBD_for_SRO(hint, 0);
	if(Print){
		printf("hint (SBD) : ");
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_hint_bit[j]));
		}
		printf("\n");
	}
	for(i=0; i<dim; i++)
	{
		ciper_coeff_bit[i] = SBD_for_SRO(coeff[i], 0);			// query(coeff) 비트 변환 수행
		psi[i] = SCMP(ciper_hint_bit, ciper_coeff_bit[i]);		// 프사이 계산 (1이면 계수가 양수인 항, 0이면 계수가 음수인 항)
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["sTopk_I_preprocessing_query"] = time_variable.find("sTopk_I_preprocessing_query")->second + duration_sec.count();


	startTime = std::chrono::system_clock::now(); // startTime check
	// 각 노드 비트 변환 수행
	cout << "\n\nnode.LL, node.RR SBD Start" <<endl;
	for(i=0; i<NumNode; i++) {	
		for(j=0; j<dim; j++) {
			ciper_nodeLL_bit[i][j] = SBD_for_SRO(node[i].LL[j], 0);				// node LL bound 변환
			ciper_nodeRR_bit[i][j] = SBD_for_SRO(node[i].RR[j], 1);	 			// node RR bound 변환
		}
		if(Print)
		{
			printf("%dth node LL bound\n", i);
			for(j=0; j<dim; j++) {
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, node[i].LL[j]));
			}
			cout<<endl;
			for(j=0; j<dim; j++) {
				for(m=0; m<size+1; m++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeLL_bit[i][j][m]));
				}
				cout<<endl;
			}
			printf("%dth node RR bound\n", i);
			for(j=0; j<dim; j++) {
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, node[i].RR[j]));
			}
			cout<<endl;
			for(j=0; j<dim; j++) {
				for(m=0; m<size+1; m++) {
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_nodeRR_bit[i][j][m]));
				}
				printf("\n");
			}
		}
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["sTopk_I_SBD_node"] = time_variable.find("sTopk_I_SBD_node")->second + duration_sec.count();


	startTime = std::chrono::system_clock::now(); // startTime check	
	for(i=0; i<dim; i++) {	
		temp_q[i] = SM_p1(max_val[i], psi[i]);	
	}
	for(i=0; i<dim; i++) {
		ciper_qLL_bit[i] = SBD_for_SRO(temp_q[i], 0);		// query LL bound 변환
		ciper_qRR_bit[i] = SBD_for_SRO(temp_q[i], 1);		// query RR bound 변환
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["sTopk_I_SBD_query"] = time_variable.find("sTopk_I_SBD_query")->second + duration_sec.count();


	startTime = std::chrono::system_clock::now(); // startTime check	

	for(i=0; i<NumNode; i++) {	
		alpha[i] = SRO(ciper_qLL_bit, ciper_qRR_bit, ciper_nodeLL_bit[i], ciper_nodeRR_bit[i]);		
		if(!verify_flag)	{	//  검증 단계에서는 수행하지 않음
			for(j=0; j<dim; j++) {
				node[i].LL[j] = SM_p1(node[i].LL[j], SBN(alpha[i]));
				node[i].RR[j] = SM_p1(node[i].RR[j], SBN(alpha[i]));
			}
		}
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["sTopk_I_SRO_node_query"] = time_variable.find("sTopk_I_SRO_node_query")->second + duration_sec.count();


	while(1) {
		cnt = 0;
		if(!verify_flag)	{	// 검증 단계가 아닐 시에는, cand에 저장
			startTime = std::chrono::system_clock::now(); // startTime check	
			cand	= sNodeRetrievalforTopk(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_I_extract_cand"] = time_variable.find("sTopk_I_extract_cand")->second + duration_sec.count();
	
		}else {		// 검증 단계일 시에는, 이전 결과와 cand를 합침
			startTime = std::chrono::system_clock::now(); // startTime check	
			temp_cand = sNodeRetrievalforTopk(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_I_extract_expand_cand"] = time_variable.find("sTopk_I_extract_expand_cand")->second + duration_sec.count();
	
			printf("cnt : %d \n", cnt);
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
		paillier_ciphertext_t*** ciper_SBD_score	= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);
		paillier_ciphertext_t*** ciper_V2			= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);
		paillier_ciphertext_t** ciper_score			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t** ciper_V				=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t** ciper_mid			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t** ciper_Smax			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		for( i = 0 ; i < cnt ; i++ ){
			ciper_score[i] 		= 	(paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			ciper_SBD_score[i]	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
			ciper_V[i]			= cipher_zero;
			ciper_V2[i]			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		}
		for( i = 0 ; i < size; i++ ){
			ciper_Smax[i] = cipher_zero;
		}

		startTime = std::chrono::system_clock::now(); // startTime check	
		for( i = 0 ; i < cnt ; i++ ){
			ciper_score[i] = computeScore(q, cand[i]);
			ciper_SBD_score[i] = SBD(ciper_score[i]);
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["sTopk_I_computescore&SBD"] = time_variable.find("sTopk_I_computescore")->second + duration_sec.count();
	
		for( s = 0 ; s < k ; s++ ){
			printf("\n%dth sTopk start \n", s+1);
			startTime = std::chrono::system_clock::now(); // startTime check	
			ciper_Smax = Smax_n(ciper_SBD_score, cnt);	// bit로 표현된 암호화 max score 추출
			if(!verify_flag){
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["sTopk_I_smax_n"] = time_variable.find("sTopk_I_smax_n")->second + duration_sec.count();
			}else{
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["sTopk_I_smax_n_expand"] = time_variable.find("sTopk_I_smax_n_expand")->second + duration_sec.count();
			}


			startTime = std::chrono::system_clock::now(); // startTime check	
			// bit로 표현된 암호화 max score를 암호화 정수로 변환
			for( j = size ; j > 0 ; j-- ){
				t = (int)pow(2, j-1);
				cipher_binary = paillier_create_enc(t);
				cipher_binary = SM_p1(cipher_binary, ciper_Smax[size-j]);
				paillier_mul(pubkey, cipher_MAX, cipher_binary, cipher_MAX);
				paillier_freeciphertext(cipher_binary);
				//gmp_printf("min : %Zd\n", cipher_MAX);
			}
			gmp_printf("MAX Value score: %Zd\n", paillier_dec(0, pubkey, prvkey, cipher_MAX));

			// 이전 iteration에서 max값으로 선택된 데이터의 score가 secure하게 MIN으로 변환되었기 때문에, 
			// 데이터의 점수 계산을 모두 다시 수행
			if( s != 0 ){
				for( i = 0 ; i < cnt; i++ ){
					temp_score = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
					for( j = size ; j > 0 ; j--){
						t = (int)pow(2, j-1);
						cipher_binary = paillier_create_enc(t);
						cipher_binary = SM_p1(cipher_binary, ciper_SBD_score[i][size-j]);
						paillier_mul(pubkey, temp_score, cipher_binary, temp_score);
					}
					ciper_score[i] = temp_score;
					temp_score = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_I_recovery_score"] = time_variable.find("sTopk_I_recovery_score")->second + duration_sec.count();

			startTime = std::chrono::system_clock::now(); // startTime check	
			// 각 데이터의 점수와 max score와의 차를 구함 (max 데이터의 경우에만 0으로 만들기 위함)
			for( i = 0 ; i < cnt ; i++ ){
				paillier_subtract(pubkey, temp_score, ciper_score[i], cipher_MAX);
				ciper_mid[i] = SM_p1(temp_score, ciper_rand);
				if(Print)
				{
					gmp_printf("dist-dist : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_mid[i]));
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_I_subtract"] = time_variable.find("sTopk_I_subtract")->second + duration_sec.count();

			ciper_V = Topk_sub(ciper_mid, cnt);


			startTime = std::chrono::system_clock::now(); // startTime check	
			// max 데이터 추출
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0 ; j < dim; j++ ){
					ciper_V2[i][j] = SM_p1(ciper_V[i], cand[i][j]);
					if(i==0) {
						ciper_result[s][j] = ciper_V2[i][j];
					}
					else {
						paillier_mul(pubkey, ciper_result[s][j], ciper_V2[i][j], ciper_result[s][j]);
					}
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_I_extract_topk"] = time_variable.find("sTopk_I_extract_topk")->second + duration_sec.count();


			startTime = std::chrono::system_clock::now(); // startTime check	
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0; j < size ; j++ ){
					ciper_SBD_score[i][j] = SM_p1(SBN(ciper_V[i]), ciper_SBD_score[i][j]);
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_I_updateScore"] = time_variable.find("sTopk_I_updateScore")->second + duration_sec.count();
			cipher_MAX	= paillier_create_enc_zero();
		}
		if(!verify_flag){
			startTime = std::chrono::system_clock::now(); // startTime check	
			if(cnt < k)	 // 노드의 fanout이 k보다 적은 경우를 처리함
				kth_score = computeScore(q, ciper_result[cnt-1]);
			else	// k개가 찾아졌다면, k번째 데이터의 score를 계산
				kth_score = computeScore(q, ciper_result[k-1]);
			gmp_printf("%dth score : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_score));
			for(i=0; i<NumNode; i++) {	
				for(j=0; j<dim; j++) {
					paillier_mul(pubkey, temp_q[j], SM_p1(node[i].RR[j], psi[j]), SM_p1(node[i].LL[j], SBN(psi[j])));	
				}
				ciper_kthscore_bit = SBD_for_SRO(kth_score, 0);		// k번째 score SBD 수행
				ciper_tempscore_bit = SBD_for_SRO(computeScore(q, temp_q), 1);  // 각 노드의 score를 SBD 수행
				alpha[i] = SCMP(ciper_kthscore_bit , ciper_tempscore_bit);	// k번째 score보다 높을 가능성이 있으면 alpha=1
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_I_check_expand_cond"] = time_variable.find("sTopk_I_check_expand_cond")->second + duration_sec.count();
		}
		if(!verify_flag)	{
			verify_flag = true;		// 검증을 수행하기 위해 flag를 true로 변경
		}
		else
			break;
	}
	for(i=0; i<k; i++){
		printf("%d final result : ", i);
		gmp_printf(" %Zd : ", paillier_dec(0, pubkey, prvkey, computeScore(q, ciper_result[i])));
		for(j=0; j<dim; j++){
			gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, ciper_result[i][j]));
			paillier_mul(pubkey,ciper_result[i][j],ciper_result[i][j],ciper_rand);
		}
		printf("\n");
	}
	
	free(coeff);
	free(ciper_coeff_bit);
	free(temp_q);
	free(ciper_hint_bit);
	free(ciper_kthscore_bit);
	free(ciper_tempscore_bit);
	free(ciper_qLL_bit);
	free(ciper_qRR_bit);
	free(ciper_nodeLL_bit);
	free(ciper_nodeRR_bit);
	free(pt);
	free(cipher_MAX);
	free(ciper_scores);
	free(ciper_rand);
	free(temp_score);
	free(alpha);
	free(psi);	
	cout << "End Line" <<endl;
	return SkNNm_Bob2(ciper_result, rand, k, dim);
}
// Proposed skNN with secure Index
int** protocol::STopk(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode) {
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;

	bool verify_flag = false;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_ciphertext_t * hint = paillier_create_enc_zero();
	paillier_ciphertext_t*** cand;
	paillier_ciphertext_t*** temp_cand ;

	paillier_ciphertext_t** coeff = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t*** ciper_coeff_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t** temp_q = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** ciper_hint_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** ciper_kthscore_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** ciper_tempscore_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t*** ciper_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** ciper_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t**** ciper_nodeLL_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t**** ciper_nodeRR_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	
	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* cipher_binary;
	paillier_ciphertext_t* cipher_MAX	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_scores	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_score = paillier_create_enc_zero();
	paillier_ciphertext_t* kth_score = 0;

	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	paillier_ciphertext_t* temp_alpha;
	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);	 // 프사이
	for(i=0; i<NumNode; i++){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);

		ciper_nodeLL_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
		ciper_nodeRR_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	}
			
	for( i = 0 ; i < k ; i++ ){
		ciper_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim; j ++ ){
			ciper_result[i][j] = paillier_create_enc_zero();
		}
	}

	for(i=0; i<dim; i++){
		coeff[i] = paillier_create_enc_zero();
		temp_q[i] = paillier_create_enc_zero();

		psi[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);
	}
	
	hint = q[dim]; 

	// query의 각 차원 데이터와 hint를 더함
	for(i=0; i<dim; i++) {
		paillier_mul(pubkey, coeff[i], q[i], hint);
	}

	// debugging
	/*
	for(i=0; i<dim; i++) {
		gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, coeff[i]));
		printf("\n");
	}
	*/

	ciper_hint_bit = SBD_for_SRO(hint, 0);		// hint SBD 수행
	
	// debugging
	printf("hint (SBD) : ");
	for(j=0; j<size+1; j++) {
		gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_hint_bit[j]));
	}
	printf("\n");
	
	for(i=0; i<dim; i++) {
		ciper_coeff_bit [i] = SBD_for_SRO(coeff[i], 0);			// query(coeff) 비트 변환 수행
		
		// debugging
		printf("coeff (SBD) : ");
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_coeff_bit [i][j]));
		}
		printf("\n");

		psi[i] = SCMP(ciper_hint_bit, ciper_coeff_bit[i]);		// 프사이 계산 (1이면 계수가 양수인 항, 0이면 계수가 음수인 항)
	}
	
	for(i=0; i<dim; i++) {
		gmp_printf("psi : %Zd\n", paillier_dec(0, pubkey, prvkey, psi[i]));
	}

	startTime = clock();

	// 각 노드 비트 변환 수행 
	for(i=0; i<NumNode; i++) {	
		for(j=0; j<dim; j++) {
			ciper_nodeLL_bit[i][j] = SBD_for_SRO(node[i].LL[j], 0);				// node LL bound 변환
			ciper_nodeRR_bit[i][j] = SBD_for_SRO(node[i].RR[j], 1);	 			// node RR bound 변환
		}

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
	}

	endTime = clock();
	gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("node SBD time: %f\n", gap);

	// 첫 탐색 노드를 찾기 위한 질의 셋팅(양수 계수는 Max, 음수 계수는 0)
	for(i=0; i<dim; i++) {	
		temp_q[i] = SM_p1(max_val[i], psi[i]);	
	}

	// debugging
	printf("===query for node searching===\n");
	for(i=0; i<dim; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, temp_q[i]));
	}
	printf("\n");

	// query 비트 변환 수행
	for(i=0; i<dim; i++) {
		ciper_qLL_bit[i] = SBD_for_SRO(temp_q[i], 0);		// query LL bound 변환
		ciper_qRR_bit[i] = SBD_for_SRO(temp_q[i], 1);		// query RR bound 변환
	}
	


	// score가 가장 높을 것으로 예상되는 노드 탐색
	// 각 노드 SRO 호출
	startTime = clock();
	for(i=0; i<NumNode; i++) {	
		alpha[i] = SRO(ciper_qLL_bit, ciper_qRR_bit, ciper_nodeLL_bit[i], ciper_nodeRR_bit[i]);		

		if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
			printf("%dth node overlaps the query region.\n", i);

		/*
		if(!verify_flag)	{	//  검증 단계에서는 수행하지 않음
			// 이미 검색이 완료된 노드가 재 탐색되는 것을 방지하기 위해, bound를 MAX로 변환
			for(j=0; j<dim; j++) {
				for(m=0; m<size+1; m++) {
					ciper_nodeLL_bit[i][j][m] = SBOR(ciper_nodeLL_bit[i][j][m], alpha[i]);
					ciper_nodeRR_bit[i][j][m] = SBOR(ciper_nodeRR_bit[i][j][m], alpha[i]);
				}
			}
		}
		*/

		if(!verify_flag)	{	//  검증 단계에서는 수행하지 않음
			// 이미 검색이 완료된 노드가 재 탐색되는 것을 방지하기 위해, bound를 0으로 변환
			for(j=0; j<dim; j++) {
				node[i].LL[j] = SM_p1(node[i].LL[j], SBN(alpha[i]));
				node[i].RR[j] = SM_p1(node[i].RR[j], SBN(alpha[i]));
			}
		}

		// debugging
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
	}
	endTime = clock();
	gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("SRO time : %f\n", gap);

	printf("test  (%d line)\n", __LINE__);


	printf("alpha\n");
	for(i=0; i<NumNode; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, alpha[i]));
	}
	printf("\n");

	while(1) {
		cnt = 0;

		if(!verify_flag)	{	// 검증 단계가 아닐 시에는, cand에 저장
			cand	= sNodeRetrievalforTopk(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
		}
		else {		// 검증 단계일 시에는, 이전 결과와 cand를 합침
			temp_cand	= sNodeRetrievalforTopk(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);

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

		// debugging
		for(i=0; i<cnt; i++) {
			for( j = 0 ; j < dim; j ++ ){
				gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
			}
			printf("\n");
		}

		paillier_ciphertext_t** ciper_score = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t*** ciper_SBD_score = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);
		paillier_ciphertext_t** ciper_V=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t*** ciper_V2 = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);

		paillier_ciphertext_t** ciper_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t** ciper_Smax = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		
		for( i = 0 ; i < size; i++ ){
			ciper_Smax[i] = cipher_zero;
		}

		for( i = 0 ; i < cnt ; i++ ){
			ciper_score[i] 	= 	(paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(alpha[i]->c);

			ciper_SBD_score[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
			ciper_V[i]	= cipher_zero;
			ciper_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		}

		printf("\nSSED & SBD start \n");
		startTime = clock();
		for( i = 0 ; i < cnt ; i++ ){
			//ciper_score[i] = computeScore2(q, cand[i], coeff, hint);
			ciper_score[i] = computeScore(q, cand[i]);
			ciper_SBD_score[i] = SBD(ciper_score[i]);
	
			gmp_printf("score : %Zd\n",  paillier_dec(0, pubkey, prvkey, ciper_score[i]));	
			for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_SBD_score[i][j]));
			}
			printf("\n");
		}
		endTime = clock();
		gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		printf("SSED & SBD time : %f\n", gap);


		for( s = 0 ; s < k ; s++ ){
			printf("\n%dth sTopk start \n", s+1);
			startTime = clock();
			ciper_Smax = Smax_n(ciper_SBD_score, cnt);	// bit로 표현된 암호화 max score 추출
			endTime = clock();
			gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("%dth sMAXn time : %f\n", s+1, gap);
		
			
			for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smax[j]));
			}
			printf("\n");	

			// bit로 표현된 암호화 max score를 암호화 정수로 변환
			for( j = size ; j > 0 ; j-- ){
				t = (int)pow(2, j-1);
				cipher_binary = paillier_create_enc(t);
				cipher_binary = SM_p1(cipher_binary, ciper_Smax[size-j]);
				paillier_mul(pubkey, cipher_MAX, cipher_binary, cipher_MAX);
				paillier_freeciphertext(cipher_binary);
				//gmp_printf("min : %Zd\n", cipher_MAX);
			}
			gmp_printf("max score: %Zd\n", paillier_dec(0, pubkey, prvkey, cipher_MAX));


			// 이전 iteration에서 max값으로 선택된 데이터의 score가 secure하게 MIN으로 변환되었기 때문에, 
			// 데이터의 점수 계산을 모두 다시 수행
			if( s != 0 ){
				for( i = 0 ; i < cnt; i++ ){
					for( j = size ; j > 0 ; j--){
						t = (int)pow(2, j-1);
						cipher_binary = paillier_create_enc(t);
						cipher_binary = SM_p1(cipher_binary, ciper_SBD_score[i][size-j]);
						paillier_mul(pubkey, temp_score, cipher_binary, temp_score);
					}
					ciper_score[i] = temp_score;
					temp_score = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}

			for( i = 0 ; i < cnt ; i++ ){
				gmp_printf("score : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_score[i]) );
			}

			// 각 데이터의 점수와 max score와의 차를 구함 (max 데이터의 경우에만 0으로 만들기 위함)
			for( i = 0 ; i < cnt ; i++ ){
				paillier_subtract(pubkey, temp_score, ciper_score[i], cipher_MAX);
				ciper_mid[i] = SM_p1(temp_score, ciper_rand);
				//gmp_printf("dist - dist : ",ciper_mid[i]);
			}

			ciper_V = SkNNm_sub(ciper_mid, cnt);
			
			
			for( i = 0 ; i < cnt ; i++ ){
				gmp_printf("%dth ciper_V : %Zd\n", i, paillier_dec(0, pubkey, prvkey, ciper_V[i]));
			}

			// max 데이터 추출
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0 ; j < dim; j++ ){
					gmp_printf("cand : %Zd (%d line)\n", paillier_dec(0, pubkey, prvkey, cand[i][j]), __LINE__);
					ciper_V2[i][j] = SM_p1(ciper_V[i], cand[i][j]);
					if(i==0) {
						ciper_result[s][j] = ciper_V2[i][j];
					}
					else {
						paillier_mul(pubkey, ciper_result[s][j], ciper_V2[i][j], ciper_result[s][j]);
					}
				}
			}

			printf("test\n");
		
			for( j = 0 ; j < dim; j++ ){
				gmp_printf("ciper_result : %Zd ", paillier_dec(0, pubkey, prvkey, ciper_result[s][j]));
			}
			printf("\n");
			
			
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0; j < size ; j++ ){
					ciper_SBD_score[i][j] = SM_p1(SBN(ciper_V[i]), ciper_SBD_score[i][j]);
				}
			}

			cipher_MAX = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
		}


		if(!verify_flag)	{	
			if(cnt < k)	 // 노드의 fanout이 k보다 적은 경우를 처리함
				kth_score = computeScore(q, ciper_result[cnt-1]);
			else	// k개가 찾아졌다면, k번째 데이터의 score를 계산
				kth_score = computeScore(q, ciper_result[k-1]);

			gmp_printf("%dth score : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_score));

			
			for(i=0; i<NumNode; i++) {	
				for(j=0; j<dim; j++) {
					printf("test\n");
					// 확장 탐색 노드를 찾기 위한 질의 셋팅(양수 계수는 해당 노드의 UB, 음수 계수는 해당 노드의 LB)
					gmp_printf("%d line: %Zd\n", __LINE__, paillier_dec(0, pubkey, prvkey, SM_p1(node[i].RR[j], psi[j])));
					gmp_printf("%d line: %Zd\n", __LINE__, paillier_dec(0, pubkey, prvkey, SM_p1(node[i].LL[j], SBN(psi[j]))));

					paillier_mul(pubkey, temp_q[j], SM_p1(node[i].RR[j], psi[j]), SM_p1(node[i].LL[j], SBN(psi[j])));	
					gmp_printf("%d line: %Zd\n	", __LINE__, paillier_dec(0, pubkey, prvkey, temp_q[j]));
				}

				ciper_kthscore_bit = SBD_for_SRO(kth_score, 0);		// k번째 score SBD 수행
				//ciper_kthscore_bit = SBD(kth_score);		// k번째 score SBD 수행

				gmp_printf("computeScore(q, temp_q) : %Zd %d line\n", paillier_dec(0, pubkey, prvkey, computeScore(q, temp_q)), __LINE__);
				ciper_tempscore_bit = SBD_for_SRO(computeScore(q, temp_q), 1);  // 각 노드의 score를 SBD 수??
				//ciper_tempscore_bit = SBD(computeScore(q, temp_q));  // 각 노드의 score를 SBD 수행


				alpha[i] = SCMP(ciper_kthscore_bit , ciper_tempscore_bit);	// k번째 score보다 높을 가능성이 있으면 alpha=1

				gmp_printf("alpha : %Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
			}
		}

		if(!verify_flag)	{
			verify_flag = true;		// 검증을 수행하기 위해 flag를 true로 변경
		}
		else
			break;	// 검증을 끝냄
	}
		
	// user(Bob)에게 결과 전송을 위해 random 값 삽입
	for(i=0; i<k; i++){
		printf("%d final result : ", i);
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,ciper_result[i][j],ciper_result[i][j],ciper_rand);
			gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, ciper_result[i][j]));
		}
		printf("\n");
	}
	paillier_freeciphertext(temp_score);	
	
	return SkNNm_Bob2(ciper_result, rand, k, dim);
}
	


// Proposed skNN with secure Index + SMSn
int** protocol::STopk_G(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode) {
	printf("\n=== STopk_G start ===\n");
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;

	//쿼리 및 hint 세팅
	paillier_ciphertext_t*	MIN				= paillier_create_enc(0);
	paillier_ciphertext_t*	TMP_alpha		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(TMP_alpha->c);
	paillier_ciphertext_t*	C_RAND			= paillier_create_enc(rand);
	paillier_ciphertext_t*	hint			= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(hint->c);
	paillier_ciphertext_t* kth_score		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(kth_score->c);	
	
	paillier_ciphertext_t** Q				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** coeff			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t** psi				= (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);	 // 프사이
	paillier_ciphertext_t** alpha			= (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	paillier_ciphertext_t***RESULT			= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t***cand;
	paillier_ciphertext_t***temp_cand ;

	for( i = 0 ; i < dim ; i++ )
	{
		Q[i]		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		psi[i]		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		coeff[i]	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(Q[i]->c);
		mpz_init(psi[i]->c);
		mpz_init(coeff[i]->c);
	}
	for( i = 0 ; i < k ; i++ )
	{
		RESULT[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j++ )
		{
			RESULT[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(RESULT[i][j]->c);
		}
	}
	for( i = 0 ; i < NumNode ; i++ )
	{
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);
	}
	
	startTime = std::chrono::system_clock::now(); // startTime check
	hint = q[dim];
	if(Print == false)
	{
		gmp_printf("hint : %Zd  /  %Zd\n", paillier_dec(0, pubkey, prvkey, hint),paillier_dec(0, pubkey, prvkey, q[dim]));
	}
	for( i = 0 ; i < dim ; i++ )
	{
		paillier_mul(pubkey, coeff[i], q[i], hint);
		if(Print == false)
		{
			cout <<"q[i] + hint : "<<endl;
			gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, coeff[i]));
			gmp_printf("max val : %Zd \n", paillier_dec(0, pub, prv, max_val[i]));
		}
	}
	for( i = 0 ; i < dim ; i++ )
	{
		psi[i] = GSCMP(hint, coeff[i]);
		if(Print == false)
		{
			gmp_printf("psi : %Zd", paillier_dec(0, pubkey, prvkey, psi[i]));
			if ( strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, psi[i]) ), paillier_plaintext_to_str(plain_one)) == 0 )
			{
				cout << " Q["<< i <<"] > 0 " <<endl; 
			}else{
				cout << " Q["<< i <<"] < 0 " <<endl; 
			}
		}
	}

	for( i = 0 ; i < dim ; i++ )
	{	
		Q[i] = SM_p1(max_val[i], psi[i]);
	}

	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["sTopk_G_preprocessing_query"] = time_variable.find("sTopk_G_preprocessing_query")->second + duration_sec.count();



	startTime = std::chrono::system_clock::now(); // startTime check
	for( i = 0 ; i < NumNode ; i++ )
	{	
		alpha[i] = GSRO(Q, Q, node[i].LL, node[i].RR);
		//gmp_printf("alpha: %Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));

		if(!verify_flag)
		{	//  검증 단계에서는 수행하지 않음
			for( j = 0 ; j < dim ; j++ )
			{
				node[i].LL[j] = SM_p1(node[i].LL[j], SBN(alpha[i]));
				node[i].RR[j] = SM_p1(node[i].RR[j], SBN(alpha[i]));
			}
		}
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["sTopk_G_GSRO_node_query"] = time_variable.find("sTopk_G_GSRO_node_query")->second + duration_sec.count();




	while(1)
	{
		cnt = 0;
		if(!verify_flag)
		{	
			startTime = std::chrono::system_clock::now(); // startTime check
			cand = GSRO_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_G_extract_cand"] = time_variable.find("sTopk_G_extract_cand")->second + duration_sec.count();
		}
		else
		{
			startTime = std::chrono::system_clock::now(); // startTime check
			temp_cand = GSRO_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_G_extract_expand_cand"] = time_variable.find("sTopk_G_extract_expand_cand")->second + duration_sec.count();

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
					cand[i][j] = RESULT[i][j];
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

		int MAX_idx = 0;
		paillier_ciphertext_t*	MAX				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(MAX->c);

		paillier_ciphertext_t**	V				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t** SCORE			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t**	SCORE_MINUS_MAX	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t***V2				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);

		for( i = 0 ; i < cnt ; i++ )
		{
			V[i]				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			SCORE[i]			= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			SCORE_MINUS_MAX[i]	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			V2[i]				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);

			mpz_init(V[i]->c);
			mpz_init(SCORE[i]->c);
			mpz_init(SCORE_MINUS_MAX[i]->c);
			for( j = 0 ; j < dim ; j++ )
			{
				V2[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
				mpz_init(V2[i][j]->c);
			}
		}

		startTime = std::chrono::system_clock::now(); // startTime check	
		for( i = 0 ; i < cnt ; i++ )
		{
			SCORE[i] = computeScore(q, cand[i]);
			if(Print)
			{
				gmp_printf("SCORE : %Zd\n",  paillier_dec(0, pubkey, prvkey, SCORE[i]));	
			}
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["sTopk_G_computeScore"] = time_variable.find("sTopk_G_computeScore")->second + duration_sec.count();

		for( s = 0 ; s < k ; s++ )
		{
			cout << s+1<<"dth MAXn_Topk start"<<endl;
			startTime = std::chrono::system_clock::now(); // startTime check	
			MAX_idx = MAXn(SCORE, cnt);
			MAX = SCORE[MAX_idx];

			gmp_printf("MAX : %Zd \n", paillier_dec(0, pubkey, prvkey, MAX));
			if(!verify_flag) {
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["sTopk_G_sMaxn"] = time_variable.find("sTopk_G_sMaxn")->second + duration_sec.count();
			}
			else {
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["sTopk_G_sMaxn_second"] = time_variable.find("sTopk_G_sMaxn_second")->second + duration_sec.count();
			}
			

			startTime = std::chrono::system_clock::now(); // startTime check	
			for( i = 0 ; i < cnt ; i++ )
			{
				paillier_subtract(pubkey, SCORE_MINUS_MAX[i], SCORE[i], MAX);
				SCORE_MINUS_MAX[i] = SM_p1(SCORE_MINUS_MAX[i], C_RAND);
				if(Print)
				{
					gmp_printf("SCORE-MAX : %Zd \n", paillier_dec(0, pubkey, prvkey, SCORE_MINUS_MAX[i]));
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_G_subtract"] = time_variable.find("sTopk_G_subtract")->second + duration_sec.count();


			
			V = Topk_sub(SCORE_MINUS_MAX, cnt);


			startTime = std::chrono::system_clock::now(); // startTime check	
			for( i = 0 ; i < cnt ; i++ )
			{
				paillier_subtract(pubkey, TMP_alpha, cipher_one, V[i]);
				paillier_mul(pubkey, SCORE[i], SM_p1(V[i], MIN), SM_p1(TMP_alpha, SCORE[i]));
				for( j = 0 ; j < dim; j++ )
				{
					V2[i][j] = SM_p1(V[i], cand[i][j]);
					if( i == 0 )
					{
						RESULT[s][j] = V2[i][j];
					}
					else
					{
						paillier_mul(pubkey, RESULT[s][j], V2[i][j], RESULT[s][j]);
					}
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_G_extract_topk"] = time_variable.find("sTopk_G_extract_topk")->second + duration_sec.count();
		}

		if(!verify_flag)
		{	
			startTime = std::chrono::system_clock::now(); // startTime check	

			if(cnt < k)	 // 노드의 fanout이 k보다 적은 경우를 처리함
				kth_score = computeScore(q, RESULT[cnt-1]);
			else	// k개가 찾아졌다면, k번째 데이터의 score를 계산
				kth_score = computeScore(q, RESULT[k-1]);
			//gmp_printf("%dth score : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_score));
			for( i = 0 ; i < NumNode ; i++ )
			{	
				for( j = 0 ; j < dim ; j++ )
				{
					paillier_mul(pubkey, Q[j], SM_p1(node[i].RR[j], psi[j]), SM_p1(node[i].LL[j], SBN(psi[j])));	
				}
				alpha[i] = GSCMP(kth_score , computeScore(q, Q));	// k번째 score보다 높을 가능성이 있으면 alpha=1
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["sTopk_G_check_expand_cond"] = time_variable.find("sTopk_G_check_expand_cond")->second + duration_sec.count();
		}
		if(!verify_flag)
		{
			verify_flag = true;
		}else
			break;
	}
	for( i = 0 ; i < k ; i++ )
	{
		printf("%d final result : ", i);
		gmp_printf(" %Zd : \n", paillier_dec(0, pubkey, prvkey, computeScore(q, RESULT[i])));
	
		for( j = 0 ; j < dim ; j++ )
		{
			paillier_mul(pubkey, RESULT[i][j], RESULT[i][j], C_RAND);
		}	
	}
	cout << "End Line" <<endl;
	return SkNNm_Bob2(RESULT, rand, k, dim);
}

int protocol::MAXn(paillier_ciphertext_t** ciper, int cnt){
	int i = 0, s = 0;
	int buff = 0;
	int R1 = 5;
	int quot = 0;
	int remain = 0;
	int iter = 0;
	int bundle = 0;
	int ciper_idx = 0;
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
/*
	if(Print)
	{
		cout << " quot : "<< quot << " remain : " << remain << "  " << buff <<endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!MAXn Start!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		for( i = 0 ; i < cnt ; i++ ){
			cout << i+1 << " ciper : ";
			gmp_printf(" %Zd\n",paillier_dec(0, pub, prv, ciper[i]));
		}
	}
*/
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
/*
	if(Print)
	{
		cout<< "~~~~~~~~~~~~~~~~~~~~~~~Data Packing After~~~~~~~~~~~~~~~~~~~~~~~~\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		cout << " quot : "<< quot << " remain : " << remain << endl;
		for( i = 0 ; i < cnt ; i++ )
		{
			cout << i+1 << " ciper : "<< R_array_result[i] + V_array_result[i] <<endl;
		}
	}
*/	
	if(!MAXn_flag){
		MAXn_flag = true;
		makeGate();  // 매번해야하는지 ? 체크
	}
	MAXn_result_idx = 0;
	MAXn_tmp_val = V_array_result[0];
	MAXn_tmp_rand = R_array_result[0];
	
	for( i = 0 ; i < cnt ; i++ )
	{
		//cout << MAXn_tmp_val+MAXn_tmp_rand << " ";
		if(G_CMP(V_array_result[i], R_array_result[i], MAXn_tmp_val, MAXn_tmp_rand)){
		//	cout << "<";
			MAXn_tmp_val    = V_array_result[i];	
			MAXn_tmp_rand   = R_array_result[i];
			MAXn_result_idx = i;
		}else{
		//	cout << ">";
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
		cout << "MAXn Value : " << MAXn_tmp_val+MAXn_tmp_rand <<"   idx : "<<MAXn_result_idx;
		gmp_printf("  %Zd\n",paillier_dec(0, pub, prv, ciper[MAXn_result_idx]));	
	}
	
	return MAXn_result_idx;
}

paillier_ciphertext_t** protocol::Topk_sub(paillier_ciphertext_t** ciper_n, int n){
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
int ** protocol::STopk_B(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, int NumData){
	int i = 0;
	int s = 0;
	int n = 0;
	int t = 0;
	int j = 0;
	int rand = 5;
	n = NumData;
	paillier_ciphertext_t * hint = paillier_create_enc_zero();

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* cipher_binary;
	paillier_ciphertext_t* cipher_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_dist	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_dist = paillier_create_enc_zero();

	paillier_ciphertext_t** cipher_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t*** cipher_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*n);
	paillier_ciphertext_t** cipher_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** cipher_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** cipher_V=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t*** cipher_V2 = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*n);

	paillier_ciphertext_t*** cipher_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	for( i = 0 ; i < size; i++ ){
		cipher_Smin[i] = cipher_zero;
	}

	for( i = 0 ; i < NumData ; i++ ){
		cipher_distance[i] 	= cipher_zero;
		cipher_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		cipher_V[i]	= cipher_zero;
		cipher_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	}

	for( i = 0 ; i < k ; i++ ){
		cipher_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim; j ++ ){
			cipher_result[i][j] = paillier_create_enc_zero();
		}
	}

	startTime = std::chrono::system_clock::now(); // startTime check
	printf("\n=== Data compute & SBD start ===\n");
	for( i = 0 ; i < NumData ; i++ ){
		cipher_distance[i] = computeScore(q, data[i]);
		cipher_SBD_distance[i] = SBD(cipher_distance[i]);
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["sTopk_B_SBD(CS(score))"] = time_variable.find("sTopk_B_SBD(CS(score))")->second + duration_sec.count();
	
	for( s = 0 ; s < k ; s++ ){
		printf("\n%dth sMAX start \n", s+1);

		startTime = std::chrono::system_clock::now(); // startTime check
		cipher_Smin = Smax_n(cipher_SBD_distance, NumData);	// bit로 표현된 암호화 min 거리 추출
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["sTopk_B_smax(score))"] = time_variable.find("sTopk_B_smax(score))")->second + duration_sec.count();

		startTime = std::chrono::system_clock::now(); // startTime check
		// bit로 표현된 암호화 max 거리를 암호화 정수로 변환
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			cipher_binary = paillier_create_enc(t);
			cipher_binary = SM_p1(cipher_binary, cipher_Smin[size-j]);
			paillier_mul(pubkey, cipher_min, cipher_binary, cipher_min);
			paillier_freeciphertext(cipher_binary);
			//paillier_print("min : ", cipher_min);
		}
		paillier_print("max dist : ", cipher_min);
		
		printf("\n== recalculate query<->data distances ===\n");
		// 이전 iteration에서 min값으로 선택된 데이터의 거리가 secure하게 MAX로 변환되었기 때문에, 
		// 질의-데이터 간 거리 계산을 모두 다시 수행
		if( s != 0 ){
			for( i = 0 ; i < NumData ; i++ ){
				for( j = size ; j > 0 ; j--){
					t = (int)pow(2, j-1);
					cipher_binary = paillier_create_enc(t);
					cipher_binary = SM_p1(cipher_binary, cipher_SBD_distance[i][size-j]);
					paillier_mul(pubkey, cipher_dist, cipher_binary, cipher_dist);
				}
				cipher_distance[i] = cipher_dist;
				cipher_dist = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
			}
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["sTopk_B_recover_MAX_Score"] = time_variable.find("sTopk_B_recover_MAX_Score")->second + duration_sec.count();
		

		startTime = std::chrono::system_clock::now(); // startTime check
		// 질의-데이터 거리와 min 거리와의 차를 구함 (min 데이터의 경우에만 0으로 만들기 위함)
		for( i = 0 ; i < NumData ; i++ ){
			paillier_subtract(pubkey, temp_dist, cipher_distance[i], cipher_min);
			cipher_mid[i] = SM_p1(temp_dist, cipher_rand);
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["sTopk_B_subtract_Score(SBD_score))"] = time_variable.find("sTopk_B_subtract_Score(SBD_score))")->second + duration_sec.count();



		cipher_V = Topk_sub(cipher_mid, NumData);

		startTime = std::chrono::system_clock::now(); // startTime check
		// min 데이터 추출
		for( i = 0 ; i < NumData ; i++ ){
			for( j = 0 ; j < dim; j++ ){
				cipher_V2[i][j] = SM_p1(cipher_V[i], data[i][j]);
				paillier_mul(pubkey, cipher_result[s][j], cipher_V2[i][j], cipher_result[s][j]);
			}
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["sTopk_B_extract_result"] = time_variable.find("sTopk_B_extract_result")->second + duration_sec.count();


		printf("cipher_result : ");	
		for( j = 0 ; j < dim; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey,cipher_result[s][j]));
		}
		printf("\n");		

		startTime = std::chrono::system_clock::now(); // startTime check
		// Data SBOR 수행
		for( i = 0 ; i < NumData ; i++ ){
			for( j = 0; j < size ; j++ ){
				cipher_SBD_distance[i][j] = SM_p1(SBN(cipher_V[i]), cipher_SBD_distance[i][j]);
			}
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["sTopk_B_update_score"] = time_variable.find("sTopk_B_update_score")->second + duration_sec.count();

		cipher_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
	}
	
	// user(Bob)에게 결과 전송을 위해 random 값 삽입
	for(i=0;i<k;i++){
		for(j=0; j<dim; j++){
			paillier_mul(pubkey, cipher_result[i][j], cipher_result[i][j], cipher_rand);
		}
	}

	paillier_freeciphertext(temp_dist);	
	
	return SkNNm_Bob2(cipher_result, rand, k, dim);
}