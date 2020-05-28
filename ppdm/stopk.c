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

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	int** node_group;

	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);
	if(*NumNodeGroup == 0)
		return 0;


	for(i=0; i<*NumNodeGroup; i++) {
		printf("%dth Node Group : ", i+1);
		for(j=1; j<=node_group[i][0]; j++) {	// 0������ �ش� ��� �׷쿡 ��� ��尡 �ִ����� ����Ǿ� ����
			printf("%d ", node_group[i][j]);
		}
		printf("\n");
	}

	int nodeId = 0;
	int dataId = 0;
	int z = 0;
	int remained = 0;	  // ��� �׷� ������ ���� ó���� �����Ͱ� �����ִ� ��尡 ����� ������

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
	
	// ��� �׷� ���� ������ ������ ����, ���� ������ �����ϴ� ��� �� ������ ����
	for(i=0; i<*NumNodeGroup; i++) {
		remained = node_group[i][0];	 // �ش� ��� �׷� ������ ���� ó���� �����Ͱ� ��?��ִ?��尡 ����� ������
		for(j=0; j<FanOut; j++) {
			cout << j <<endl;
			if(remained == 0)	// �ش� ��� �׷� ������ �� �̻� ó���� �����Ͱ� ���ٸ�, ���� ���� �Ѿ
				break;

			if( (i*FanOut + j) / (float)(*NumNodeGroup*FanOut) >= progress)	{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			for(m=1; m<=node_group[i][0]; m++) {	 // 0������ �ش� ��� �׷쿡 ��� ��尡 �ִ����� ����Ǿ� ����
				nodeId = node_group[i][m];	 // ��� �׷쿡�� ��� ID�� �ϳ��� ����
				if(node[nodeId].NumData >= j+1) 
				{
					dataId = node[nodeId].indata_id[j];	// �ش� ��忡 ����� ������ ID�� �ϳ��� ����

					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(data[dataId][z], alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
						} else {
							tmp[z] = SM_p1(data[dataId][z], alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}

					if(node[nodeId].NumData == j+1)		// �ش� ��尡 ������ �����͸� ó���Ѵٸ�, remained�� 1 ���ҽ�Ŵ
						remained--;
				}
				else {		// �ش� ��忡�� �����Ͱ� ������, ���� ��� �׷� �� �ٸ� ��忡�� ��??ó���� �����Ͱ� �ִ� ��츦 �ڵ鸵
					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(ciper_zero, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
						} else {
							tmp[z] = SM_p1(ciper_zero, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}

				}
			}
			(*cnt)++;	// ��� �׷��� ������ �ѹ��� ������, ������ �ϳ��� �ϼ���
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
	printf("111set_num : %d\n", *NumNodeGroup);

	if(*NumNodeGroup == 0)
		return 0;

	for(i=0; i<*NumNodeGroup; i++) {
		printf("%dth Node Group : ", i);
		for(j=1; j<=node_group[i][0]; j++) {	// 0������ �ش� ��� �׷쿡 ��� ��尡 �ִ����� ����Ǿ� ����
			printf("%d ", node_group[i][j]);
		}
		printf("\n");
	}
	int nodeId = 0;
	int dataId = 0;
	int z = 0;
	int remained = 0;	  // ��� �׷� ������ ���� ó���� �����Ͱ� �����ִ� ��尡 ����� ������

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

	// ��� �׷� ���� ������ ������ ����, ���� ������ �����ϴ� ��� �� ������ ����
	for(i=0; i<*NumNodeGroup; i++) {
		remained = node_group[i][0];	 // �ش� ��� �׷� ������ ���� ó���� �����Ͱ� �����ִ� ��尡 ����� ������

		for(j=0; j<FanOut; j++) {
			//cout << "extract data s" << endl; 
			if(remained == 0)	// �ش� ��� �׷� ������ �� �̻� ó���� �����Ͱ� ���ٸ�, ���� ���� �Ѿ
				break;

			if( (i*FanOut + j) / (float)(*NumNodeGroup*FanOut) >= progress)	{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			for(m=1; m<=node_group[i][0]; m++) {	 // 0������ �ش� ��� �׷쿡 ��� ��尡 �ִ����� ����Ǿ� ����
				nodeId = node_group[i][m];	 // ��� �׷쿡�� ��� ID�� �ϳ��� ����
				if (nodeId == -1){
					if(node[nodeId].NumData == j+1)		// �ش� ��尡 ������ �����͸� ó���Ѵٸ�, remained�� 1 ���ҽ�Ŵ
						remained--;
					continue;
				}
				
				//printf("selected node ID : %d\n", nodeId);

				if(node[nodeId].NumData >= j+1) 
				{
					dataId = node[nodeId].indata_id[j];	// �ش� ��忡 ����� ������ ID�� �ϳ��� ����
					//printf("selected data ID : %d\n", dataId);

					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(data[dataId][z], alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
						} else {
							tmp[z] = SM_p1(data[dataId][z], alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}

					if(node[nodeId].NumData == j+1)		// �ش� ��尡 ������ �����͸� ó���Ѵٸ�, remained�� 1 ���ҽ�Ŵ
						remained--;
				}
				else {		// �ش� ��忡�� �����Ͱ� ������, ���� ��� �׷� �� �ٸ� ��忡�� ���� ó���� �����Ͱ� �ִ� ��츦 �ڵ鸵
					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(ciper_zero, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
						} else {
							tmp[z] = SM_p1(ciper_zero, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}
				}
			}
			(*cnt)++;	// ��� �׷��� ������ �ѹ��� ������, ������ �ϳ��� �ϼ���
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
				copy_ciper[e] = Smax_basic1(copy_ciper[e],copy_ciper[q]);
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
	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;
	time_t startTime = 0; time_t endTime = 0; float gap = 0.0;
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
	paillier_ciphertext_t* ciper_binary;
	paillier_ciphertext_t* ciper_max	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_scores	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_score = paillier_create_enc_zero();
	paillier_ciphertext_t* kth_score = 0;

	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	paillier_ciphertext_t* temp_alpha;
	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);	 // ������
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
		ciper_coeff_bit[i] = SBD_for_SRO(coeff[i], 0);			// query(coeff) ��Ʈ ��ȯ ����
		psi[i] = SCMP(ciper_hint_bit, ciper_coeff_bit[i]);		// ������ ��� (1�̸� ����� ����� ��, 0�̸� ����� ������ ��)
		// debugging
		if(Print)
		{
			printf("coeff (SBD) : ");
			for(j=0; j<size+1; j++) 
			{
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_coeff_bit[i][j]));
			}
			printf("\n");
			gmp_printf("psi : %Zd ", paillier_dec(0, pubkey, prvkey, psi[i]));
			if ( strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, psi[i]) ), paillier_plaintext_to_str(plain_one)) == 0 )
			{
				cout << " Q["<< i <<"] > 0 " <<endl; 
			}else{
				cout << " Q["<< i <<"] < 0 " <<endl; 
			}
		}
	}
	startTime = clock();
	// �� ��� ��Ʈ ��ȯ ����
	cout << "\n\nnode.LL, node.RR SBD Start" <<endl;
	for(i=0; i<NumNode; i++) {	
		for(j=0; j<dim; j++) {
			ciper_nodeLL_bit[i][j] = SBD_for_SRO(node[i].LL[j], 0);				// node LL bound ��ȯ
			ciper_nodeRR_bit[i][j] = SBD_for_SRO(node[i].RR[j], 1);	 			// node RR bound ��ȯ
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
	endTime = clock();
	node_Processing_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("node SBD time: %f\n", node_Processing_time);
	
	for(i=0; i<dim; i++) {	
		temp_q[i] = SM_p1(max_val[i], psi[i]);	
	}
	for(i=0; i<dim; i++) {
		ciper_qLL_bit[i] = SBD_for_SRO(temp_q[i], 0);		// query LL bound ��ȯ
		ciper_qRR_bit[i] = SBD_for_SRO(temp_q[i], 1);		// query RR bound ��ȯ
	}
	// debugging
	if(Print)
	{
		printf("===query for node searching===\n");
		cout << "Max_val : ";
		for(i=0; i<dim; i++) {
			gmp_printf("%Zd\t\t", paillier_dec(0, pubkey, prvkey, max_val[i]));
		}cout<<endl;
		for(i=0; i<dim; i++) {
			gmp_printf("%Zd\t\t", paillier_dec(0, pubkey, prvkey, temp_q[i]));
		}
		cout<<endl;
		for(i=0; i<dim; i++) {
			for(m=0; m<size+1; m++) {
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_qLL_bit[i][m]));
			}
			cout <<" ";
		}
		printf("\n");
	}
	cout << "FIND Range\tFIND Range\tFIND Range\tFIND Range\tFIND Range\t"<<endl;
	startTime = clock();
	for(i=0; i<NumNode; i++) {	
		alpha[i] = SRO(ciper_qLL_bit, ciper_qRR_bit, ciper_nodeLL_bit[i], ciper_nodeRR_bit[i]);		
		if(!verify_flag)	{	//  ���� �ܰ迡���� �������� ����
			for(j=0; j<dim; j++) {
				node[i].LL[j] = SM_p1(node[i].LL[j], SBN(alpha[i]));
				node[i].RR[j] = SM_p1(node[i].RR[j], SBN(alpha[i]));
			}
		}
		if(1)
		{
			gmp_printf("alpha %d : %Zd  ", i+1 , paillier_dec(0, pubkey, prvkey, alpha[i]));
			if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) {
				cout << " OverLaps "<<endl;
			}else{
				cout << " Not overlaps "<<endl;
			}
		}
	}
	endTime = clock();
	node_SRO_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("SRO time : %f\n", node_SRO_time);
	printf("test  (%d line)\n", __LINE__);
	while(1) {
		cnt = 0;
		if(!verify_flag)	{	// ���� �ܰ谡 �ƴ� �ÿ���, cand�� ����
			startTime = clock();
			cand	= sNodeRetrievalforTopk(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_first_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("data_extract_first_time : %f\n", data_extract_first_time);
	
			if(Print)
			{
				cout << "cnt : " <<cnt <<endl;
				cout << "######cand######"<<endl;
				for( i = 0 ; i < cnt ; i ++){
					cout << i+1 <<" : ";
					for( j = 0 ; j < dim ; j++){
						gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
					}
					cout<<endl;
				}
			}
		}else {		// ���� �ܰ��� �ÿ���, ���� ����� cand�� ��ħ
			startTime = clock();
			temp_cand = sNodeRetrievalforTopk(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_second_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("data_extract_first_time : %f\n", data_extract_second_time);
	
			printf("cnt : %d \n", cnt);
			if(cnt == 0) {	// ������ ���� �߰� Ž���� �ʿ��� ��尡 ���� ��츦 ó����
				break;
			}
			
			cand =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ ){
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = paillier_create_enc_zero();
				}
			}
			// ���� k���� ����� ����
			for(i=0; i<k; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = ciper_result[i][j];
				}
			}
			// ���� ã�� �ĺ� ����� ����
			for(i=0; i<cnt; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i+k][j] = temp_cand[i][j];
				}
			}
			cnt = cnt + k;
			if(Print)
			{
				cout << "######ciper_result + TEMP_cand######"<<endl;
				for( i = 0 ; i < cnt ; i ++){
					cout << i+1 <<" : ";
					for( j = 0 ; j < dim ; j++){
						gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
					}
					cout<<endl;
				}
			}
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
			ciper_V[i]			= ciper_zero;
			ciper_V2[i]			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		}
		for( i = 0 ; i < size; i++ ){
			ciper_Smax[i] = ciper_zero;
		}
		/*
		if(1)
		{
			cout << "ciper_Smax : "; 
			for( j = 0 ; j < dim ; j++){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smax[j]));
			}cout<<endl;
		}
		*/
		startTime = clock();
		for( i = 0 ; i < cnt ; i++ ){
			//ciper_score[i] = computeScore2(q, cand[i], coeff, hint);
			ciper_score[i] = computeScore(q, cand[i]);
			ciper_SBD_score[i] = SBD(ciper_score[i]);
			/*
			if(Print)
			{
				gmp_printf("score : %Zd\n",  paillier_dec(0, pubkey, prvkey, ciper_score[i]));	
				for( j = 0 ; j < size ; j++ ){
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_SBD_score[i][j]));
				}
				printf("\n");
			}
			*/
		}
		endTime = clock();
		data_SSED_SBD_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		printf("SSED & SBD time : %f\n", data_SSED_SBD_time);

		for( s = 0 ; s < k ; s++ ){
			printf("\n%dth sTopk start \n", s+1);
			startTime = clock();
			ciper_Smax = Smax_n(ciper_SBD_score, cnt);	// bit�� ǥ���� ��ȣȭ max score ����
			endTime = clock();
			if(!verify_flag){
				sMINn_first_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("%dth sMAXn first time : %f\n", s+1, sMINn_first_time);
			}else{
				sMINn_second_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("%dth sMAXn second time : %f\n", s+1, sMINn_second_time);
			}
			/*
			if(Print)
			{
				cout<< "MAX : ";
				for( j = 0 ; j < size ; j++ ){
					gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smax[j]));
				}
				printf("\n");	
			}
			*/
			// bit�� ǥ���� ��ȣȭ max score�� ��ȣȭ ������ ��ȯ
			for( j = size ; j > 0 ; j-- ){
				t = (int)pow(2, j-1);
				ciper_binary = paillier_create_enc(t);
				ciper_binary = SM_p1(ciper_binary, ciper_Smax[size-j]);
				paillier_mul(pubkey, ciper_max, ciper_binary, ciper_max);
				paillier_freeciphertext(ciper_binary);
				//gmp_printf("min : %Zd\n", ciper_max);
			}
			gmp_printf("MAX Value score: %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_max));

			if(Print)
			{
				gmp_printf("MAX Value score: %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_max));
			}

			// ���� iteration���� max������ ���õ� �������� score�� secure�ϰ� MIN���� ��ȯ�Ǿ��� ������, 
			// �������� ���� ����� ��� �ٽ� ����
			if( s != 0 ){
				for( i = 0 ; i < cnt; i++ ){
					temp_score = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
					for( j = size ; j > 0 ; j--){
						t = (int)pow(2, j-1);
						ciper_binary = paillier_create_enc(t);
						ciper_binary = SM_p1(ciper_binary, ciper_SBD_score[i][size-j]);
						paillier_mul(pubkey, temp_score, ciper_binary, temp_score);
					}
					ciper_score[i] = temp_score;
					temp_score = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}
			if(Print)
			{
				for( i = 0 ; i < cnt ; i++ ){
					gmp_printf("score : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_score[i]) );
				}
			}

			// �� �������� ������ max score���� ���� ���� (max �������� ��쿡�� 0���� ����� ����)
			for( i = 0 ; i < cnt ; i++ ){
				paillier_subtract(pubkey, temp_score, ciper_score[i], ciper_max);
				ciper_mid[i] = SM_p1(temp_score, ciper_rand);
				if(Print)
				{
					gmp_printf("dist-dist : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_mid[i]));
				}
			}

			ciper_V = Topk_sub(ciper_mid, cnt);
			
			if(Print)
			{
				for( i = 0 ; i < cnt ; i++ ){
					gmp_printf("%d ciper_V : %Zd\n", i+1, paillier_dec(0, pubkey, prvkey, ciper_V[i]));
				}
			}

			// max ������ ����
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
			startTime = clock();
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0; j < size ; j++ ){
					ciper_SBD_score[i][j] = SM_p1(SBN(ciper_V[i]), ciper_SBD_score[i][j]);
				}
			}
			ciper_max = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);
			endTime = clock();
			data_SBOR_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			if(Print)
			{
				for( i = 0 ; i < cnt ; i++ ){
					cout<<"cand : ";
					for( j = 0 ; j < dim; j++ ){
						gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
					}
					cout<<endl;
				}
				printf("result : ");
				for( j = 0 ; j < dim; j++ ){
					gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, ciper_result[s][j]));
				}
				printf("\n");
				cout<<"SBD score"<<endl;
				for( i = 0 ; i < cnt ; i++ ){
					cout<< i+1 <<" : ";
					for( j = 0; j < size ; j++ ){
						gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_SBD_score[i][j]));
					}
					cout<<endl;
				}
				cout<<endl;
			}
		}
		if(Print)
		{
			cout<< "!!!!!!!!!!!!!MAX !!!!!!!!!!!!!"<<endl;
			for( s = 0 ; s < k ; s ++){
				cout << s << " : ";
				for( i = 0 ; i < dim ; i ++){
					gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_result[s][i]));
				}
				cout<<endl;
			}
		}
		if(!verify_flag){
			startTime = clock();
			if(cnt < k)	 // ����� fanout�� k���� ���� ��츦 ó����
				kth_score = computeScore(q, ciper_result[cnt-1]);
			else	// k���� ã�����ٸ�, k��° �������� score�� ���
				kth_score = computeScore(q, ciper_result[k-1]);
			gmp_printf("%dth score : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_score));
			for(i=0; i<NumNode; i++) {	
				for(j=0; j<dim; j++) {
					paillier_mul(pubkey, temp_q[j], SM_p1(node[i].RR[j], psi[j]), SM_p1(node[i].LL[j], SBN(psi[j])));	
				}
				ciper_kthscore_bit = SBD_for_SRO(kth_score, 0);		// k��° score SBD ����
				//ciper_kthscore_bit = SBD(kth_score);		// k��° score SBD ����
				ciper_tempscore_bit = SBD_for_SRO(computeScore(q, temp_q), 1);  // �� ����� score�� SBD ����
				//ciper_tempscore_bit = SBD(computeScore(q, temp_q));  // �� ����� score�� SBD ����
				alpha[i] = SCMP(ciper_kthscore_bit , ciper_tempscore_bit);	// k��° score���� ���� ���ɼ��� ������ alpha=1
				//alpha[i] = SCMP_for_SBD(ciper_kthscore_bit , ciper_tempscore_bit);	// k��° score���� ���� ���ɼ��� ������ alpha=1
	
				if(1){
					gmp_printf("computeScore(q, temp_q) : %Zd %d line\n", paillier_dec(0, pubkey, prvkey, computeScore(q, temp_q)), __LINE__);
					gmp_printf("alpha : %Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
				}
			}
			endTime = clock();
			node_expansion_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		}
		if(!verify_flag)	{
			verify_flag = true;		// ������ �����ϱ� ���� flag�� true�� ����
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
	free(ciper_max);
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
	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��
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

	paillier_ciphertext_t* ciper_binary;
	paillier_ciphertext_t* ciper_max	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_scores	= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_score = paillier_create_enc_zero();
	paillier_ciphertext_t* kth_score = 0;

	paillier_ciphertext_t*** ciper_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	paillier_ciphertext_t* temp_alpha;
	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);	 // ������
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

	// query�� �� ���� �����Ϳ� hint�� ����
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

	ciper_hint_bit = SBD_for_SRO(hint, 0);		// hint SBD ����
	
	// debugging
	printf("hint (SBD) : ");
	for(j=0; j<size+1; j++) {
		gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_hint_bit[j]));
	}
	printf("\n");
	
	for(i=0; i<dim; i++) {
		ciper_coeff_bit [i] = SBD_for_SRO(coeff[i], 0);			// query(coeff) ��Ʈ ��ȯ ����
		
		// debugging
		printf("coeff (SBD) : ");
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_coeff_bit [i][j]));
		}
		printf("\n");

		psi[i] = SCMP(ciper_hint_bit, ciper_coeff_bit[i]);		// ������ ��� (1�̸� ����� ����� ��, 0�̸� ����� ������ ��)
	}
	
	for(i=0; i<dim; i++) {
		gmp_printf("psi : %Zd\n", paillier_dec(0, pubkey, prvkey, psi[i]));
	}

	startTime = clock();

	// �� ��� ��Ʈ ��ȯ ���� 
	for(i=0; i<NumNode; i++) {	
		for(j=0; j<dim; j++) {
			ciper_nodeLL_bit[i][j] = SBD_for_SRO(node[i].LL[j], 0);				// node LL bound ��ȯ
			ciper_nodeRR_bit[i][j] = SBD_for_SRO(node[i].RR[j], 1);	 			// node RR bound ��ȯ
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

	// ù Ž�� ��带 ã�� ���� ���� ����(��� ����� Max, ���� ����� 0)
	for(i=0; i<dim; i++) {	
		temp_q[i] = SM_p1(max_val[i], psi[i]);	
	}

	// debugging
	printf("===query for node searching===\n");
	for(i=0; i<dim; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, temp_q[i]));
	}
	printf("\n");

	// query ��Ʈ ��ȯ ����
	for(i=0; i<dim; i++) {
		ciper_qLL_bit[i] = SBD_for_SRO(temp_q[i], 0);		// query LL bound ��ȯ
		ciper_qRR_bit[i] = SBD_for_SRO(temp_q[i], 1);		// query RR bound ��ȯ
	}
	


	// score�� ���� ���� ������ ����Ǵ� ��� Ž��
	// �� ��� SRO ȣ��
	startTime = clock();
	for(i=0; i<NumNode; i++) {	
		alpha[i] = SRO(ciper_qLL_bit, ciper_qRR_bit, ciper_nodeLL_bit[i], ciper_nodeRR_bit[i]);		

		if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
			printf("%dth node overlaps the query region.\n", i);

		/*
		if(!verify_flag)	{	//  ���� �ܰ迡���� �������� ����
			// �̹� �˻��� �Ϸ�� ��尡 �� Ž���Ǵ� ���� �����ϱ� ����, bound�� MAX�� ��ȯ
			for(j=0; j<dim; j++) {
				for(m=0; m<size+1; m++) {
					ciper_nodeLL_bit[i][j][m] = SBOR(ciper_nodeLL_bit[i][j][m], alpha[i]);
					ciper_nodeRR_bit[i][j][m] = SBOR(ciper_nodeRR_bit[i][j][m], alpha[i]);
				}
			}
		}
		*/

		if(!verify_flag)	{	//  ���� �ܰ迡���� �������� ����
			// �̹� �˻��� �Ϸ�� ��尡 �� Ž���Ǵ� ���� �����ϱ� ����, bound�� 0���� ��ȯ
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

		if(!verify_flag)	{	// ���� �ܰ谡 �ƴ� �ÿ���, cand�� ����
			cand	= sNodeRetrievalforTopk(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
		}
		else {		// ���� �ܰ��� �ÿ���, ���� ����� cand�� ��ħ
			temp_cand	= sNodeRetrievalforTopk(data, ciper_qLL_bit, ciper_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);

			printf("cnt : %d \n", cnt);
			if(cnt == 0) {	// ������ ���� �߰� Ž���� �ʿ��� ��尡 ���� ��츦 ó����
				break;
			}

			cand =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ ){
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = paillier_create_enc_zero();
				}
			}

			// ���� k���� ����� ����
			for(i=0; i<k; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = ciper_result[i][j];
				}
			}

			// ���� ã�� �ĺ� ����� ����
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
			ciper_Smax[i] = ciper_zero;
		}

		for( i = 0 ; i < cnt ; i++ ){
			ciper_score[i] 	= 	(paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(alpha[i]->c);

			ciper_SBD_score[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
			ciper_V[i]	= ciper_zero;
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
			ciper_Smax = Smax_n(ciper_SBD_score, cnt);	// bit�� ǥ���� ��ȣȭ max score ����
			endTime = clock();
			gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("%dth sMAXn time : %f\n", s+1, gap);
		
			
			for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smax[j]));
			}
			printf("\n");	

			// bit�� ǥ���� ��ȣȭ max score�� ��ȣȭ ������ ��ȯ
			for( j = size ; j > 0 ; j-- ){
				t = (int)pow(2, j-1);
				ciper_binary = paillier_create_enc(t);
				ciper_binary = SM_p1(ciper_binary, ciper_Smax[size-j]);
				paillier_mul(pubkey, ciper_max, ciper_binary, ciper_max);
				paillier_freeciphertext(ciper_binary);
				//gmp_printf("min : %Zd\n", ciper_max);
			}
			gmp_printf("max score: %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_max));


			// ���� iteration���� max������ ���õ� �������� score�� secure�ϰ� MIN���� ��ȯ�Ǿ��� ������, 
			// �������� ���� ����� ��� �ٽ� ����
			if( s != 0 ){
				for( i = 0 ; i < cnt; i++ ){
					for( j = size ; j > 0 ; j--){
						t = (int)pow(2, j-1);
						ciper_binary = paillier_create_enc(t);
						ciper_binary = SM_p1(ciper_binary, ciper_SBD_score[i][size-j]);
						paillier_mul(pubkey, temp_score, ciper_binary, temp_score);
					}
					ciper_score[i] = temp_score;
					temp_score = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}

			for( i = 0 ; i < cnt ; i++ ){
				gmp_printf("score : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_score[i]) );
			}

			// �� �������� ������ max score���� ���� ���� (max �������� ��쿡�� 0���� ����� ����)
			for( i = 0 ; i < cnt ; i++ ){
				paillier_subtract(pubkey, temp_score, ciper_score[i], ciper_max);
				ciper_mid[i] = SM_p1(temp_score, ciper_rand);
				//gmp_printf("dist - dist : ",ciper_mid[i]);
			}

			ciper_V = SkNNm_sub(ciper_mid, cnt);
			
			
			for( i = 0 ; i < cnt ; i++ ){
				gmp_printf("%dth ciper_V : %Zd\n", i, paillier_dec(0, pubkey, prvkey, ciper_V[i]));
			}

			// max ������ ����
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

			ciper_max = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
		}


		if(!verify_flag)	{	
			if(cnt < k)	 // ����� fanout�� k���� ���� ��츦 ó����
				kth_score = computeScore(q, ciper_result[cnt-1]);
			else	// k���� ã�����ٸ�, k��° �������� score�� ���
				kth_score = computeScore(q, ciper_result[k-1]);

			gmp_printf("%dth score : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_score));

			
			for(i=0; i<NumNode; i++) {	
				for(j=0; j<dim; j++) {
					printf("test\n");
					// Ȯ�� Ž�� ��带 ã�� ���� ���� ����(��� ����� �ش� ����� UB, ���� ����� �ش� ����� LB)
					gmp_printf("%d line: %Zd\n", __LINE__, paillier_dec(0, pubkey, prvkey, SM_p1(node[i].RR[j], psi[j])));
					gmp_printf("%d line: %Zd\n", __LINE__, paillier_dec(0, pubkey, prvkey, SM_p1(node[i].LL[j], SBN(psi[j]))));

					paillier_mul(pubkey, temp_q[j], SM_p1(node[i].RR[j], psi[j]), SM_p1(node[i].LL[j], SBN(psi[j])));	
					gmp_printf("%d line: %Zd\n	", __LINE__, paillier_dec(0, pubkey, prvkey, temp_q[j]));
				}

				ciper_kthscore_bit = SBD_for_SRO(kth_score, 0);		// k��° score SBD ����
				//ciper_kthscore_bit = SBD(kth_score);		// k��° score SBD ����

				gmp_printf("computeScore(q, temp_q) : %Zd %d line\n", paillier_dec(0, pubkey, prvkey, computeScore(q, temp_q)), __LINE__);
				ciper_tempscore_bit = SBD_for_SRO(computeScore(q, temp_q), 1);  // �� ����� score�� SBD ��??
				//ciper_tempscore_bit = SBD(computeScore(q, temp_q));  // �� ����� score�� SBD ����


				alpha[i] = SCMP(ciper_kthscore_bit , ciper_tempscore_bit);	// k��° score���� ���� ���ɼ��� ������ alpha=1

				gmp_printf("alpha : %Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
			}
		}

		if(!verify_flag)	{
			verify_flag = true;		// ������ �����ϱ� ���� flag�� true�� ����
		}
		else
			break;	// ������ ����
	}
		
	// user(Bob)���� ��� ������ ���� random �� ����
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
	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;
	time_t startTime = 0;	time_t endTime = 0;	float gap = 0.0;

	//���� �� hint ����
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
	paillier_ciphertext_t** psi				= (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*dim);	 // ������
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
	
	
	hint = q[dim];
	if(Print)
	{
		gmp_printf("hint : %Zd  /  %Zd\n", paillier_dec(0, pubkey, prvkey, hint),paillier_dec(0, pubkey, prvkey, q[dim]));
	}
	for( i = 0 ; i < dim ; i++ )
	{
		paillier_mul(pubkey, coeff[i], q[i], hint);
		if(Print)
		{
			cout <<"q[i] + hint : "<<endl;
			gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, coeff[i]));
			gmp_printf("max val : %Zd \n", paillier_dec(0, pub, prv, max_val[i]));
		}
	}
	for( i = 0 ; i < dim ; i++ )
	{
		psi[i] = GSCMP(hint, coeff[i]);
		if(Print)
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
	if(Print)
	{
		printf("===query for node searching===\n");
		cout << "Max_val\n";
		for( i = 0 ; i < dim ; i++ )
		{
			gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, max_val[i]));
		}cout<<endl;
		cout << "Q : "<<endl;
		for( i = 0 ; i < dim ; i++ ) 
		{
			gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, Q[i]));
		}
		cout<<endl;
	}
	cout << "FIND Range\tFIND Range\tFIND Range\tFIND Range\tFIND Range\t"<<endl;
	startTime = clock();
	for( i = 0 ; i < NumNode ; i++ )
	{	
		alpha[i] = DP_GSRO(Q, Q, node[i].LL, node[i].RR);
		
		if(!verify_flag)
		{	//  ���� �ܰ迡���� �������� ����
			for( j = 0 ; j < dim ; j++ )
			{
				node[i].LL[j] = SM_p1(node[i].LL[j], SBN(alpha[i]));
				node[i].RR[j] = SM_p1(node[i].RR[j], SBN(alpha[i]));
			}
		}
		if(Print)
		{
			gmp_printf("alpha %d : %Zd  ", i+1 , paillier_dec(0, pubkey, prvkey, alpha[i]));
			if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) {
				cout << " OverLaps "<<endl;
			}else{
				cout << " Not overlaps "<<endl;
			}
		}
	}
	endTime = clock();
	node_SRO_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("GSRO time : %f\n", node_SRO_time);
	while(1)
	{
		cnt = 0;
		if(!verify_flag)
		{	
			startTime = clock();
			cand = GSRO_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_first_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_first_time);
			
			if(Print)
			{
				cout << "cnt : " <<cnt <<endl;
				cout << "######cand######"<<endl;
				for( i = 0 ; i < cnt ; i++ )
				{
					cout << i+1 <<" : ";
					for( j = 0 ; j < dim ; j++ )
					{
						gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
					}
					cout<<endl;
				}
			}
		}
		else
		{
			startTime = clock();
			temp_cand = GSRO_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_second_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node expansion time : %f\n", data_extract_second_time);
			printf("cnt : %d \n", cnt);
			if(cnt == 0) {	// ������ ���� �߰� Ž���� �ʿ��� ��尡 ���� ��츦 ó����
				break;
			}
			cand =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ ){
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = paillier_create_enc_zero();
				}
			}
			// ���� k���� ����� ����
			for(i=0; i<k; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = RESULT[i][j];
				}
			}
			// ���� ã�� �ĺ� ����� ����
			for(i=0; i<cnt; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i+k][j] = temp_cand[i][j];
				}
			}
			cnt = cnt + k;
			if(Print)
			{
				cout << "######ciper_result + TEMP_cand######"<<endl;
				for( i = 0 ; i < cnt ; i ++){
					cout << i+1 <<" : ";
					for( j = 0 ; j < dim ; j++){
						gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
					}
					cout<<endl;
				}
			}
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

		startTime = clock();
		for( i = 0 ; i < cnt ; i++ )
		{
			SCORE[i] = computeScore(q, cand[i]);
			if(Print)
			{
				gmp_printf("SCORE : %Zd\n",  paillier_dec(0, pubkey, prvkey, SCORE[i]));	
			}
		}
		endTime = clock();
		data_SSED_SBD_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		printf("ComputeScore & SBD time : %f\n", data_SSED_SBD_time);

		for( s = 0 ; s < k ; s++ )
		{
			cout << s+1<<"dth MAXn_Topk start"<<endl;
			startTime = clock();
			MAX_idx = MAXn(SCORE, cnt);
			MAX = SCORE[MAX_idx];
			if(Print)
			{
				gmp_printf("MAX VALUE : %Zd\n",  paillier_dec(0, pubkey, prvkey, MAX));	
			}
			endTime = clock();
			if(!verify_flag) {
				sMINn_first_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("ComputeScore & SBD time : %f\n", sMINn_first_time);
			}
			else {
				sMINn_second_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
				printf("ComputeScore & SBD time : %f\n", sMINn_second_time);
			}
			
			for( i = 0 ; i < cnt ; i++ )
			{
				paillier_subtract(pubkey, SCORE_MINUS_MAX[i], SCORE[i], MAX);
				SCORE_MINUS_MAX[i] = SM_p1(SCORE_MINUS_MAX[i], C_RAND);
				if(Print)
				{
					gmp_printf("SCORE-MAX : %Zd \n", paillier_dec(0, pubkey, prvkey, SCORE_MINUS_MAX[i]));
				}
			}
			
			V = Topk_sub(SCORE_MINUS_MAX, cnt);

			for( i = 0 ; i < cnt ; i++ )
			{
				if(Print)
				{
					gmp_printf("V : %Zd \n", paillier_dec(0, pubkey, prvkey, V[i]));

				}
				paillier_subtract(pubkey, TMP_alpha, ciper_one, V[i]);
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
		}
		if(Print)
		{
			cout<< "!!!!!!!!!!!!!MAXn!!!!!!!!!!!!!"<<endl;
			for( s = 0 ; s < k ; s++ )
			{
				cout << s << " : ";
				for( i = 0 ; i < dim ; i++ )
				{
					gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, RESULT[s][i]));
				}
				cout<<endl;
			}
		}
		if(!verify_flag)
		{	
			startTime = clock();
			if(cnt < k)	 // ����� fanout�� k���� ���� ��츦 ó����
				kth_score = computeScore(q, RESULT[cnt-1]);
			else	// k���� ã�����ٸ�, k��° �������� score�� ���
				kth_score = computeScore(q, RESULT[k-1]);
			//gmp_printf("%dth score : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_score));
			for( i = 0 ; i < NumNode ; i++ )
			{	
				for( j = 0 ; j < dim ; j++ )
				{
					paillier_mul(pubkey, Q[j], SM_p1(node[i].RR[j], psi[j]), SM_p1(node[i].LL[j], SBN(psi[j])));	
				}
				alpha[i] = GSCMP(kth_score , computeScore(q, Q));	// k��° score���� ���� ���ɼ��� ������ alpha=1
				
				if(Print)
				{
					gmp_printf("computeScore(q, Q) : %Zd %d line\n", paillier_dec(0, pubkey, prvkey, computeScore(q, Q)), __LINE__);
					gmp_printf("alpha : %Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
				}
			}
			endTime = clock();
			node_expansion_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
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
		gmp_printf(" %Zd : ", paillier_dec(0, pubkey, prvkey, computeScore(q, RESULT[i])));
	
		for( j = 0 ; j < dim ; j++ )
		{
	//		gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, RESULT[i][j]));
			paillier_mul(pubkey, RESULT[i][j], RESULT[i][j], C_RAND);
		}
	
	//	printf("\n");
	
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
		makeGate();  // �Ź��ؾ��ϴ��� ? üũ
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
/*
paillier_ciphertext_t* protocol::DP_computeScore(paillier_ciphertext_t** cipher1, paillier_ciphertext_t** cipher2, int col_num)
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
	
	paillier_ciphertext_t* mul =  paillier_create_enc_zero();

	paillier_ciphertext_t** result_array=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*col_num);
	for(int i=0;i<col_num;i++){
		result_array[i]	= paillier_create_enc_zero();
	}


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

		paillier_exp(pubkey,mul,cipher1[i],cipher2[i]); 
		gmp_printf("mul : %Zd\n",paillier_dec(0, pubkey, prvkey, mul));

				
		paillier_exp(pubkey, ciper_tmp, mul, shift); //multi		
		gmp_printf("ciper_tmp : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_tmp));

		paillier_mul(pubkey,ciper_value,ciper_value,ciper_tmp); //E_add
		//gmp_printf("value : %Zd\n", paillier_dec(0, pubkey, prvkey, ciper_value));
	}
	
		
	//result_array=unDP_SSED(ciper_value,col_num);
	paillier_ciphertext_t* test = unDP_computeScore(ciper_value,col_num);

	result = unDP_computeScore(ciper_value,col_num);

	//gmp_printf("unDP_SSED result : %Zd\n",paillier_dec(0, pubkey, prvkey, test));

	for(int i=0; i<col_num; i++){
		mpz_mul(tmp->m,rand[i]->m,rand[i]->m);
		ciper_tmp=paillier_enc(0,pub,tmp,paillier_get_rand_devurandom);
		paillier_subtract(pub,result,result,ciper_tmp );

		mpz_mul_ui (rand[i]->m,rand[i]->m,2);
		paillier_subtract(pub,ciper_tmp,cipher1[i],cipher2[i]);
		paillier_exp(pubkey, sub, ciper_tmp, rand[i]);
		paillier_subtract(pub,result,result,sub);
	}

	//gmp_printf("result : %Zd\n",paillier_dec(0, pubkey, prvkey, result));

	return result;	

}

paillier_ciphertext_t* protocol::unDP_computeScore(paillier_ciphertext_t* cipher1, int col_num){
	paillier_ciphertext_t** result_array=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*col_num);

	paillier_plaintext_t* shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* tmp = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* result = paillier_dec(0, pub, prv, cipher1);

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
	
	for(int i=0; i<col_num; i++){
		gmp_printf("unDP_SSED re : %Zd\n",paillier_dec(0, pub, prv, result_array[i]));
	}
	
	return enc_test;
	//return result_array;
}
*/
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

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_ciphertext_t * hint = paillier_create_enc_zero();

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


	printf("\n=== Data compute & SBD start ===\n");
	startTime = clock();
	for( i = 0 ; i < n ; i++ ){
		ciper_distance[i] = computeScore(q, data[i]);
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
	printf("data compute & SBD time : %f\n", (float)(endTime-startTime)/(CLOCKS_PER_SEC));

	for( s = 0 ; s < k ; s++ ){
		printf("\n%dth sMAX start \n", s+1);
		startTime = clock();
		ciper_Smin = Smax_n(ciper_SBD_distance, n);	// bit�� ǥ���� ��ȣȭ min �Ÿ� ����
		endTime = clock();
		gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		sMINn_first_time += gap;
		printf("%dth sMAX time : %f\n", s+1, gap);
	
		/*
		for( j = 0 ; j < size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, ciper_Smin[j]));
		}
		printf("\n");	
		*/
		
		// bit�� ǥ���� ��ȣȭ max �Ÿ��� ��ȣȭ ������ ��ȯ
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			ciper_binary = paillier_create_enc(t);
			ciper_binary = SM_p1(ciper_binary, ciper_Smin[size-j]);
			paillier_mul(pubkey, ciper_min, ciper_binary, ciper_min);
			paillier_freeciphertext(ciper_binary);
			//paillier_print("min : ", ciper_min);
		}
		paillier_print("max dist : ", ciper_min);
		
		printf("\n== recalculate query<->data distances ===\n");
		// ���� iteration���� min������ ���õ� �������� �Ÿ��� secure�ϰ� MAX�� ��ȯ�Ǿ��� ������, 
		// ����-������ �� �Ÿ� ����� ��� �ٽ� ����
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
		
		// ����-������ �Ÿ��� min �Ÿ����� ���� ���� (min �������� ��쿡�� 0���� ����� ����)
		for( i = 0 ; i < n ; i++ ){
			paillier_subtract(pubkey, temp_dist, ciper_distance[i], ciper_min);
			ciper_mid[i] = SM_p1(temp_dist, ciper_rand);
		}

		ciper_V = Topk_sub(ciper_mid, n);

		// min ������ ����
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
		
		// Data SBOR ����
		startTime = clock();
		for( i = 0 ; i < n ; i++ ){
			for( j = 0; j < size ; j++ ){
				ciper_SBD_distance[i][j] = SM_p1(SBN(ciper_V[i]), ciper_SBD_distance[i][j]);
			}
		}
		endTime = clock();
		data_SBOR_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);

		ciper_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
	}
	
	// user(Bob)���� ��� ������ ���� random �� ����
	for(i=0;i<k;i++){
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,ciper_result[i][j],ciper_result[i][j],ciper_rand);
		}
	}

	paillier_freeciphertext(temp_dist);	
	
	return SkNNm_Bob2(ciper_result, rand, k, dim);
}