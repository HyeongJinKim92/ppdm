#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <iostream>
#include "protocol.h"

using namespace std;

int** protocol::sRange_B(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode,int* result_num)
{
	printf("\n=====now sRange_B start=====\n");

	int i=0, j=0, m=0;
	
	int rand = 5;
	int NumNodeGroup = 0;

	paillier_ciphertext_t*** cipher_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** cipher_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);

	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	
	cout << "\n=====================SBD QUERY=====================\n";
	// query ��Ʈ ��ȯ ����
	startTime = std::chrono::system_clock::now(); // startTime check
	for( i = 0 ; i < dim ; i++ ) {
		cipher_qLL_bit[i] = SBD_for_SRO(q.LL[i], 0);		// query LL bound ��ȯ
		cipher_qRR_bit[i] = SBD_for_SRO(q.RR[i], 1);		// query RR bound ��ȯ
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	//time_variable.insert(make_pair("sRange_B_SBD_for_SRO(QUERY)", duration_sec.count() ));
	time_variable["SBD_Q"] = time_variable.find("SBD_Q")->second + duration_sec.count();



	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumData);
	paillier_ciphertext_t*** candLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** candRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	cout << "\n=====================SRO===================\n";

	for( i = 0 ; i < NumData ; i++ ) {
		startTime = std::chrono::system_clock::now(); // startTime check
		for( j = 0 ; j < dim ; j++ ) {
			candLL_bit[j] = SBD_for_SRO(data[i][j], 0);			// query cand ��ȯ
			candRR_bit[j] = SBD_for_SRO(data[i][j], 1);		
		}		
		endTime = std::chrono::system_clock::now();  // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["SBD_D"] = time_variable.find("SBD_D")->second + duration_sec.count();
		startTime = std::chrono::system_clock::now(); // startTime check
		alpha[i] = SRO(candLL_bit, candRR_bit, cipher_qLL_bit, cipher_qRR_bit);
		endTime = std::chrono::system_clock::now();  // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["SRO_D"] = time_variable.find("SRO_D")->second + duration_sec.count();
	}
	cout << endl;
	paillier_ciphertext_t*** result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData); //�Ҵ��� �����ؾ��ϰ�

	//�� Ȯ�� && If ��i = 1, E(t��i)�� result�� ����
	for( i = 0 ; i < NumData ; i++ ){
		if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pub, prv, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
		{
			result[(*result_num)]=data[i]; //result�� ����
			(*result_num)++;
		}	
	}

	// user(Bob)���� ��� ������ ���� random �� ����
	for( i = 0 ; i < (*result_num) ; i++ )
	{
		for( j = 0 ; j < dim ; j++ ) 
		{	
				paillier_mul(pubkey, result[i][j], result[i][j], cipher_rand);
		}
	}
	
	return FsRange_Bob(result, rand, (*result_num), dim);
}

int** protocol::sRange_I(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num)
{
	printf("\n===== Now sRange_I starts =====\n");
	int i=0, j=0, m=0;

	float gap = 0.0;
	int rand = 5;
	int NumNodeGroup = 0;

	paillier_ciphertext_t*** cipher_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** cipher_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);

	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	
	startTime = std::chrono::system_clock::now(); // check startTime
	// query ��Ʈ ��ȯ ����
	for(i=0; i<dim; i++) {
		cipher_qLL_bit[i] = SBD_for_SRO(q.LL[i], 0);			// query LL bound ��ȯ
		cipher_qRR_bit[i] = SBD_for_SRO(q.RR[i], 1);		// query RR bound ��ȯ
	}
	endTime = std::chrono::system_clock::now(); // check endTime
	duration_sec = endTime - startTime;
	time_variable["SBD_Q"] = time_variable.find("SBD_Q")->second + duration_sec.count();

	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��

	paillier_ciphertext_t*** cand ;
	
	cand = sNodeRetrievalforRange(data, cipher_qLL_bit, cipher_qRR_bit, node, NumData, NumNode, &cnt, &NumNodeGroup, 0);
	totalNumOfRetrievedNodes += NumNodeGroup;


	paillier_ciphertext_t*** candLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** candRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);

	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);

	startTime = std::chrono::system_clock::now();
	for( i = 0 ; i < cnt ; i++ ) {
		startTime = std::chrono::system_clock::now();
		for(j=0; j<dim; j++) {
			candLL_bit[j] = SBD_for_SRO(cand[i][j], 0);			// query cand ��ȯ
			candRR_bit[j] = SBD_for_SRO(cand[i][j], 1);		
		}
		endTime = std::chrono::system_clock::now();
		duration_sec = endTime - startTime;
		time_variable["SBD_D"] = time_variable.find("SBD_D")->second + duration_sec.count();

		startTime = std::chrono::system_clock::now();
		alpha[i] = SRO(candLL_bit, candRR_bit, cipher_qLL_bit, cipher_qRR_bit);		
		endTime = std::chrono::system_clock::now();
		duration_sec = endTime - startTime;
		time_variable["SRO_D"] = time_variable.find("SRO_D")->second + duration_sec.count();

	}


	paillier_ciphertext_t*** result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumNodeGroup*FanOut); //�Ҵ��� �����ؾ��ϰ�

	//�� Ȯ�� && If ��i = 1, E(t��i)�� result�� ����
	for(i=0; i<cnt; i++){
		if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pub, prv, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
		{
			result[(*result_num)]=cand[i]; //result�� ����
			(*result_num)++;
		}	
	}


	// user(Bob)���� ��� ������ ���� random �� ����
	for(i=0;i<(*result_num);i++){
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,result[i][j],result[i][j],cipher_rand);
			
		}
	}

	return FsRange_Bob(result, rand, (*result_num), dim);
}

paillier_ciphertext_t*** protocol::sNodeRetrievalforRange(paillier_ciphertext_t*** data, paillier_ciphertext_t*** cipher_qLL_bit, paillier_ciphertext_t*** cipher_qRR_bit, boundary* node, int NumData, int NumNode, int* cnt, int* NumNodeGroup, int type)
{
	printf("\n===== Now sNodeRetrievalforRange starts =====\n");
	int i=0, j=0, m=0;

	std::chrono::system_clock::time_point startTime;
	std::chrono::system_clock::time_point endTime;
	std::chrono::duration<float> duration_sec;		




	int** node_group;

	paillier_ciphertext_t*** cipher_nodeLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** cipher_nodeRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);

	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	for( i = 0 ; i < NumNode ; i++ ){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);
	}

	float progress = 0.1;
	

	printf("\n=== Node SBD & SRO start ===\n");
	startTime = std::chrono::system_clock::now();
	// �� ��� ��Ʈ ��ȯ ���� �� SRO ȣ��
	for( i = 0 ; i < NumNode ; i++ ) {
		startTime = std::chrono::system_clock::now();	
		for( j = 0 ; j < dim ; j++ ) {
			cipher_nodeLL_bit[j] = SBD_for_SRO(node[i].LL[j], 0);				// node LL bound ��ȯ
			cipher_nodeRR_bit[j] = SBD_for_SRO(node[i].RR[j], 1);	 			// node RR bound ��ȯ
		}
		endTime = std::chrono::system_clock::now();
		duration_sec = endTime - startTime;
		time_variable["nodeRetrievalSBD"] = time_variable.find("nodeRetrievalSBD")->second + duration_sec.count();

		startTime = std::chrono::system_clock::now();	
		if(type == 0) {
			alpha[i] = SRO(cipher_qLL_bit, cipher_qRR_bit, cipher_nodeLL_bit, cipher_nodeRR_bit);		
		}
		else if(type == 1) {
			alpha[i] = faster_SRO(cipher_qLL_bit, cipher_qRR_bit, cipher_nodeLL_bit, cipher_nodeRR_bit);		
		}
		endTime = std::chrono::system_clock::now();
		duration_sec = endTime - startTime;
		time_variable["nodeRetrievalSRO"] = time_variable.find("nodeRetrievalSRO")->second + duration_sec.count();
	}

	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);
	
	for( i = 0 ; i < *NumNodeGroup ; i++ ) {
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


	startTime = std::chrono::system_clock::now();
	// ��� �׷� ���� ������ ������ ����, ���� ������ �����ϴ� ��� �� ������ ����
	for(i=0; i<*NumNodeGroup; i++) {
		remained = node_group[i][0];	 // �ش� ��� �׷� ������ ���� ó���� �����Ͱ� �����ִ� ��尡 ����� ������

		for(j=0; j<FanOut; j++) {

			if(remained == 0)	// �ش� ��� �׷� ������ �� �̻� ó���� �����Ͱ� ���ٸ�, ���� ���� �Ѿ
				break;

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
							paillier_mul(pub, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}
					if(node[nodeId].NumData == j+1)		// �ش� ��尡 ������ �����͸� ó���Ѵٸ�, remained�� 1 ���ҽ�Ŵ
						remained--;
				}
				else {		// �ش� ��忡�� �����Ͱ� ������, ���� ��� �׷� �� �ٸ� ��忡�� ���� ó���� �����Ͱ� �ִ� ��츦 �ڵ鸵
					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(cipher_MAX, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
						} else {
							tmp[z] = SM_p1(cipher_MAX, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							paillier_mul(pub, cand[*cnt][z], cand[*cnt][z], tmp[z]);
						}
					}
				}
			}
			(*cnt)++;	// ��� �׷��� ������ �ѹ��� ������, ������ �ϳ��� �ϼ���
		}
	}
	endTime = std::chrono::system_clock::now();
	duration_sec = endTime - startTime;
	time_variable["nodeRetrieval"] = time_variable.find("nodeRetrieval")->second + duration_sec.count();
	return cand;
}

//KHJ 
int ** protocol::sRange_sub(paillier_ciphertext_t** alpha, int node_num, int * set_num){
	int i = 0, j = 0, k = 0, s = 0;
	int cnt = 0;
	int node_set_size = 0;
	int ** sRange_sub_result;
	paillier_plaintext_t** plain_alpha = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*node_num);
	for( i = 0 ; i < node_num ; i++ )
	{
		plain_alpha[i] = paillier_dec(0, pub, prv, alpha[i]);
		if (strcmp( paillier_plaintext_to_str( plain_alpha[i] ), paillier_plaintext_to_str(plain_one)) == 0) {
			cnt++;
		}
	}
	//printf("cnt : %d\n", cnt);
	*set_num = cnt;
	//printf("set_num : %d\n", *set_num);
	
	if(cnt == 0)
		return 0;

	// cnt ��  dec(alpha) == 1 �� ����
	// node_num %cnt >0 �̸� �ϳ��� �迭�� �߰��ؼ� �����ؾ���
	if( (node_num % cnt) > 0 ){
		node_set_size = node_num/cnt + 1;
	}else{
		node_set_size = node_num/cnt;
	}
	//�ʱ� �迭 -1 �� ����
	sRange_sub_result = (int **)malloc(sizeof(int *)*cnt);
	for( i = 0 ; i < cnt ; i++ )
	{
		sRange_sub_result[i] = (int *)malloc(sizeof(int)*(node_set_size+1));
		for( j = 0 ; j < node_set_size+1 ; j++)
		{
			sRange_sub_result[i][j] = -1;
		}
	}
	j = 0;
	for( i = 0 ; i < node_num ; i++ )
	{
		if (strcmp( paillier_plaintext_to_str( plain_alpha[i] ), paillier_plaintext_to_str(plain_one)) == 0) {
			//printf("equal  :  %d \t input : %d\n", i, j);
			sRange_sub_result[j][1] = i;
			sRange_sub_result[j][0] = 1;
			j++;
		}
	}
	j = 0;
	k = 2;
	s = 2;
	for( i = 0 ; i < node_num ; i++ )
	{
		if (strcmp( paillier_plaintext_to_str( plain_alpha[i] ), paillier_plaintext_to_str(plain_one)) != 0) {
			sRange_sub_result[j][k] = i;
			sRange_sub_result[j][0] = s;
			j++;
			if( j == cnt )
			{
				j = 0;
				k++;
				s++;
			}
		}
	}
	
	return sRange_sub_result;
}



int ** protocol::FsRange_Bob(paillier_ciphertext_t*** cipher_result, int rand, int result_num, int col_num){
	int i = 0 , j = 0;
	paillier_plaintext_t* plain;

	int ** kNN = (int **)malloc(sizeof(int*)*result_num);
	for( int i = 0 ; i < result_num ; i++ )
		kNN[i] = (int*)malloc(sizeof(int)*col_num);
	
	for( int i = 0 ; i < result_num; i++ ){
		for( int j = 0 ; j < col_num ; j++ ){
			plain = paillier_dec(0, pubkey, prvkey, cipher_result[i][j]);
			kNN[i][j] = mpz_get_ui(plain->m);
			kNN[i][j] = kNN[i][j] - rand;
		}
	}	

	paillier_freeplaintext(plain);

	for( int i = 0 ; i < result_num ; i++ ) {
		for( int j = 0 ; j < dim ; j++ ) {
			paillier_freeciphertext(cipher_result[i][j]);
		}
		free(cipher_result[i]);
	}
	free(cipher_result);

	return kNN;
}

// added by mc. 0731
paillier_ciphertext_t* protocol::faster_SRO(paillier_ciphertext_t*** qLL, paillier_ciphertext_t*** qRR, paillier_ciphertext_t*** nodeLL, paillier_ciphertext_t*** nodeRR)
{
	//printf("\n===== Now faster_SRO starts ===== \n");
	
	int i, j = 0;
	
	bool func = false;
	paillier_plaintext_t* Rand_value=paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* cipher_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** cipher_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_G	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	paillier_plaintext_t* alpha = paillier_plaintext_from_ui(0);

	paillier_plaintext_t* final_alpha = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	for(i=0; i<size+1; i++){
		cipher_W[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_G[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_H[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_L[i]	= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_M[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_PI[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp[i]	= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp2[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));

		mpz_init(cipher_W[i]->c);
		mpz_init(cipher_G[i]->c);
		mpz_init(cipher_H[i]->c);
		mpz_init(cipher_L[i]->c);
		mpz_init(cipher_M[i]->c);
		mpz_init(cipher_PI[i]->c);
		mpz_init(temp[i]->c);
		mpz_init(temp2[i]->c);
	}

	for(i=0; i<dim; i++) {		
		for(j=0; j<size+1; j++) {
			// compute W
			if(func){	// true :  F : u>v	
				paillier_subtract(pubkey, cipher_W[j], qLL[i][j], SM_p1(qLL[i][j], nodeRR[i][j]));	
			}else{
				paillier_subtract(pubkey, cipher_W[j], nodeRR[i][j], SM_p1(qLL[i][j], nodeRR[i][j]));
			}

			// compute G
			cipher_G[j]=SBXOR(qLL[i][j],nodeRR[i][j]);

			// compute H
			if(j==0){
				paillier_exp(pubkey,cipher_H[j], cipher_zero, Rand_value);
				paillier_mul(pubkey,cipher_H[j], cipher_H[j], cipher_G[j]);
			}else{
				paillier_exp(pubkey,temp[j],cipher_H[j-1],Rand_value);
				paillier_mul(pubkey,cipher_H[j],temp[j],cipher_G[j]);
			}
		}

		for(j=0; j<size+1; j++){
			// compute PI
			paillier_mul(pubkey, cipher_PI[j], cipher_H[j], cipher_minus);	// PI

			// compute L with enhanced privacy (new IDEA)
			mpz_init(temp[i]->c);
			paillier_exp(pubkey,temp[j],cipher_PI[j],Rand_value);	// randomize PI
			paillier_exp(pubkey,temp2[j],cipher_W[j],Rand_value); // randomize W
			paillier_mul(pubkey,cipher_L[j],temp[j], temp2[j]);
		}

		alpha = faster_SRO2(cipher_L, alpha);

		if(func){
			if(strcmp( paillier_plaintext_to_str(alpha), paillier_plaintext_to_str(plain_zero)) == 0)
				alpha=plain_one;
			else
				alpha=plain_zero;
		}
		
		//�� �κп� ���İ� 0�̸� break!!
		if(strcmp( paillier_plaintext_to_str(alpha), paillier_plaintext_to_str(plain_zero)) == 0)
			return cipher_zero;	

		// AND operation (using SM Protocol) 
		if(i == 0)			
			plain_mul(pubkey,final_alpha,plain_one,alpha);//sm -> �׳� ���� ������ !!
		else
			plain_mul(pubkey,final_alpha,final_alpha,alpha);

	}

	// q.RR �� node.LL�� ���ؼ� ����
	for(i=0; i<dim; i++) {
		for(j=0; j<size+1; j++) {
			// compute W
			if(func){	// true :  F : u>v	
				paillier_subtract(pubkey, cipher_W[j], nodeLL[i][j], SM_p1(nodeLL[i][j], qRR[i][j]));	
			}else{
				paillier_subtract(pubkey, cipher_W[j], qRR[i][j], SM_p1(nodeLL[i][j], qRR[i][j]));
			}

			// compute G
			cipher_G[j]=SBXOR(nodeLL[i][j], qRR[i][j]);

			// compute H
			if(j==0){
				paillier_exp(pubkey,cipher_H[j], cipher_zero, Rand_value);
				paillier_mul(pubkey,cipher_H[j], cipher_H[j], cipher_G[j]);
			}else{
				paillier_exp(pubkey,temp[j],cipher_H[j-1],Rand_value);
				paillier_mul(pubkey,cipher_H[j],temp[j],cipher_G[j]);
			}
		}

		for(j=0; j<size+1; j++){
			// compute PI
			paillier_mul(pubkey, cipher_PI[j], cipher_H[j], cipher_minus);	// PI
			mpz_init(temp[i]->c);
			paillier_exp(pubkey,temp[j],cipher_PI[j],Rand_value);	// randomize PI
			paillier_exp(pubkey,temp2[j],cipher_W[j],Rand_value); // randomize W
			paillier_mul(pubkey,cipher_L[j],temp[j], temp2[j]);
		}

		alpha = faster_SRO2(cipher_L, alpha);
		
		if(func){
			//alpha = SBN(alpha);
			if(strcmp( paillier_plaintext_to_str(alpha), paillier_plaintext_to_str(plain_zero)) == 0)
				alpha=plain_one;
			else
				alpha=plain_zero;
		}
		
		//�� �κп� ���İ� 0�̸� break!!
		if(strcmp( paillier_plaintext_to_str(alpha), paillier_plaintext_to_str(plain_zero)) == 0)
			return cipher_zero;		
	
		// AND operation (using SM Protocol) 
		//final_alpha = SM_p1(final_alpha, alpha);
		plain_mul(pubkey,final_alpha,final_alpha,alpha);

	}

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(cipher_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size+1; i++){
		paillier_freeciphertext(cipher_W[i]);	paillier_freeciphertext(cipher_H[i]);	paillier_freeciphertext(cipher_L[i]);	paillier_freeciphertext(cipher_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return paillier_enc(0, pubkey, final_alpha, paillier_get_rand_devurandom);		
}

// added by mc. 0731   //������ �ٲٰ�
paillier_plaintext_t*  protocol::faster_SRO2(paillier_ciphertext_t** cipher_L, paillier_plaintext_t* alpha)
{
	int i;

	// compute alpha
	for(i=0; i<size+1; i++)	{
		if (strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, cipher_L[i])), paillier_plaintext_to_str(plain_zero)) == 0) {
			//alpha = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
			alpha = plain_zero;
			//plain_print("alpha 0 : ", alpha);
			break;
		}
		else {
			//alpha = paillier_enc(0, pubkey, plain_one, paillier_get_rand_devurandom);
			alpha = plain_one;
			//plain_print("alpha 1 : ", alpha);
		}
	}

	//plain_print("alpha : ", alpha);
	
	return alpha;
}


int** protocol::faster_sRange(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num)
{
	printf("\n===== Now sRange starts =====\n");
	//printf("NumNode : %d\n", NumNode);
	int i=0, j=0, m=0;
	int rand = 5;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	int NumNodeGroup = 0;

	paillier_ciphertext_t*** cipher_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** cipher_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** cipher_nodeLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** cipher_nodeRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);

	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	paillier_print("rand : ",cipher_rand);

	// query ��Ʈ ��ȯ ����
	for(i=0; i<dim; i++) {
		//gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, q.LL[i]));

		cipher_qLL_bit[i] = SBD_for_SRO(q.LL[i], 0);			// query LL bound ��ȯ
		//gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, q.RR[i]));

		cipher_qRR_bit[i] = SBD_for_SRO(q.RR[i], 1);		// query RR bound ��ȯ
	}

	/*	
	printf("Query LL bound\n");
	for(i=0; i<dim; i++) {
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, cipher_qLL_bit[i][j]));
		}
		printf("\n");
	}
	printf("Query RR bound\n");
	for(i=0; i<dim; i++) {	
		for(j=0; j<size+1; j++) {
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, cipher_qRR_bit[i][j]));
		}
		printf("\n");
	}
	*/	

	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��

	paillier_ciphertext_t*** cand ;

	cand = sNodeRetrievalforRange(data, cipher_qLL_bit, cipher_qRR_bit, node, NumData, NumNode, &cnt, &NumNodeGroup, 1);
	totalNumOfRetrievedNodes += NumNodeGroup;

	paillier_ciphertext_t*** candLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** candRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);

	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);

	// cand ��Ʈ ��ȯ ���� �� SPEȣ��
	for(i=0; i<cnt; i++) {
		startTime = clock();
		
		for(j=0; j<dim; j++) {
			candLL_bit[j] = SBD_for_SRO(cand[i][j], 0);			// query cand ��ȯ
			candRR_bit[j] = SBD_for_SRO(cand[i][j], 1);		
		}		
		endTime = clock();
		data_SSED_SBD_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		printf("data SSED & SBD time : %f\n", (float)(endTime-startTime)/(CLOCKS_PER_SEC));

		startTime = clock();

		alpha[i]=faster_SRO(candLL_bit,candRR_bit,cipher_qLL_bit, cipher_qRR_bit);

		endTime = clock();
		data_SPE_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		printf("data SPE time : %f\n", (float)(endTime-startTime)/(CLOCKS_PER_SEC));
	}


	paillier_ciphertext_t*** result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumNodeGroup*FanOut); //�Ҵ��� �����ؾ��ϰ�

	//�� Ȯ�� && If ��i = 1, E(t��i)�� result�� ����
	for(i=0; i<cnt; i++){
		if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pub, prv, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
		{
			result[(*result_num)]=cand[i]; //result�� ����
			(*result_num)++;
		}	
	}

	// user(Bob)���� ��� ������ ���� random �� ����
	for(i=0;i<(*result_num);i++){
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,result[i][j],result[i][j],cipher_rand);
		}
	}
	

	return FsRange_Bob(result, rand, (*result_num), dim);
}

paillier_ciphertext_t* protocol::SPE(paillier_ciphertext_t*** qBound, paillier_ciphertext_t*** nodeLL, paillier_ciphertext_t*** nodeRR)
{
	printf("\n===== Now SPE starts ===== \n");
	
	int i, j = 0;
	
	bool func = false;
	paillier_plaintext_t* Rand_value=paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* cipher_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** cipher_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_G	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);
	paillier_ciphertext_t* final_alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(final_alpha->c);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	for(i=0; i<size+1; i++){
		cipher_W[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_G[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_H[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_L[i]	= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_M[i]	 = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		cipher_PI[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp[i]	= (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		temp2[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));

		mpz_init(cipher_W[i]->c);
		mpz_init(cipher_G[i]->c);
		mpz_init(cipher_H[i]->c);
		mpz_init(cipher_L[i]->c);
		mpz_init(cipher_M[i]->c);
		mpz_init(cipher_PI[i]->c);
		mpz_init(temp[i]->c);
		mpz_init(temp2[i]->c);
	}

	

	// q.LL �� node.RR�� ���ؼ� ����
	for(i=0; i<dim; i++) {
		printf("\n%dth dimension\n", i);
		
		for(j=0; j<size+1; j++) {
			// compute W
			if(func){	// true :  F : u>v	
				paillier_subtract(pubkey, cipher_W[j], qBound[i][j], SM_p1(qBound[i][j], nodeRR[i][j]));					
			}else{
				paillier_subtract(pubkey, cipher_W[j], nodeRR[i][j], SM_p1(qBound[i][j], nodeRR[i][j]));				
			}

			// compute G
			cipher_G[j]=SBXOR(qBound[i][j],nodeRR[i][j]);

			// compute H
			if(j==0){
				paillier_exp(pubkey,cipher_H[j], cipher_zero, Rand_value);
				paillier_mul(pubkey,cipher_H[j], cipher_H[j], cipher_G[j]);
			}else{
				paillier_exp(pubkey,temp[j],cipher_H[j-1],Rand_value);
				paillier_mul(pubkey,cipher_H[j],temp[j],cipher_G[j]);
			}

		}

		for(j=0; j<size+1; j++){
			// compute PI
			paillier_mul(pubkey, cipher_PI[j], cipher_H[j], cipher_minus);	// PI

			// compute L with enhanced privacy (new IDEA)
			mpz_init(temp[i]->c);
			paillier_exp(pubkey,temp[j],cipher_PI[j],Rand_value);	// randomize PI
			paillier_exp(pubkey,temp2[j],cipher_W[j],Rand_value); // randomize W
			paillier_mul(pubkey,cipher_L[j],temp[j], temp2[j]);
		}

		alpha = SRO2(cipher_L, alpha);

		if(func){
			alpha = SBN(alpha);
		}
		//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));
	
		// AND operation (using SM Protocol) 
		if(i == 0)
			final_alpha = SM_p1(cipher_one, alpha);	// ���ʿ��� 1�� AND ������ �ؾ���. �ѹ��̶� 0�� ������ 0���� �ٲ�� ��
		else
			final_alpha = SM_p1(final_alpha, alpha);

		//gmp_printf("final alpha (after) %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));
	}


	// q.RR �� node.LL�� ���ؼ� ����
	for(i=0; i<dim; i++) {
		//printf("\n%dth dimension\n", i);
		
		for(j=0; j<size+1; j++) {
			// compute W
			if(func){	// true :  F : u>v	
				paillier_subtract(pubkey, cipher_W[j], nodeLL[i][j], SM_p1(nodeLL[i][j], qBound[i][j]));	
			}else{
				paillier_subtract(pubkey, cipher_W[j], qBound[i][j], SM_p1(nodeLL[i][j], qBound[i][j]));
			}

			// compute G
			cipher_G[j]=SBXOR(nodeLL[i][j], qBound[i][j]);

			// compute H
			if(j==0){
				paillier_exp(pubkey,cipher_H[j], cipher_zero, Rand_value);
				paillier_mul(pubkey,cipher_H[j], cipher_H[j], cipher_G[j]);
			}else{
				paillier_exp(pubkey,temp[j],cipher_H[j-1],Rand_value);
				paillier_mul(pubkey,cipher_H[j],temp[j],cipher_G[j]);
			}
		}

		for(j=0; j<size+1; j++){
			// compute PI
			paillier_mul(pubkey, cipher_PI[j], cipher_H[j], cipher_minus);	// PI

			// compute L with enhanced privacy (new IDEA)
			mpz_init(temp[i]->c);
			paillier_exp(pubkey,temp[j],cipher_PI[j],Rand_value);	// randomize PI
			paillier_exp(pubkey,temp2[j],cipher_W[j],Rand_value); // randomize W
			paillier_mul(pubkey,cipher_L[j],temp[j], temp2[j]);
		}

		alpha = SRO2(cipher_L, alpha);

		if(func){
			alpha = SBN(alpha);
		}
		//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));
	
		// AND operation (using SM Protocol) 
		final_alpha = SM_p1(final_alpha, alpha);

		//gmp_printf("final alpha (after) %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));
	}

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(cipher_Rand_value);		paillier_freeciphertext(tmp);

	paillier_freeciphertext(alpha);
	
	for(i=0; i<size+1; i++){
		paillier_freeciphertext(cipher_W[i]);	paillier_freeciphertext(cipher_H[i]);	paillier_freeciphertext(cipher_L[i]);	paillier_freeciphertext(cipher_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return final_alpha;
		
}