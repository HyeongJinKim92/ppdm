#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <iostream>
#include "protocol.h"

using namespace std;

// ������ �迭 SkNN
int** protocol::SkNN_B(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, int k, int row_number){
	int i=0, s=0, n=0, t=0, j=0;
	int rand = 5;
	n = row_number;
	
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

	for( i = 0 ; i < n ; i++ ){
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


	printf("\n=== Data SSED & SBD start ===\n");

	startTime = std::chrono::system_clock::now(); // startTime check		
	for( i = 0 ; i < n ; i++ ){
		cipher_distance[i] = SSEDm(query, data[i], dim);
		cipher_SBD_distance[i] = SBD(cipher_distance[i]);		
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["SSED&SBD"] = time_variable.find("SSED&SBD")->second + duration_sec.count();


	for( s = 0 ; s < k ; s++ ){
		printf("\n%dth sMINn start \n", s+1);

		startTime = std::chrono::system_clock::now(); // startTime check		
		cipher_Smin = Smin_n(cipher_SBD_distance, n);	// bit�� ǥ���� ��ȣȭ min �Ÿ� ����	
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["SMIN"] = time_variable.find("SMIN")->second + duration_sec.count();
		
		startTime = std::chrono::system_clock::now(); // startTime check		
		// bit�� ǥ���� ��ȣȭ min �Ÿ��� ��ȣȭ ������ ��ȯ
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			cipher_binary = paillier_create_enc(t);
			cipher_binary = SM_p1(cipher_binary, cipher_Smin[size-j]);
			paillier_mul(pubkey, cipher_min, cipher_binary, cipher_min);
			paillier_freeciphertext(cipher_binary);
			//paillier_print("min : ", cipher_min);
		}
		paillier_print("min dist : ", cipher_min);
		
		printf("\n== recalculate query<->data distances ===\n");
		// ���� iteration���� min������ ���õ� �������� �Ÿ��� secure�ϰ� MAX�� ��ȯ�Ǿ��� ������, 
		// ����-������ �� �Ÿ� ����� ��� �ٽ� ����
		if( s != 0 ){
			for( i = 0 ; i < n ; i++ ){
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
		time_variable["Bi to De"] = time_variable.find("Bi to De")->second + duration_sec.count();


		// ����-������ �Ÿ��� min �Ÿ����� ���� ���� (min �������� ��쿡�� 0���� ����� ����)
		startTime = std::chrono::system_clock::now(); // startTime check		
		for( i = 0 ; i < n ; i++ ){
			paillier_subtract(pubkey, temp_dist, cipher_distance[i], cipher_min);
			cipher_mid[i] = SM_p1(temp_dist, cipher_rand);
			//paillier_print("dist - dist : ",cipher_mid[i]);
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["subtract"] = time_variable.find("subtract")->second + duration_sec.count();

		cipher_V = SkNNm_sub(cipher_mid, n);

		startTime = std::chrono::system_clock::now(); // startTime check		
		// min ������ ����
		for( i = 0 ; i < n ; i++ ){
			for( j = 0 ; j < dim; j++ ){
				cipher_V2[i][j] = SM_p1(cipher_V[i], data[i][j]);
				paillier_mul(pubkey, cipher_result[s][j], cipher_V2[i][j], cipher_result[s][j]);
			}
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["kNN"] = time_variable.find("kNN")->second + duration_sec.count();

		// Data SBOR ����
		startTime = std::chrono::system_clock::now(); // startTime check		
		for( i = 0 ; i < n ; i++ ){
			for( j = 0; j < size ; j++ ){
				cipher_SBD_distance[i][j] = SBOR(cipher_V[i], cipher_SBD_distance[i][j]);
			}
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["update"] = time_variable.find("update")->second + duration_sec.count();

		cipher_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
	}
	
	// user(Bob)���� ��� ������ ���� random �� ����
	for(i=0;i<k;i++){
		for(j=0; j<dim; j++){
			paillier_mul(pubkey, cipher_result[i][j], cipher_result[i][j], cipher_rand);
		}
	}

	paillier_freeciphertext(temp_dist);	
	
	return SkNNm_Bob2(cipher_result, rand, k, dim);
}

// Proposed skNN with secure Index
int** protocol::SkNN_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode) {
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��
	int NumNodeGroup = 0;

	bool verify_flag = false;

	paillier_ciphertext_t*** cand;
	paillier_ciphertext_t*** temp_cand ;
	paillier_ciphertext_t* temp_coord1 ;
	paillier_ciphertext_t* temp_coord2 ;
	paillier_ciphertext_t* temp_coord3 ;
	
	paillier_ciphertext_t*** cipher_nodedist_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumNode);
	paillier_ciphertext_t** temp_nodedist_bit = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t*** cipher_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** cipher_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t**** cipher_nodeLL_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t**** cipher_nodeRR_bit = (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t** shortestPoint = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	for(i=0; i<dim; i++){
		shortestPoint[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(shortestPoint[i]->c);
	}

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* cipher_binary;
	paillier_ciphertext_t* cipher_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_dist	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_dist = paillier_create_enc_zero();
	paillier_ciphertext_t* kth_dist = 0;
	paillier_ciphertext_t** kth_dist_bit = 0;

	paillier_ciphertext_t*** cipher_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*3);
	paillier_ciphertext_t* temp_alpha;
	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	for(i=0; i<NumNode; i++){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);

		cipher_nodedist_bit[i] = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*(size+1));
		for( j = 0 ; j < size+1 ; j++ ){ 
			cipher_nodedist_bit[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cipher_nodedist_bit[i][j]->c);
		}

		cipher_nodeLL_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
		cipher_nodeRR_bit[i] =  (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	}

	for(i = 0 ; i < size+1 ; i++ ){ 
		temp_nodedist_bit[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(temp_nodedist_bit[i]->c);
	}

	for( i = 0 ; i < k ; i++ ){
		cipher_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim; j ++ ){
			cipher_result[i][j] = paillier_create_enc_zero();
		}
	}
		

	startTime = std::chrono::system_clock::now(); // startTime check		
	// query ��Ʈ ��ȯ ����
	for(i=0; i<dim; i++) {
		cipher_qLL_bit[i] = SBD_for_SRO(q[i], 0);		// query LL bound ��ȯ
		cipher_qRR_bit[i] = SBD_for_SRO(q[i], 1);		// query RR bound ��ȯ
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["SBD_Q"] = time_variable.find("SBD_Q")->second + duration_sec.count();


	printf("\n=== Node SBD start ===\n");

	startTime = std::chrono::system_clock::now(); // startTime check		
	// �� ��� ��Ʈ ��ȯ ���� 
	for(i=0; i<NumNode; i++) {	
		for(j=0; j<dim; j++) {
			cipher_nodeLL_bit[i][j] = SBD_for_SRO(node[i].LL[j], 0);				// node LL bound ��ȯ
			cipher_nodeRR_bit[i][j] = SBD_for_SRO(node[i].RR[j], 1);	 			// node RR bound ��ȯ
		}
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["nodeRetrievalSBD"] = time_variable.find("nodeRetrievalSBD")->second + duration_sec.count();



	float progress;
	
	while(1) 
	{
		progress = 0.1;
		startTime = std::chrono::system_clock::now(); // startTime check		
		if(!verify_flag)
			printf("\n=== Node SRO start  ===\n");
		else
			printf("\n=== Node expansion start  ===\n");
		
		for(i=0; i<NumNode; i++) 
		{	
			if(i/(float)NumNode >= progress)
			{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			if(!verify_flag)	{	//  ���� �ܰ迡���� �������� ����
				alpha[i] = SRO(cipher_qLL_bit, cipher_qRR_bit, cipher_nodeLL_bit[i], cipher_nodeRR_bit[i]);		
			}
			else {	// ���� �ܰ迡�� k��° �Ÿ����� ������ �����ϴ� ��带 ã�� ����
				// �� ��忡�� ���ǿ��� �ִ� ���� ã��
				for(j=0; j<dim; j++) {	
					psi[0] = SCMP(cipher_qLL_bit[j], cipher_nodeLL_bit[i][j]);		// ������ ��� (1�̸� q�� ��ǥ�� ����� LL bound���� ũ��, 0�̸� q�� ��ǥ�� ����)
					psi[1] = SCMP(cipher_qLL_bit[j], cipher_nodeRR_bit[i][j]);	// ������ ��� (1�̸� q�� ��ǥ�� ����� UR bound���� �۰�, 0�̸� q�� ��ǥ�� ŭ)
					psi[2] = SBXOR(psi[0], psi[1]);	// psi[2]�� 1�̸�, �ش� ���������� ���� ��ǥ�� ����� LL bound �� UR bound ���̿� ������ �ǹ�

					temp_coord1 = SM_p1(psi[2], q[j]);		// psi[2]�� 1�̸�, ���ǿ����� ������ ���� �ִܰŸ� �̹Ƿ�, ������ ��ǥ�� ����� ��
					temp_coord2 = SM_p1(psi[0], node[i].LL[j]);		// psi[0]�� 0�̸�, LL bound�� ��ǥ�� �츲
					temp_coord3 = SM_p1(SBN(psi[0]), node[i].RR[j]);		// psi[0]�� 1�̸�, RR bound�� ��ǥ�� �츲
					paillier_mul(pubkey, temp_coord3, temp_coord2, temp_coord3);
					temp_coord3 = SM_p1(SBN(psi[2]), temp_coord3);

					paillier_mul(pubkey, shortestPoint[j], temp_coord1, temp_coord3);
				}
				
				temp_nodedist_bit = SBD_for_SRO(SSEDm(q, shortestPoint, dim), 0);
				for(j=0; j<size+1; j++) {
					cipher_nodedist_bit[i][j] = SBOR(cipher_nodedist_bit[i][j], temp_nodedist_bit[j]);
				}

				alpha[i] = SCMP(cipher_nodedist_bit[i], kth_dist_bit);	 // �������� �ִ� �Ÿ��� ����ؼ� k��° ����� �Ÿ��� ��. �������� �Ÿ��� �� ������ 1 ����				
			}

			
			if(!verify_flag)	{	//  ���� �ܰ迡���� �������� ����
				// �̹� �˻��� �Ϸ�� ��尡 �� Ž���Ǵ� ���� �����ϱ� ����, bound�� MAX�� ��ȯ
				for(m=0; m<size+1; m++) {
					cipher_nodedist_bit[i][m] = alpha[i];
				}
			}
		}

		if(!verify_flag)	{
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["nodeRetrievalSRO"] = time_variable.find("nodeRetrievalSRO")->second + duration_sec.count();
		}
		else {
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["excheck"] = time_variable.find("excheck")->second + duration_sec.count();
		}


		cnt = 0;

		if(!verify_flag)	{	// ���� �ܰ谡 �ƴ� �ÿ���, cand�� ����
			startTime = std::chrono::system_clock::now(); // startTime check				
			cand = sNodeRetrievalforkNN(data, cipher_qLL_bit, cipher_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
		
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["nodeRetrieval"] = time_variable.find("nodeRetrieval")->second + duration_sec.count();
		}
		else {		// ���� �ܰ��� �ÿ���, ���� ����� cand�� ��ħ
			startTime = std::chrono::system_clock::now(); // startTime check				
			temp_cand	= sNodeRetrievalforkNN(data, cipher_qLL_bit, cipher_qRR_bit, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["exnodeRetrieval"] = time_variable.find("exnodeRetrieval")->second + duration_sec.count();

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
					cand[i][j] = cipher_result[i][j];
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

		paillier_ciphertext_t** cipher_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t*** cipher_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);
		paillier_ciphertext_t** cipher_V=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t*** cipher_V2 = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);

		paillier_ciphertext_t** cipher_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t** cipher_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		
		for( i = 0 ; i < size; i++ ){
			cipher_Smin[i] = cipher_zero;
		}

		for( i = 0 ; i < cnt ; i++ ){
			cipher_distance[i] 	= cipher_zero;
			cipher_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
			cipher_V[i]	= cipher_zero;
			cipher_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		}

		printf("\n=== Data SSED & SBD start ===\n");
		startTime = std::chrono::system_clock::now(); // startTime check				
		for( i = 0 ; i < cnt ; i++ ){
			cipher_distance[i] = SSEDm(q, cand[i], dim);
			cipher_SBD_distance[i] = SBD(cipher_distance[i]);			
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["SSED&SBD"] = time_variable.find("SSED&SBD")->second + duration_sec.count();

		for( s = 0 ; s < k ; s++ ){
			//printf("\n%dth sMINn start \n", s+1);
			startTime = std::chrono::system_clock::now(); // startTime check				
			cipher_Smin = Smin_n(cipher_SBD_distance, cnt);	// bit�� ǥ���� ��ȣȭ min �Ÿ� ����

			if(!verify_flag)	{
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["SMIN"] = time_variable.find("SMIN")->second + duration_sec.count();
			}
			else	{
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["exSMIN"] = time_variable.find("exSMIN")->second + duration_sec.count();
			}
		
			startTime = std::chrono::system_clock::now(); // startTime check				
			// bit�� ǥ���� ��ȣȭ min �Ÿ��� ��ȣȭ ������ ��ȯ
			for( j = size ; j > 0 ; j-- ){
				t = (int)pow(2, j-1);
				cipher_binary = paillier_create_enc(t);
				cipher_binary = SM_p1(cipher_binary, cipher_Smin[size-j]);
				paillier_mul(pubkey, cipher_min, cipher_binary, cipher_min);
				paillier_freeciphertext(cipher_binary);
				//paillier_print("min : ", cipher_min);
			}
			paillier_print("min dist : ", cipher_min);

			printf("\n== recalculate query<->data distances ===\n");
			// ���� iteration���� min������ ���õ� �������� �Ÿ��� secure�ϰ� MAX�� ��ȯ�Ǿ��� ������, 
			// ����-������ �� �Ÿ� ����� ��� �ٽ� ����
			if( s != 0 ){
				for( i = 0 ; i < cnt; i++ ){
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
			time_variable["Bi to De"] = time_variable.find("Bi to De")->second + duration_sec.count();
			
			startTime = std::chrono::system_clock::now(); // startTime check				
			// ����-������ �Ÿ��� min �Ÿ����� ���� ���� (min �������� ��쿡�� 0���� ����� ����)
			for( i = 0 ; i < cnt ; i++ ){
				paillier_subtract(pubkey, temp_dist, cipher_distance[i], cipher_min);
				cipher_mid[i] = SM_p1(temp_dist, cipher_rand);
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["subtract"] = time_variable.find("subtract")->second + duration_sec.count();


			cipher_V = SkNNm_sub(cipher_mid, cnt);
			

			startTime = std::chrono::system_clock::now(); // startTime check				
			// min ������ ����
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0 ; j < dim; j++ ){
					//gmp_printf("cand : %Zd(%d line)\n", paillier_dec(0, pubkey, prvkey, cand[i][j]), __LINE__);
					cipher_V2[i][j] = SM_p1(cipher_V[i], cand[i][j]);
					if(i==0) {
						cipher_result[s][j] = cipher_V2[i][j];
					}
					else {
						paillier_mul(pubkey, cipher_result[s][j], cipher_V2[i][j], cipher_result[s][j]);
					}
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["kNN"] = time_variable.find("kNN")->second + duration_sec.count();

			
			// Data SBOR ����
			startTime = std::chrono::system_clock::now(); // startTime check				
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0; j < size ; j++ ){
					cipher_SBD_distance[i][j] = SBOR(cipher_V[i], cipher_SBD_distance[i][j]);
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["update"] = time_variable.find("update")->second + duration_sec.count();
		
			cipher_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
		}

		if(!verify_flag)	{	
			kth_dist = SSEDm(q, cipher_result[k-1], dim); // k��° ��������� �Ÿ�
			//gmp_printf("%dth dist : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_dist));

			kth_dist_bit = SBD_for_SRO(kth_dist, 1);
		}

		if(!verify_flag)	{
			verify_flag = true;		// ������ �����ϱ� ���� flag�� true�� ����
		}
		else
			break;	// ������ ����
	}
		

	// user(Bob)���� ��� ������ ���� random �� ����
	for(i=0; i<k; i++){
		//printf("%d final result : ", i);
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,cipher_result[i][j],cipher_result[i][j],cipher_rand);
			//gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, cipher_result[i][j]));
		}
		//printf("\n");
	}

	paillier_freeciphertext(temp_dist);	
	
	return SkNNm_Bob2(cipher_result, rand, k, dim);
}

// Proposed skNN with secure Index + SMSn
int** protocol::SkNN_G(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode) {
	printf("\n=== SkNN_G start ===\n");
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;
	
	char Blink[10];
	
	float progress;
	

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
	paillier_ciphertext_t***	cipher_result	= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	
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
		cipher_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j++ )
		{
			cipher_result[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cipher_result[i][j]->c);
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
		startTime = std::chrono::system_clock::now(); // startTime check				

		if(!verify_flag)
		{	
			{
				for( i = 0 ; i < NumNode ; i ++ )
				{
					alpha[i] = GSRO(q, q, node[i].LL, node[i].RR);
				}
			}
		}else{
			for( i  = 0 ; i < NumNode ; i++ )
			{	
				for(j=0; j<dim; j++) 
				{	
					psi[0]			= GSCMP(q[j], node[i].LL[j]);		// ������ ��� (1�̸� q�� ��ǥ�� ����� LL bound���� ũ��, 0�̸� q�� ��ǥ�� ����)
					psi[1]			= GSCMP(q[j], node[i].RR[j]);	// ������ ��� (1�̸� q�� ��ǥ�� ����� UR bound���� �۰�, 0�̸� q�� ��ǥ�� ŭ)
					psi[2]			= SBXOR(psi[0], psi[1]);	// psi[2]�� 1�̸�, �ش� ���������� ���� ��ǥ�� ����� LL bound �� UR bound ����??������ �ǹ�
					temp_coord1		= SM_p1(psi[2], q[j]);		// psi[2]�� 1�̸�, ���ǿ����� ������ ���� �ִܰŸ� �̹Ƿ�, ������ ��ǥ�� ����� ��
					temp_coord2		= SM_p1(psi[0], node[i].LL[j]);		// psi[0]�� 0�̸�, LL bound�� ��ǥ�� �츲
					temp_coord3		= SM_p1(SBN(psi[0]), node[i].RR[j]);		// psi[0]�� 1�̸�, RR bound�� ��ǥ�� �츲
					paillier_mul(pubkey, temp_coord3, temp_coord2, temp_coord3);
					temp_coord3		= SM_p1(SBN(psi[2]), temp_coord3);
					paillier_mul(pubkey, shortestPoint[j], temp_coord1, temp_coord3);
				}
				TMP_NodeDIST		= SSEDm(q, shortestPoint, dim);
				paillier_subtract(pubkey, TMP_alpha, cipher_one, alpha[i]);
				paillier_mul(pubkey, NodeDIST[i], SM_p1(alpha[i], MAX), SM_p1(TMP_alpha, TMP_NodeDIST));				
				alpha[i]			= GSCMP(NodeDIST[i], K_DIST);
			}
		}

		if(!verify_flag)	
		{
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["nodeRetrievalSRO"] = time_variable.find("nodeRetrievalSRO")->second + duration_sec.count();
		}
		else 
		{
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["excheck"] = time_variable.find("excheck")->second + duration_sec.count();
		}
		cnt = 0;

		if(!verify_flag)
		{	// ��??�ܰ谡 �ƴ� �ÿ���, cand�� ����
			startTime = std::chrono::system_clock::now(); // startTime check				
			cand = GSRO_sNodeRetrievalforkNN(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["nodeRetrieval"] = time_variable.find("nodeRetrieval")->second + duration_sec.count();
		}else{
			startTime = std::chrono::system_clock::now(); // startTime check				
			temp_cand = GSRO_sNodeRetrievalforkNN(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["exnodeRetrieval"] = time_variable.find("exnodeRetrieval")->second + duration_sec.count();

			if(cnt == 0)
			{	// ������ ���� �߰� Ž���� �ʿ��� ��尡 ���� ��츦 ó����
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

			// ���� k���� ����� ����
			//cout << "cand : result "<<endl;
			for(i=0; i<k; i++)
			{
				for( j = 0 ; j < dim; j ++ )
				{
					cand[i][j] = cipher_result[i][j];
				}
			}
			// ���� ã�� �ĺ� ����� ����
			for(i=0; i<cnt; i++)
			{
				for( j = 0 ; j < dim; j ++ )
				{
					cand[i+k][j] = temp_cand[i][j];
				}
			}
			cnt = cnt + k;
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


		startTime = std::chrono::system_clock::now(); // startTime check				
		for( i = 0 ; i < cnt ; i ++)
		{
			DIST[i]			= SSEDm(cand[i], q, dim);
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["SSED"] = time_variable.find("SSED")->second + duration_sec.count();


		for(s = 0 ; s < k ; s ++ )
		{
			cout << s+1 <<" th KNN Start !!!!!!!!!!!!"<<endl;

			
			startTime = std::chrono::system_clock::now(); // startTime check				

			idx = SMSn(DIST, cnt);

			if(!verify_flag) {
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["SMIN"] = time_variable.find("SMIN")->second + duration_sec.count();
			}
			else {
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["exSMIN"] = time_variable.find("exSMIN")->second + duration_sec.count();
			}
			

			MIN = DIST[idx];
			gmp_printf("MIN : %Zd \n", paillier_dec(0, pubkey, prvkey, MIN));


			startTime = std::chrono::system_clock::now(); // startTime check				
			for( i = 0 ; i < cnt ; i++ )
			{
				paillier_subtract(pubkey, DIST_MINUS_MIN[i], DIST[i], MIN);
				DIST_MINUS_MIN[i] = SM_p1(DIST_MINUS_MIN[i], C_RAND);		
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["subtract"] = time_variable.find("subtract")->second + duration_sec.count();
		
			V = SkNNm_sub(DIST_MINUS_MIN, cnt);

			startTime = std::chrono::system_clock::now(); // startTime check				
			for( i = 0 ; i < cnt ; i++ )
			{			
				paillier_subtract(pubkey, TMP_alpha, cipher_one, V[i]);
				paillier_mul(pubkey, DIST[i], SM_p1(V[i], MAX), SM_p1(TMP_alpha, DIST[i]));
				for( j = 0 ; j < dim; j++ )
				{
					V2[i][j] = SM_p1(V[i], cand[i][j]);
					if( i == 0 )
					{
						cipher_result[s][j] = V2[i][j];
					}
					else
					{
						paillier_mul(pubkey, cipher_result[s][j], V2[i][j], cipher_result[s][j]);
					}
				}
			}
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["kNN"] = time_variable.find("kNN")->second + duration_sec.count();

		}

		if(!verify_flag)
		{	
			K_DIST = SSEDm(q, cipher_result[k-1], dim); // k��° ��������� �Ÿ�
		}

		if(!verify_flag)
			verify_flag = true;
		else
			break;
	}

	for(i=0; i<k; i++)
	{
		for(j=0; j<dim; j++)
		{
			paillier_mul(pubkey,cipher_result[i][j],cipher_result[i][j],C_RAND);
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

	return SkNNm_Bob2(cipher_result, rand, k, dim);
}

paillier_ciphertext_t*** protocol::sNodeRetrievalforkNN(paillier_ciphertext_t*** data, paillier_ciphertext_t*** cipher_qLL_bit, paillier_ciphertext_t*** cipher_qRR_bit, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now sNodeRetrievalforkNN starts =====\n");
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
			if(remained == 0)	// �ش� ��� �׷� ������ �� �̻� ó���� �����Ͱ� ���ٸ�, ���� ���� �Ѿ
				break;

			if( (i*FanOut + j) / (float)(*NumNodeGroup*FanOut) >= progress)	{
				printf("%.0f%%... ", progress*100);
				progress += 0.1;
			}

			for(m=1; m<=node_group[i][0]; m++) {	 // 0������ �ش� ��� �׷쿡 ��� ��?�?�ִ�?��?����Ǿ� ����
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
				else {		// �ش� ��忡�� �����Ͱ� ������, ���� ��� �׷� �� �ٸ� ��忡�� ���� ó���� �����Ͱ� �ִ� ��츦 �ڵ鸵
					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(cipher_MAX, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
						} else {
							tmp[z] = SM_p1(cipher_MAX, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
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




int ** protocol::GSRO_SkNNb(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode)
{
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 10;
	bool verify_flag = false;

	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��
	int NumNodeGroup = 0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;
	
	int * idx = (int *)malloc(sizeof(int)*k);
	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	paillier_ciphertext_t*** cand;
	paillier_ciphertext_t* cipher_Max = paillier_create_enc(10000);
	paillier_ciphertext_t* temp_coord1 ;
	paillier_ciphertext_t* temp_coord2 ;
	paillier_ciphertext_t* temp_coord3 ;
	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*3);
	paillier_ciphertext_t* temp_alpha = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp_alpha->c);
	paillier_ciphertext_t*** temp_cand ;
	paillier_ciphertext_t** cipher_nodedist = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumNode);

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
		cipher_nodedist[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(cipher_nodedist[i]->c);	
	}

	paillier_ciphertext_t*** cipher_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t)*k);
	for( i = 0 ; i < k ; i++){
		cipher_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j++){
			cipher_result[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cipher_result[i][j]->c);
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

			if(!verify_flag)	{	//  ���� �ܰ迡���� �������� ����
				for( i = 0 ; i < NumNode; i++){
					alpha[i] = GSRO( q, q, node[i].LL, node[i].RR);
				}
			}

			else {	// ���� �ܰ迡�� k��° �Ÿ����� ������ �����ϴ� ��带 ã�� ����
				// �� ��忡�� ���ǿ��� �ִ� ���� ã��
				for(j=0; j<dim; j++) {	
					psi[0] = GSCMP(q[j], node[i].LL[j]);		// ������ ��� (1�̸� q�� ��ǥ�� ����� LL bound���� ũ��, 0�̸� q�� ��ǥ�� ����)
					psi[1] = GSCMP(q[j], node[i].RR[j]);	// ������ ��� (1�̸� q�� ��ǥ�� ����� UR bound���� �۰�, 0�̸� q�� ��ǥ�� ŭ)
					psi[2] = SBXOR(psi[0], psi[1]);	// psi[2]�� 1�̸�, �ش� ���������� ���� ��ǥ�� ����� LL bound �� UR bound ���̿� ������ �ǹ�
					temp_coord1 = SM_p1(psi[2], q[j]);		// psi[2]�� 1�̸�, ���ǿ����� ������ ���� �ִܰŸ� �̹Ƿ�, ������ ��ǥ�� ����� ��
					temp_coord2 = SM_p1(psi[0], node[i].LL[j]);		// psi[0]�� 0�̸�, LL bound�� ��ǥ�� �츲
					temp_coord3 = SM_p1(SBN(psi[0]), node[i].RR[j]);		// psi[0]�� 1�̸�, RR bound�� ��ǥ�� �츲
					paillier_mul(pubkey, temp_coord3, temp_coord2, temp_coord3);
					temp_coord3 = SM_p1(SBN(psi[2]), temp_coord3);
					paillier_mul(pubkey, shortestPoint[j], temp_coord1, temp_coord3);
				}
				temp_nodedist = SSEDm(q, shortestPoint, dim);
				paillier_subtract(pubkey, temp_alpha, cipher_one, alpha[i]);
				paillier_mul(pubkey, cipher_nodedist[i], SM_p1(alpha[i], cipher_Max), SM_p1(temp_alpha, temp_nodedist));
				alpha[i] = GSCMP(cipher_nodedist[i], kth_dist);	 // �������� �ִ� �Ÿ��� ����ؼ� k��° ����� �Ÿ��� ��. �������� �Ÿ��� �� ������ 1 ����
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

		if(!verify_flag)	{	// ���� �ܰ谡 �ƴ� �ÿ���, cand�� ����
			startTime = clock();
			
			cand = GSRO_sNodeRetrievalforkNN(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
		
			endTime = clock();
			data_extract_first_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_first_time);

		}
		else {
			// ���� �ܰ��� �ÿ���, ���� ����� cand�� ��ħ
			startTime = clock();
			cout << "second GSRO_sNodeRetrievalforkNN start" <<endl;
			temp_cand = GSRO_sNodeRetrievalforkNN(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = clock();
			data_extract_second_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_second_time);
			cout << "0" <<endl;
			
			if(cnt == 0) {	// ������ ���� �߰� Ž���� �ʿ��� ��尡 ���� ��츦 ó����
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
			// ���� k���� ����� ����
			for(i=0; i<k; i++) {
				for( j = 0 ; j < dim; j ++ ){
					cand[i][j] = cipher_result[i][j];
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
		paillier_ciphertext_t** cipher_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		for( i = 0 ; i < cnt ; i++ ){
			cipher_distance[i] = SSEDm(q, cand[i], dim);
			gmp_printf("%d  distance : %Zd\n", i , paillier_dec(0, pub, prv, cipher_distance[i]));
		}
		idx = SkNNb_C2(cipher_distance, k, cnt);
		
		for( i = 0 ; i < k ; i++ ){
			printf("idx :  %d\n", idx[i]);
			for( j = 0 ; j < dim ; j++ ){
				cipher_result[i][j] = cand[idx[i]][j];
			}
		}
		
		if(!verify_flag)	{	
			kth_dist = SSEDm(q, cipher_result[k-1], dim); // k��° ��������� �Ÿ�
		}
		
		if(!verify_flag)	{
			verify_flag = true;		// ������ �����ϱ� ���� flag�� true�� ����
		}
		else{
			break;	// ������ ����
		}
	}
	cout<<"break1"<<endl;
	for( i = 0 ; i < k ; i++ ){
		for( j = 0 ; j < dim ; j++ ){
			paillier_mul(pub, cipher_result[i][j], cipher_result[i][j], cipher_rand);
			gmp_printf("%Zd\t", paillier_dec(0, pub, prv, cipher_result[i][j]));
		}
		cout <<endl;
	}

	free(cipher_Max);
	free(idx);
	
	return SkNNb_Bob(cipher_result, rand, k);
}
// added by KHI. 150808
paillier_ciphertext_t* protocol::GSCMP(paillier_ciphertext_t* u, paillier_ciphertext_t* v){
	//cout << "GSCMP Start" << endl;
	//random value setting = 5
	paillier_ciphertext_t* cipher_Rand_value = paillier_create_enc(5);
	//���� ��츦 �Ǻ��� �� �ʿ�.
	paillier_plaintext_t* shift = paillier_plaintext_from_ui(2);
	
	//random value �� arg �� u , v �� ���ؼ� ������ ru , rv�� �����Ҵ�
	paillier_ciphertext_t* ru = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	paillier_ciphertext_t* rv = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(ru->c);
	mpz_init(rv->c);

	//ru, rv, rand_value �� ����� ������� ������ alpha�� ��?? �޸𸮸� �Ҵ��� �ʿ�� ���� �ֳ��ϸ� GSCMP_sub���� alpha��� ������ �޸𸮸� �Ҵ��� �ְ� �װ��� �ּҰ��� �����ϱ� �����̴�.
	paillier_ciphertext_t* alpha;
	
	//u * 2, v * 2 + 1
	
	paillier_exp(pubkey, ru, u, shift);
	
	paillier_exp(pubkey, rv, v, shift);
	paillier_mul(pubkey, rv, rv, cipher_one); 
	
	paillier_mul(pubkey, ru, ru, cipher_Rand_value); 
	paillier_mul(pubkey, rv, rv, cipher_Rand_value); 
	
	alpha = GSCMP_sub( ru, rv, cipher_Rand_value);
	//�޸� ����
	free(cipher_Rand_value);
	free(ru);
	free(rv);
	//gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, alpha));
	
	return alpha;
}
paillier_ciphertext_t* protocol::GSCMP_sub(paillier_ciphertext_t* ru, paillier_ciphertext_t* rv, paillier_ciphertext_t* cipher_Rand){
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
	plain_R	= paillier_dec(0, pub, prv, cipher_Rand);
	
	int a = mpz_get_ui(plain_U->m);
	int b = mpz_get_ui(plain_V->m);
	int R = mpz_get_ui(plain_R->m);

	mtx.lock();
	makeGate();
	mtx.unlock();

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
	//printf("\n===== Now SCMP starts ===== \n");
	
	int i, j = 0;
	
	bool func = false;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* cipher_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** cipher_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_G		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

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
		
	for(i=0; i<size+1; i++) {
		// compute W
		if(func){	// true :  F : u>v	
			paillier_subtract(pubkey, cipher_W[i], u[i], SM_p1(u[i], v[i]));	
		}else{
			paillier_subtract(pubkey, cipher_W[i], v[i], SM_p1(u[i], v[i]));
		}

		// compute G
		cipher_G[i]=SBXOR(u[i], v[i]);

		// compute H
		if(i==0){
			paillier_exp(pubkey,cipher_H[i], cipher_zero, Rand_value);
			paillier_mul(pubkey,cipher_H[i], cipher_H[i], cipher_G[i]);
		}else{
			paillier_exp(pubkey,temp[i],cipher_H[i-1],Rand_value);
			paillier_mul(pubkey,cipher_H[i],temp[i],cipher_G[i]);
		}
	}

	// debugging
	/*
	printf("\nW : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_W[i]));
	}
	printf("\nG : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_G[i]));
	}

	printf("\nH : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_H[i]));
	}
	*/

	for(i=0; i<size+1; i++){
		// compute PI
		paillier_mul(pubkey, cipher_PI[i], cipher_H[i], cipher_minus);	// PI

		// compute L with enhanced privacy (new IDEA)
		mpz_init(temp[i]->c);
		paillier_exp(pubkey,temp[i],cipher_PI[i],Rand_value);	// randomize PI
		paillier_exp(pubkey,temp2[i],cipher_W[i],Rand_value); // randomize W
		paillier_mul(pubkey,cipher_L[i],temp[i], temp2[i]);
	}
	
	// debugging
	/*
	printf("\nPI : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_PI[i]));
	}

	printf("\nL : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_L[i]));
	}
	printf("\n");
	*/

	alpha = SRO2(cipher_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}
	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(cipher_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size+1; i++){
		paillier_freeciphertext(cipher_W[i]);	paillier_freeciphertext(cipher_H[i]);	paillier_freeciphertext(cipher_L[i]);	paillier_freeciphertext(cipher_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return alpha;
	
}

paillier_ciphertext_t* protocol::SCMP_M(paillier_ciphertext_t** u, paillier_ciphertext_t** v)
{
	int i, j = 0;
	
	bool func = false;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* cipher_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** cipher_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_G		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));

	for(i=0; i<size; i++){
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

	for(i=0; i<size; i++) {
		// compute W
		if(func){	// true :  F : u>v	
			paillier_subtract(pubkey, cipher_W[i], u[i], SM_p1(u[i], v[i]));	
		}else{
			paillier_subtract(pubkey, cipher_W[i], v[i], SM_p1(u[i], v[i]));
		}
		// compute G
		cipher_G[i]=SBXOR(u[i], v[i]);
		// compute H
		if(i==0){
			paillier_exp(pubkey,cipher_H[i], cipher_zero, Rand_value);
			paillier_mul(pubkey,cipher_H[i], cipher_H[i], cipher_G[i]);
		}else{
			paillier_exp(pubkey,temp[i],cipher_H[i-1],Rand_value);
			paillier_mul(pubkey,cipher_H[i],temp[i],cipher_G[i]);
		}
	}


	for(i=0; i<size; i++){
		// compute PI
		paillier_mul(pubkey, cipher_PI[i], cipher_H[i], cipher_minus);	// PI

		// compute L with enhanced privacy (new IDEA)
		mpz_init(temp[i]->c);
		paillier_exp(pubkey,temp[i],cipher_PI[i],Rand_value);	// randomize PI
		paillier_exp(pubkey,temp2[i],cipher_W[i],Rand_value); // randomize W
		paillier_mul(pubkey,cipher_L[i],temp[i], temp2[i]);
	}
	

	alpha = SRO2(cipher_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}
	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(cipher_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size; i++){
		paillier_freeciphertext(cipher_W[i]);	paillier_freeciphertext(cipher_H[i]);	paillier_freeciphertext(cipher_L[i]);	paillier_freeciphertext(cipher_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return alpha;
	
}
paillier_ciphertext_t* protocol::SCMP_Clustering(paillier_ciphertext_t** u, paillier_ciphertext_t** v)
{
	int i, j = 0;
	
	bool func = true;
	paillier_plaintext_t* Rand_value = paillier_plaintext_from_ui(5);
	paillier_ciphertext_t* cipher_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** cipher_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_G		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** cipher_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

	paillier_ciphertext_t* tmp		= paillier_create_enc_zero();
	paillier_ciphertext_t** temp	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
	paillier_ciphertext_t** temp2	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));

	for(i=0; i<size; i++){
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

	for(i=0; i<size; i++) {
		// compute W
		if(!func){	// true :  F : u>v	
			paillier_subtract(pubkey, cipher_W[i], u[i], SM_p1(u[i], v[i]));	
		}else{
			paillier_subtract(pubkey, cipher_W[i], v[i], SM_p1(u[i], v[i]));
		}
		// compute G
		cipher_G[i]=SBXOR(u[i], v[i]);
		// compute H
		if(i==0){
			paillier_exp(pubkey,cipher_H[i], cipher_zero, Rand_value);
			paillier_mul(pubkey,cipher_H[i], cipher_H[i], cipher_G[i]);
		}else{
			paillier_exp(pubkey,temp[i],cipher_H[i-1],Rand_value);
			paillier_mul(pubkey,cipher_H[i],temp[i],cipher_G[i]);
		}
	}


	for(i=0; i<size; i++){
		// compute PI
		paillier_mul(pubkey, cipher_PI[i], cipher_H[i], cipher_minus);	// PI

		// compute L with enhanced privacy (new IDEA)
		mpz_init(temp[i]->c);
		paillier_exp(pubkey,temp[i],cipher_PI[i],Rand_value);	// randomize PI
		paillier_exp(pubkey,temp2[i],cipher_W[i],Rand_value); // randomize W
		paillier_mul(pubkey,cipher_L[i],temp[i], temp2[i]);
	}
	
	alpha = SRO2(cipher_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}
	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(cipher_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size; i++){
		paillier_freeciphertext(cipher_W[i]);	paillier_freeciphertext(cipher_H[i]);	paillier_freeciphertext(cipher_L[i]);	paillier_freeciphertext(cipher_M[i]);
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
	paillier_ciphertext_t* cipher_Rand_value = paillier_enc(0, pubkey, Rand_value, paillier_get_rand_devurandom);	

	paillier_ciphertext_t** cipher_W		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_G		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_H		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_L		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_M		 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t** cipher_PI	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));

	paillier_ciphertext_t* alpha = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(alpha->c);

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
		
	for(i=0; i<size; i++) {
		// compute W
		if(func){	// true :  F : u>v	
			paillier_subtract(pubkey, cipher_W[i], u[i], SM_p1(u[i], v[i]));	
		}else{
			paillier_subtract(pubkey, cipher_W[i], v[i], SM_p1(u[i], v[i]));
		}

		// compute G
		cipher_G[i]=SBXOR(u[i], v[i]);

		// compute H
		if(i==0){
			paillier_exp(pubkey,cipher_H[i], cipher_zero, Rand_value);
			paillier_mul(pubkey,cipher_H[i], cipher_H[i], cipher_G[i]);
		}else{
			paillier_exp(pubkey,temp[i],cipher_H[i-1],Rand_value);
			paillier_mul(pubkey,cipher_H[i],temp[i],cipher_G[i]);
		}
	}

	// debugging
	/*
	printf("\nW : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_W[i]));
	}
	printf("\nG : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_G[i]));
	}

	printf("\nH : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_H[i]));
	}
	*/

	for(i=0; i<size; i++){
		// compute PI
		paillier_mul(pubkey, cipher_PI[i], cipher_H[i], cipher_minus);	// PI

		// compute L with enhanced privacy (new IDEA)
		mpz_init(temp[i]->c);
		paillier_exp(pubkey,temp[i],cipher_PI[i],Rand_value);	// randomize PI
		paillier_exp(pubkey,temp2[i],cipher_W[i],Rand_value); // randomize W
		paillier_mul(pubkey,cipher_L[i],temp[i], temp2[i]);
	}
	
	// debugging
	/*
	printf("\nPI : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_PI[i]));
	}

	printf("\nL : ");
	for(i=0; i<size+1; i++) {
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cipher_L[i]));
	}
	printf("\n");
	*/

	alpha = SRO2(cipher_L, alpha);

	if(func){
		alpha = SBN(alpha);
	}
	//gmp_printf("alpha %Zd\n", paillier_dec(0, pubkey, prvkey, alpha));

	// free memory
	paillier_freeplaintext(Rand_value);	 
	paillier_freeciphertext(cipher_Rand_value);		paillier_freeciphertext(tmp);

	for(i=0; i<size+1; i++){
		paillier_freeciphertext(cipher_W[i]);	paillier_freeciphertext(cipher_H[i]);	paillier_freeciphertext(cipher_L[i]);	paillier_freeciphertext(cipher_M[i]);
		paillier_freeciphertext(temp[i]); paillier_freeciphertext(temp2[i]);	
	}

	return alpha;
	
}

int protocol::SMSn(paillier_ciphertext_t** cipher, int cnt){
	int i = 0, s = 0;
	int buff = 0;
	int R1 = 5;
	int quot = 0;
	int remain = 0;
	int iter = 0;
	int bundle = 0;
	int cipher_idx = 0;
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
			cout << i+1 << " cipher : ";
			gmp_printf(" %Zd\n",paillier_dec(0, pub, prv, cipher[i]));
		}
	}

	paillier_plaintext_t* shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* tmp= paillier_plaintext_from_ui(0);
	paillier_plaintext_t* p1_array = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* R1_array = paillier_plaintext_from_ui(R1);
	paillier_ciphertext_t* c1_array ;
	paillier_ciphertext_t* cR1_array;
	paillier_ciphertext_t* cipher_tmp = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));	
	mpz_init(cipher_tmp->c);

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
			paillier_exp(pub, cipher_tmp, cipher[cipher_idx++], shift);
			paillier_mul(pubkey, c1_array, c1_array, cipher_tmp); //E_add
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
		mpz_init(cipher_tmp->c);
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
			cout << i+1 << " cipher : "<< R_array_result[i] + V_array_result[i] <<endl;
		}
	}
	
	if(!smsn_flag){
		smsn_flag = true;
		makeGate();  // �Ź��ؾ��ϴ��� ? üũ
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
	free(cipher_tmp);
	
	if(Print)
	{
	
		cout << "Min Value : " << smsn_tmp_val+smsn_tmp_rand <<"   idx : "<<smsn_result_idx;
		gmp_printf("  %Zd\n",paillier_dec(0, pub, prv, cipher[smsn_result_idx]));	
	}
	
	return smsn_result_idx;
}

paillier_ciphertext_t* protocol::SMSn_sub(paillier_ciphertext_t* cipher, paillier_ciphertext_t* rand ,int cnt){
	int i = 0;
	int * R1 = (int *)malloc(sizeof(int)*cnt);
	for( i = 0 ; i < cnt ; i ++){
		R1[i] = 4;
	}

	paillier_plaintext_t* shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* tmp = paillier_plaintext_from_ui(0);
	paillier_plaintext_t* rand_tmp = paillier_plaintext_from_ui(0);

	paillier_plaintext_t** R1_array = (paillier_plaintext_t**)malloc(sizeof(paillier_plaintext_t*)*cnt);

	paillier_plaintext_t* value_array = paillier_dec(0, pub, prv, cipher);  //A���� ���޵Ǵ� packing ��
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

int ** protocol::SkNN_plain(paillier_ciphertext_t*** cipher, paillier_ciphertext_t** query, int k, int row_number){
	int i = 0, j = 0;
	int rand = 10;
	int * idx = (int *)malloc(sizeof(int)*k);
	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	paillier_ciphertext_t** cipher_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*row_number);
	paillier_ciphertext_t*** cipher_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t)*k);

	paillier_ciphertext_t** cipher_query = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	
	for( i = 0 ; i < dim ; i++ ){
		cipher_query[i] = query[i];
	}
	for( i = 0 ; i < k ; i++){
		cipher_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	}
	
	
	for( i = 0 ; i < row_number ; i++ ){
		cipher_distance[i] = SSEDm(cipher[i] ,cipher_query, dim);
		//paillier_print("cipher_dist : ", cipher_distance[i]);
	}
	idx = SkNNb_C2(cipher_distance, k, row_number);
	
	for( i = 0 ; i < k ; i++ ){
		//printf("idx :  %d\n", idx[i]);
		for( j = 0 ; j < dim ; j++ ){
			cipher_result[i][j] = cipher[idx[i]][j];
			paillier_mul(pub, cipher_result[i][j], cipher_result[i][j], cipher_rand);
			//paillier_print(" ",cipher_result[i][j]);
		}
		//printf("\n");
	}

	free(idx);
	free(cipher_distance);
	free(cipher_query);
	
	return SkNNb_Bob(cipher_result, rand, k);
}
int ** protocol::SkNNb_Bob(paillier_ciphertext_t*** cipher_result, int rand, int k){
	int i = 0, j = 0;
	int ** result = (int **)malloc(sizeof(int*)*k);
	paillier_plaintext_t* plain;
	for( i = 0 ; i < k ; i ++){
		result[i] = (int *)malloc(sizeof(int)*dim);
		for( j = 0 ; j < dim ; j ++){
			plain = paillier_dec(0, pub, prv, cipher_result[i][j]);
			result[i][j] = mpz_get_si(plain->m);
			result[i][j] = result[i][j] - rand;
			//printf("%d ", result[i][j]);
			free(cipher_result[i][j]);
		}
		//printf("\n");
	}

	return result;
}
int * protocol::SkNNb_C2(paillier_ciphertext_t** cipher_dist, int k, int row_number){
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
		plain[i] = paillier_dec(0, pub, prv, cipher_dist[i]);
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

	for( i = 0 ; i < k ; i++ ){
		result[i] = idx[i];
		//printf("result :  %d\n", result[i]);
	}
	free(plain);
	free(idx);
	free(toint);
	return result;
}



int * protocol::SkNN_B(paillier_ciphertext_t** cipher, int Q, int k, int row_number){
	int i=0, s=0, n=0, t=0, j=0;
	int rand = 5;
	n = row_number;
	
	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* cipher_binary;
	paillier_ciphertext_t* cipher_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_dist	= paillier_create_enc_zero();
	
	printf("%d\n", Q);
	printf("%d\n", k);
	printf("%d\n", row_number);
	
	paillier_ciphertext_t* cipher_query = paillier_create_enc(Q);

	paillier_print("query : ",cipher_query);
	for( i = 0 ; i < n ; i++){
		paillier_print("i : ",cipher[i]);
	}

	paillier_ciphertext_t** cipher_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** cipher_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** cipher_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
	paillier_ciphertext_t** cipher_V;
	paillier_ciphertext_t** cipher_V2 = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	
	paillier_ciphertext_t*  cipher_result = paillier_create_enc_zero();
	paillier_ciphertext_t** cipher_result_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	
	paillier_ciphertext_t*** cipher_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*n);
	for( i = 0 ; i < n ; i++ ){
		cipher_mid[i] = paillier_create_enc_zero();	// �ʱ�ȭ�� ���൵ �ɵ�!!??
		cipher_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		for( j = 0 ; j < size ; j ++ ){
			cipher_SBD_distance[i][j] = paillier_create_enc_zero();	// �ʱ�ȭ�� ���൵ �ɵ�!!??
		}
	}

	for( i = 0 ; i < size; i++ ){
		cipher_Smin[i] = cipher_zero;
	}
	
	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	
	paillier_print("rand : ",cipher_rand);
	
	for( i = 0 ; i < n ; i++ ){
		cipher_distance[i] = SSED(cipher_query, cipher[i]);
		paillier_print("dist : ", cipher_distance[i]);	
		cipher_SBD_distance[i] = SBD(cipher_distance[i]);
		for( j = 0 ; j < size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, cipher_SBD_distance[i][j]));
		}
		printf("\n");
	}
	for( s = 0 ; s < k ; s++ ){
		
		cipher_Smin = Smin_n(cipher_SBD_distance, n);

		for( j = 0 ; j < size ; j++ ){
			gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, cipher_Smin[j]));
		}
		printf("\n");	
		
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			cipher_binary = paillier_create_enc(t);
			cipher_binary = SM_p1(cipher_binary, cipher_Smin[size-j]);
			paillier_mul(pubkey, cipher_min, cipher_binary, cipher_min);

			
		}
		paillier_print("min : ",cipher_min);
		
		if( s != 0 ){
			for( i = 0 ; i < n ; i++ ){
				for( j = size ; j > 0 ; j--){
					t = (int)pow(2, j-1);
					cipher_binary = paillier_create_enc(t);
					cipher_binary = SM_p1(cipher_binary, cipher_SBD_distance[i][size-j]);
					paillier_mul(pubkey, cipher_dist, cipher_binary, cipher_dist);
				}
				cipher_distance[i] = cipher_dist;
				cipher_dist = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	 // ���� cipher_zero ����ص� ���� ������!!??
			}
		}

	
		for( i = 0 ; i < n ; i++ ){
			paillier_print("distance : ",cipher_distance[i]);
		}
		
		for( i = 0 ; i < n ; i++ ){
			paillier_ciphertext_t* cipher_temp;
			paillier_subtract(pubkey, cipher_temp, cipher_distance[i], cipher_min);
			cipher_mid[i] = SM_p1(cipher_temp, cipher_rand);
			paillier_print("dist - dist : ",cipher_temp);
		}
		cipher_V = SkNNm_sub(cipher_mid, n);
		for( i = 0 ; i < n ; i++ ){
			paillier_print("cipher_V : ",cipher_V[i]);
		}
		

		for( i = 0 ; i < n ; i++ ){
			cipher_V2[i] = SM_p1(cipher_V[i], cipher[i]);
			paillier_mul(pubkey, cipher_result, cipher_V2[i], cipher_result);
		}
		paillier_print("cipher_result : ",cipher_result);
		
		cipher_result_array[s] = cipher_result;
		for( i = 0 ; i < n ; i++ ){
			for( j = 0; j < size ; j++ ){
				cipher_SBD_distance[i][j] = SBOR(cipher_V[i], cipher_SBD_distance[i][j]);
			}
		}
		cipher_result = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);
		cipher_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);
	}

	for( i = 0 ; i < n ; i++){
		for( j = 0 ; j < size ; j++ ){
				gmp_printf("%Zd", paillier_dec(0, pubkey, prvkey, cipher_SBD_distance[i][j]));
			}
		printf("\n");
	}

	for( i = 0 ; i < k ; i++ ){
		paillier_mul(pubkey, cipher_result_array[i],cipher_result_array[i], cipher_rand);
		paillier_print("cipher_result_array : ",cipher_result_array[i]);
	}
	free(cipher_binary);
	free(cipher_min);
	free(pt);
	free(cipher_dist);
	free(cipher_mid);
	free(cipher_Smin);
	free(cipher_V);
	free(cipher_V2);
	free(cipher_result);
	free(cipher_SBD_distance);

	return SkNNm_Bob(cipher_result_array, rand, k);
}

// SkNN subroutine
paillier_ciphertext_t** protocol::SkNNm_sub(paillier_ciphertext_t** cipher_n, int n){
	int i = 0;
	int cnt = 0;
	paillier_plaintext_t* plain = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	paillier_ciphertext_t** cipher_U = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	
	for( i = 0 ; i < n ; i++ ){
		plain = paillier_dec(0, pubkey, prvkey, cipher_n[i]);
		if(  (mpz_cmp(plain_zero->m, plain->m) == 0) && cnt == 0 ){
			cipher_U[i] = paillier_create_enc(1);
			cnt++;
		}else if ((mpz_cmp(plain_zero->m, plain->m) == 0))
		{
			cipher_U[i] = paillier_create_enc(0);
			cnt++;
		}else if((mpz_cmp(plain_zero->m, plain->m) != 0)){
			cipher_U[i] = paillier_create_enc(0);
		}
	}
	if(cnt == n){
		for( i = 0 ; i < n ; i++ ){
			cipher_U[i] = paillier_create_enc(0);
		}
	}
	
	paillier_freeplaintext(plain);

	return cipher_U;
}

// ������ SkNN ��� return �Լ�
int * protocol::SkNNm_Bob(paillier_ciphertext_t** cipher_result_array, int rand, int k){
	int * kNN = (int *)malloc(sizeof(int)* k);
	int i = 0;
	paillier_plaintext_t* plain;
	for( i = 0 ; i < k ; i++ ){
		plain = paillier_dec(0, pubkey, prvkey, cipher_result_array[i]);
		kNN[i] = mpz_get_ui(plain->m);
		kNN[i] = kNN[i] - rand;
	}
	free(plain);
	free(cipher_result_array);
	for( i = 0 ; i< k ; i++ ){
		printf("result : %d\n", kNN[i]);
	}
	return kNN;
}

// ������ SkNN ��� return �Լ�
int ** protocol::SkNNm_Bob2(paillier_ciphertext_t*** cipher_result, int rand, int k, int col_num){
	int i = 0 , j = 0;
	paillier_plaintext_t* plain;

	int ** kNN = (int **)malloc(sizeof(int*)*k);
	for( i = 0 ; i < k ; i++ )
		kNN[i] = (int*)malloc(sizeof(int)*col_num);
	
	for( i = 0 ; i < k ; i++ ){
		for( j = 0 ; j < col_num ; j++ ){
			//paillier_subtract(pubkey, cipher_result[i][j], cipher_result[i][j], cipher_rand);
			plain = paillier_dec(0, pubkey, prvkey, cipher_result[i][j]);
			kNN[i][j] = mpz_get_ui(plain->m);
			kNN[i][j] = kNN[i][j] - rand ; 
		}
	}

	paillier_freeplaintext(plain);
	for(i=0; i<k; i++) {
		for(j=0; j<dim; j++) {
			paillier_freeciphertext(cipher_result[i][j]);
		}
		free(cipher_result[i]);
	}
	free(cipher_result);

	return kNN;
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
	printf("set_num : %d\n", *NumNodeGroup);

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
	int remained = 0;	  // ��� �׷� ������ ���� ó���� �����Ͱ� �����ִ� ��尡 ���??������

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
					dataId = node[nodeId].indata_id[j];	// �ش� ��忡 ����� ��?��?ID�� �ϳ��� ����
					//printf("selected data ID : %d\n", dataId);

					for(z=0; z<dim; z++) {
						if(m == 1) {
							cand[*cnt][z] = SM_p1(data[dataId][z], alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
						} else {
							tmp[z] = SM_p1(data[dataId][z], alpha[nodeId]);  // �ش� ��???ID�� ���� �����Ϳ� ����
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
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

