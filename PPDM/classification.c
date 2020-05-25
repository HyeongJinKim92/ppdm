#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <iostream>
#include "protocol.h"

using namespace std;

paillier_ciphertext_t*** protocol::sNodeRetrievalforClassification(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now sNodeRetrievalforClassification starts =====\n");
	//printf("NumNode : %d\n", NumNode);
	int i=0, j=0, m=0;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	int** node_group;

	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);
	if(*NumNodeGroup == 0)
		return 0;

	//printf("\nset_num : %d\n", *NumNodeGroup);
/*
	for(i=0; i<*NumNodeGroup; i++) {
		printf("%dth Node Group : ", i+1);
		for(j=1; j<=node_group[i][0]; j++) {	// 0������ �ش� ��� �׷쿡 ��� ��尡 �ִ����� ����Ǿ� ����
			printf("%d ", node_group[i][j]);
		}
		printf("\n");
	}
*/
	int nodeId = 0;
	int dataId = 0;
	int z = 0;
	int remained = 0;	  // ��� �׷� ������ ���� ó���� �����Ͱ� �����ִ� ��尡 ����� ������
	
	paillier_ciphertext_t** tmp = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
	for( i = 0 ; i < (dim+1); i++ ){
		tmp[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(tmp[i]->c);
	}

	paillier_ciphertext_t*** cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(*NumNodeGroup)*FanOut);
		
	for( i = 0 ; i < *NumNodeGroup*FanOut ; i++ ){
		cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim+1);
		for( j = 0 ; j < dim+1 ; j ++ ){
			cand[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cand[i][j]->c);
		}
	}

	float progress = 0.1;
	
	// ��� �׷� ���� ������ ������ ����, ���� ������ �����ϴ� ��� �� ������ ����
	//printf("extract data\n");
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
				//printf("selected node ID : %d\n", nodeId);

				if(node[nodeId].NumData >= j+1) 
				{
					dataId = node[nodeId].indata_id[j];	// �ش� ��忡 ����� ������ ID�� �ϳ��� ����
					//printf("selected data ID : %d\n", dataId);

					for(z=0; z<dim+1; z++) {
						if(m == 1) {
							//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, pubkey, prvkey,data[dataId][z]), paillier_dec(0, pubkey, prvkey, alpha[nodeId]));
							//printf("cnt : %d,  dim : %d\n", *cnt, z);
							cand[*cnt][z] = SM_p1(data[dataId][z], alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
						} else {
							//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, pubkey, prvkey,data[dataId][z]), paillier_dec(0, pubkey, prvkey, alpha[nodeId]));
							//printf("cnt : %d,  dim : %d\n", *cnt, z);
							tmp[z] = SM_p1(data[dataId][z], alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, tmp[z]));
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
							//gmp_printf("cnt : %d, (%Zd)\n", *cnt, paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
						}
					}
					//printf("\n");

					if(node[nodeId].NumData == j+1)		// �ش� ��尡 ������ �����͸� ó���Ѵٸ�, remained�� 1 ���ҽ�Ŵ
						remained--;
				}
				else {		// �ش� ��忡�� �����Ͱ� ������, ���� ��� �׷� �� �ٸ� ��忡�� ���� ó���� �����Ͱ� �ִ� ��츦 �ڵ鸵
					for(z=0; z<dim+1; z++) {
						if(m == 1) {
							//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, pubkey, prvkey, ciper_MAX), paillier_dec(0, pubkey, prvkey, alpha[nodeId]));
							//printf("cnt : %d,  dim : %d\n", *cnt, z);
							cand[*cnt][z] = SM_p1(ciper_MAX, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
						} else {
							//gmp_printf("data : %Zd, alpha : %Zd \n", paillier_dec(0, pubkey, prvkey, ciper_MAX), paillier_dec(0, pubkey, prvkey, alpha[nodeId]));
							//printf("cnt : %d,  dim : %d\n", *cnt, z);
							tmp[z] = SM_p1(ciper_MAX, alpha[nodeId]);  // �ش� ������ ID�� ���� �����Ϳ� ����
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
							//gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, tmp[z]));
							paillier_mul(pubkey, cand[*cnt][z], cand[*cnt][z], tmp[z]);
							//gmp_printf("(%Zd)\n", paillier_dec(0, pubkey, prvkey, cand[*cnt][z]));
						}
					}
					//printf("\n");
				}
			}
			(*cnt)++;	// ��� �׷��� ������ �ѹ��� ������, ������ �ϳ��� �ϼ���
		}
	}
	//printf("\n");
	/*
	for(i=0; i<*cnt; i++) {
		gmp_printf("%dth data -> coord : ", i);
		for(j=0; j<dim+1; j++) {
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, cand[i][j]));
		}
		printf("\n");
	}
	*/
	return cand;
}

int** protocol::Classification_G(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, paillier_ciphertext_t** Entire_set, boundary* node, int k, int NumData, int NumNode, int Entire_num) {
	printf("\n=== Classification_G start ===\n");
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;
	
	
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
	paillier_ciphertext_t**		Rabel_Result	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	paillier_ciphertext_t***	Result			= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	
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
		Rabel_Result[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));

		Result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
		for( j = 0 ; j < dim+1 ; j++ )
		{
			Result[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(Result[i][j]->c);
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
				TMP_NodeDIST		= DP_SSED(q, shortestPoint, dim);
				paillier_subtract(pubkey, TMP_alpha, ciper_one, alpha[i]);
				paillier_mul(pubkey, NodeDIST[i], SM_p1(alpha[i], MAX), SM_p1(TMP_alpha, TMP_NodeDIST));				
				alpha[i]			= GSCMP(NodeDIST[i], K_DIST);
				if(Print)
				{
					cout<< i+1 <<"Node Expansion alpha : ";
					gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
				}
			}
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
		{	
			startTime = clock();
			cand = sNodeRetrievalforClassification(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_first_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_first_time);
		}else{
			startTime = clock();
			cout << "second GSRO_sNodeRetrievalforkNN start" <<endl;
			temp_cand = sNodeRetrievalforClassification(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			endTime = clock();
			data_extract_second_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_second_time);

			if(cnt == 0)
			{	// ������ ���� �߰� Ž���� �ʿ��� ��尡 ���� ��츦 ó����
				break;
			}

			cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ )
			{
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim+1);
				for( j = 0 ; j < dim+1 ; j ++ )
				{
					cand[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
					mpz_init(cand[i][j]->c);
				}
			}

			// ���� k���� ����� ����
			//cout << "cand : result "<<endl;
			for( i = 0 ; i < k ; i++ )
			{
				for( j = 0 ; j < dim+1 ; j ++ )
				{
					cand[i][j] = Result[i][j];
				}
			}
			// ���� ã�� �ĺ� ����� ����
			for( i = 0 ; i < cnt ; i++ )
			{
				for( j = 0 ; j < dim+1 ; j ++ )
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
			for( i = 0 ; i < cnt ; i ++)
			{
				for( j = 0 ; j < dim+1 ; j ++ )
				{
					gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, cand[i][j]));
				}
				cout<<endl;
			}
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
			V2[i]				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
			
			mpz_init(DIST_MINUS_MIN[i]->c);
			mpz_init(DIST[i]->c);
			mpz_init(V[i]->c);
			for( j = 0 ; j < dim+1 ; j++ ){
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
			cout << s+1 <<" th Classification Start !!!!!!!!!!!!"<<endl;
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
				for( j = 0 ; j < dim+1 ; j++ )
				{
					V2[i][j] = SM_p1(V[i], cand[i][j]);
					if( i == 0 )
					{
						Result[s][j] = V2[i][j];
					}
					else
					{
						paillier_mul(pubkey, Result[s][j], V2[i][j], Result[s][j]);
					}
				}
			}
		}
		if(Print)
		{
			for( i = 0 ; i < k ; i++ )
			{
				cout << "MIDDLE RESULT : ";
				for( j = 0 ; j < dim+1 ; j++ )
				{
					gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, Result[i][j]));
				}
				cout <<endl;
			}
		}
		if(!verify_flag)
		{	
			K_DIST = DP_SSED(q, Result[k-1], dim); // k��° ��������� �Ÿ�
		}

		if(!verify_flag)
			verify_flag = true;
		else
			break;
	}
	protocol_3_printf(Result, k, dim+1, "min_record", true);
	for( i = 0 ; i < k ; i++ ){
		*Rabel_Result[i] = *Result[i][dim];
	}
	
	paillier_plaintext_t* A = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	int** sknn = (int**)malloc(sizeof(int*)*k);
	for( i = 0 ; i < k ; i++ ){
		sknn[i] = (int*)malloc(sizeof(int)*(dim+1));
		for( j = 0 ; j < dim+1 ; j++ ){
			A = paillier_dec(0, pubkey, prvkey, Result[i][j]);
			sknn[i][j] = mpz_get_ui(A->m);
			printf("%d ", sknn[i][j]);
		}
		printf("\n");
	}
	SCMC(Entire_set, Rabel_Result, Entire_num, k);
	return sknn;
}

int** protocol::Classification_I(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, paillier_ciphertext_t** Entire_set, boundary* node, int k, int NumData, int NumNode, int Entire_num){
	printf("Classification_I start\n");
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // ���� ������ ��ġ�� ��� ���� �����ϴ� �� �������� ��
	int NumNodeGroup = 0;

	bool verify_flag = false;

	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_ciphertext_t***	cand;
	paillier_ciphertext_t***	tmp_cand;
	paillier_ciphertext_t*		tmp_coord1;
	paillier_ciphertext_t*		tmp_coord2;
	paillier_ciphertext_t*		tmp_coord3;
	paillier_ciphertext_t*		kth_dist;
	paillier_ciphertext_t**		kth_dist_bit;
	paillier_ciphertext_t*		Rand				= paillier_create_enc(5);
	paillier_ciphertext_t**		Rabel_Result		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	paillier_ciphertext_t**		tmp_nodedist_bit	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
	paillier_ciphertext_t**		alpha				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	paillier_ciphertext_t**		shortestPoint		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	paillier_ciphertext_t***	min_record			= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t***	node_dist_bit		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumNode);
	paillier_ciphertext_t***	qLL_bit				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t***	qRR_bit				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t****	nodeLL_bit			= (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	paillier_ciphertext_t****	nodeRR_bit			= (paillier_ciphertext_t****)malloc(sizeof(paillier_ciphertext_t***)*NumNode);
	
	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*3);

	for( i = 0 ; i < k ; i++) {
		Rabel_Result[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(Rabel_Result[i]->c);
		min_record[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
		for( j = 0 ; j < dim+1 ; j++ ){
			min_record[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(min_record[i][j]->c);
		}
	}
	for( i = 0 ; i < dim ; i++ ){
		shortestPoint[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(shortestPoint[i]->c);
	}
	for( i = 0 ; i < NumNode ; i++ ){
		alpha[i]			= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);

		node_dist_bit[i]	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size+1));
		
		nodeLL_bit[i]		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
		nodeRR_bit[i]		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
		
		for( j = 0 ; j < size+1 ; j++ ){
			node_dist_bit[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(node_dist_bit[i][j]->c);
		}
	}
	
	printf("===== Query -> SBD || Node -> SBD =====\n");
	
	startTime = clock();
	for( i = 0 ; i < dim ; i++ ){
		qLL_bit[i] = SBD_for_SRO(q[i], 0);
		qRR_bit[i] = SBD_for_SRO(q[i], 1);
	}

	for( i = 0 ; i < NumNode ; i++ ){
		for( j = 0 ; j < dim ; j++ ){
			nodeLL_bit[i][j] = SBD_for_SRO(node[i].LL[j], 0);
			nodeRR_bit[i][j] = SBD_for_SRO(node[i].RR[j], 1);
		}
	}
	endTime = clock();
	node_SBD_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("node SBD time: %f\n", node_SBD_time);
	
	while(1){
		startTime = clock();
		if(!verify_flag){
			printf("===== Node SRO start =====\n");
		}else{
			printf("===== Node Expansion start =====\n");
		}
		for( i = 0 ; i < NumNode ; i++ ){
			if(!verify_flag){
				alpha[i] = SRO(qLL_bit, qRR_bit, nodeLL_bit[i], nodeRR_bit[i]);
				//protocol_1_printf(alpha[i], "alpha", true);
			}else{
				for( j = 0 ; j < dim ; j++ ){
					psi[0] = SCMP(qLL_bit[j], nodeLL_bit[i][j]);
					psi[1] = SCMP(qLL_bit[j], nodeRR_bit[i][j]);
					psi[2] = SBXOR(psi[0], psi[1]);
					tmp_coord1 = SM_p1(psi[2], q[j]);
					tmp_coord2 = SM_p1(psi[0], node[i].LL[j]);
					tmp_coord3 = SM_p1(SBN(psi[0]), node[i].RR[j]);
					paillier_mul(pubkey, tmp_coord3, tmp_coord2, tmp_coord3);
					tmp_coord3 = SM_p1(SBN(psi[2]), tmp_coord3);
					paillier_mul(pubkey, shortestPoint[j], tmp_coord1, tmp_coord3);
				}
				tmp_nodedist_bit = SBD_for_SRO(SSEDm(q, shortestPoint, dim), 0);
				for(j=0; j<size+1; j++) {
					node_dist_bit[i][j] = SBOR(node_dist_bit[i][j], tmp_nodedist_bit[j]);
				}
				alpha[i] = SCMP(node_dist_bit[i], kth_dist_bit);
			}
			
			if( strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0)
				printf("%d node overlaps the query \n", i);
			
			if(!verify_flag){
				for( m = 0 ; m < size+1 ; m++ ){
					*node_dist_bit[i][m] = *alpha[i];
				}
			}
		}
		
		protocol_3_printf(node_dist_bit, NumNode, size+1, "node_dist_bit", false);
		
		if(!verify_flag){
			endTime = clock();
			node_SRO_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("SRO time : %f\n", node_SRO_time);
		}else{
			endTime = clock();
			node_expansion_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Expansion time : %f\n", node_expansion_time);
		}
		cnt = 0;
		if(!verify_flag){
			startTime = clock();
			cand = sNodeRetrievalforClassification(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_first_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node Retrieval time : %f\n", data_extract_first_time);
		}else{
			protocol_3_free(cand, cnt, dim+1);
			startTime = clock();
			tmp_cand = sNodeRetrievalforClassification(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_second_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("\nNode Retrieval time : %f\n", data_extract_second_time);

			if(cnt == 0) {	// ������ ���� �߰� Ž���� �ʿ��� ��尡 ���� ��츦 ó����
				break;
			}
			protocol_3_printf(tmp_cand, cnt, dim+1, "tmp_cand", false);
			protocol_3_printf(min_record, k, dim+1, "min_record", false);
			
			cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for( i = 0 ; i < cnt+k ; i++ ){
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
				for( j = 0 ; j < dim+1 ; j ++ ){
					cand[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
					mpz_init(cand[i][j]->c);
				}
			}
			protocol_3_printf(min_record, k, dim+1, "min_record", false);
			
			
			// ���� k���� ����� ����
			for( i = 0 ; i < k ; i++ ) {
				for( j = 0 ; j < dim+1 ; j++ ){
					cand[i][j] = min_record[i][j];
				}
			}
			// ���� ã�� �ĺ� ����� ����
			for( i = 0 ; i < cnt ; i++ ) {
				for( j = 0 ; j < dim+1 ; j++ ){
					cand[i+k][j] = tmp_cand[i][j];
				}
			}
			cnt = cnt + k;
		}
		//protocol_3_printf(cand, cnt, dim+1, "cand", true);
		
		paillier_ciphertext_t*		tmp;
		paillier_ciphertext_t*		D_min;
		paillier_ciphertext_t*		tmp_dst		= paillier_create_enc_zero();
		paillier_ciphertext_t**		dist		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t**		smin		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		paillier_ciphertext_t**		mid			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t**		V			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);
		paillier_ciphertext_t***	V2			= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);
		paillier_ciphertext_t***	SBD_dist	= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*cnt);
		

		for( i = 0 ; i < cnt ; i++ ){
			V[i]		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mid[i]		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			dist[i]		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			V2[i]		= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t)*(dim+1));
			SBD_dist[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(size));
			mpz_init(dist[i]->c);
			mpz_init(mid[i]->c);
			for( j = 0 ; j < size; j++ ){
				if( i == 0 ){
					smin[j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
				}
				SBD_dist[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
				mpz_init(SBD_dist[i][j]->c);
			}
		}
		
		printf("===== cand SSED & SBD start =====\n");
		startTime = clock();
		for( i = 0 ; i < cnt ; i++ ){
			dist[i] = SSEDm(q, cand[i], dim);
			SBD_dist[i] = SBD(dist[i]);
		}
		endTime = clock();
		data_SSED_SBD_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		printf("data SSED & SBD time : %f\n", (float)(endTime-startTime)/(CLOCKS_PER_SEC));

		protocol_2_printf(dist, cnt, "dist", false);
		protocol_3_printf(SBD_dist, cnt, size, "SBD_dist", false);
		
		for( s = 0 ; s < k ; s++ ){
			printf("\n%dth sMINn start \n", s+1);
			startTime = clock();
			smin = Smin_n(SBD_dist, cnt);
			endTime = clock();
			gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			if(!verify_flag){
				sMINn_first_time += gap;
				printf("%dth sMINn time : %f\n", s+1, gap);
			}
			else{
				sMINn_second_time += gap;
				printf("%dth sMINn time : %f\n", s+1, gap);
			}
			protocol_2_printf(smin, size, "smin", false);
			if( s == 0 ){
				D_min = paillier_create_enc_zero();
			}else{
				protocol_1_free(D_min);
				D_min = paillier_create_enc_zero();	
			}
			for( j = size ; j > 0 ; j-- ){
				t = (int)pow(2, j-1);
				tmp = paillier_create_enc(t);
				tmp = SM_p1(tmp, smin[size-j]);
				paillier_mul(pubkey, D_min, tmp, D_min);
				protocol_1_free(tmp);
			}
			protocol_1_printf(D_min, "D_min", false);
			printf("== recalculate query<->data distances ===\n");
			if( s != 0 ){
				for( i = 0 ; i < cnt; i++ ){
					for( j = size ; j > 0 ; j--){
						t = (int)pow(2, j-1);
						tmp = paillier_create_enc(t);
						tmp = SM_p1(tmp, SBD_dist[i][size-j]);
						paillier_mul(pubkey, tmp_dst, tmp, tmp_dst);
					}
					dist[i] = tmp_dst;
					tmp_dst = paillier_enc(0, pubkey, plain_zero, paillier_get_rand_devurandom);
				}
			}
			for( i = 0 ; i < cnt ; i++ ){
				mpz_init(mid[i]->c);
				paillier_subtract(pubkey, mid[i], dist[i], D_min);
				mid[i] = SM_p1(mid[i], Rand);
			}
			protocol_2_printf(mid, cnt, "mid", false);
			V = SkNNm_sub(mid, cnt);

			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0; j < dim+1 ; j++){
					V2[i][j] = SM_p1(V[i], cand[i][j]);
					if( i == 0 ){
						min_record[s][j] = V2[i][j];
					}else{
						paillier_mul(pubkey, min_record[s][j], V2[i][j], min_record[s][j]);
					}
				}
			}
			startTime = clock();
			for( i = 0 ; i < cnt ; i++ ){
				for( j = 0; j < size ; j++ ){
					SBD_dist[i][j] = SBOR(V[i], SBD_dist[i][j]);
				}
			}
			endTime = clock();
			data_SBOR_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			protocol_3_printf(SBD_dist, cnt, size, "SBD dist", false);
		}

		if(!verify_flag){
			kth_dist	= SSEDm(q, min_record[k-1], dim);
			kth_dist_bit = SBD_for_SRO(kth_dist, 0);
			protocol_1_printf(kth_dist, "kth_dst", true);
			protocol_2_printf(kth_dist_bit, size, "kth_dist_bit", false);
		}
		

		protocol_1_free(tmp);
		protocol_1_free(D_min);
		protocol_1_free(tmp_dst);
		protocol_2_free(dist, cnt);
		protocol_2_free(smin, size);
		protocol_2_free(mid, cnt);
		protocol_2_free(V, cnt);
		//protocol_3_free(V2, cnt, dim+1);
		protocol_3_free(SBD_dist, cnt, size);

		if(!verify_flag)	{
			verify_flag = true;		// ������ �����ϱ� ���� flag�� true�� ����
		}else{
			break;
		}
	}
	protocol_3_printf(min_record, k, dim+1, "min_record", true);
	for( i = 0 ; i < k ; i++ ){
		*Rabel_Result[i] = *min_record[i][dim];
	}

	paillier_plaintext_t* A = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	int** sknn = (int**)malloc(sizeof(int*)*k);
	for( i = 0 ; i < k ; i++ ){
		sknn[i] = (int*)malloc(sizeof(int)*(dim+1));
		for( j = 0 ; j < dim+1 ; j++ ){
			A = paillier_dec(0, pubkey, prvkey, min_record[i][j]);
			sknn[i][j] = mpz_get_ui(A->m);
			printf("%d ", sknn[i][j]);
		}
		printf("\n");
	}

	protocol_2_free(alpha, NumNode);
	protocol_2_free(shortestPoint, dim);
	protocol_3_free(cand, cnt, dim+1);
	//protocol_3_free(node_dist_bit, NumNode, size+1);
	//protocol_3_free(qLL_bit, dim, size+1);
	//protocol_3_free(qRR_bit, dim, size+1);
	//protocol_4_free(nodeLL_bit, NumNode, dim, size+1);
	//protocol_4_free(nodeRR_bit, NumNode, dim, size+1);
	
	protocol_2_printf(Rabel_Result, k, "Rabel", true);
	
	SCMC(Entire_set, Rabel_Result, Entire_num, k);
	return sknn;
}

int** protocol::Classification_M(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, paillier_ciphertext_t** Entire_set, int k, int row_number, int Entire_num){
	printf("Classification_m start\n");
	int i=0, s=0, n=0, t=0, j=0;
	int rand = 5;
	n = row_number;
	
	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* ciper_binary;
	paillier_ciphertext_t* ciper_min				= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_dist				= paillier_create_enc_zero();
	paillier_ciphertext_t* ciper_rand				= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_dist				= paillier_create_enc_zero();
	
	paillier_ciphertext_t*** min_record				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t** ciper_distance			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t*** ciper_SBD_distance		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*n);
	paillier_ciphertext_t** ciper_mid				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** ciper_Smin				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t** ciper_V					= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*n);
	paillier_ciphertext_t*** ciper_V2				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*n);

	
	paillier_ciphertext_t** ciper_rabel_result = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	for( i = 0 ; i < size; i++ ){
		ciper_Smin[i] = ciper_zero;
	}

	for( i = 0 ; i < n ; i++ ){
		ciper_distance[i] 	= ciper_zero;
		ciper_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		ciper_V[i]	= ciper_zero;
		ciper_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
	}

	for( i = 0 ; i < k ; i++ ){
		min_record[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim ; j++ ){
			min_record[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(min_record[i][j]->c);
		}
		ciper_rabel_result[i] = paillier_create_enc_zero();
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
		//printf("\n%dth sMINn start \n", s+1);
		startTime = clock();
		ciper_Smin = Smin_n(ciper_SBD_distance, n);	// bit�� ǥ���� ��ȣȭ min �Ÿ� ����
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
		
		// bit�� ǥ���� ��ȣȭ min �Ÿ��� ��ȣȭ ������ ��ȯ
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			ciper_binary = paillier_create_enc(t);
			ciper_binary = SM_p1(ciper_binary, ciper_Smin[size-j]);
			paillier_mul(pubkey, ciper_min, ciper_binary, ciper_min);
			paillier_freeciphertext(ciper_binary);
			//paillier_print("min : ", ciper_min);
		}
		paillier_print("min dist : ", ciper_min);
		
		//printf("\n== recalculate query<->data distances ===\n");
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
		
		/*
		for( i = 0 ; i < n ; i++ ){
			paillier_print("distance : ", ciper_distance[i]);
		}
		*/

		// ����-������ �Ÿ��� min �Ÿ����� ���� ���� (min �������� ��쿡�� 0���� ����� ����)
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

		// min ������ ����
		for( i = 0 ; i < n ; i++ ){
			for( j = 0 ; j < dim+1 ; j++ ){
				ciper_V2[i][j] = SM_p1(ciper_V[i], data[i][j]);
				if( i == 0 ){
					min_record[s][j] = ciper_V2[i][j];
				}else{
					paillier_mul(pubkey, min_record[s][j], ciper_V2[i][j], min_record[s][j]);
				}
			}
		}
		// Data SBOR ����
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

	protocol_3_printf(min_record, k, dim+1, "min_record", true);
	for( i = 0 ; i < k ; i++ ){
		*ciper_rabel_result[i] = *min_record[i][dim];
	}

	paillier_plaintext_t* A = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	int** sknn = (int**)malloc(sizeof(int*)*k);
	for( i = 0 ; i < k ; i++ ){
		sknn[i] = (int*)malloc(sizeof(int)*(dim+1));
		for( j = 0 ; j < dim+1 ; j++ ){
			A = paillier_dec(0, pubkey, prvkey, min_record[i][j]);
			sknn[i][j] = mpz_get_ui(A->m);
			printf("%d ", sknn[i][j]);
		}
		printf("\n");
	}

	paillier_freeciphertext(temp_dist);	
	SCMC(Entire_set, ciper_rabel_result, Entire_num, k);

	return sknn;
}

paillier_ciphertext_t** protocol::SF_P1(paillier_ciphertext_t** Entire_set, paillier_ciphertext_t** K_set, int w, int k){
	printf("Secure Frequency start\n");
	int i = 0;
	int j = 0;
	int Random = 2;
	
	paillier_ciphertext_t** T = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(k));
	paillier_ciphertext_t*** S = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(k));
	paillier_ciphertext_t** V = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(w));
	for( i = 0 ; i < w ; i++ )
	{
		V[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(V[i]->c);
	}
	for( i = 0 ; i < k ; i++ )
	{
		T[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(T[i]->c);
		S[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(w));
		for( j = 0 ; j < w ; j++ )
		{
			S[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(S[i][j]->c);
		}
	}
	
	for( i = 0 ; i < k ; i++ ){
		for ( j = 0 ; j < w ; j++ )
		{
			paillier_subtract(pub, S[i][j], Entire_set[j], K_set[i]);
		}
	}

	S = SF_P2(S, w, k);
	
	for( i = 0 ; i < w ; i++ )
	{
		//printf("%d : ", i+1);
		for( j = 0 ; j < k ; j++ )
		{
			if( j == 0 )
			{
				V[i] = S[j][i];
			}
			else
			{
				paillier_mul(pub, V[i], V[i], S[j][i]);
			}
		}
		//gmp_printf( " %Zd \n", paillier_dec(0, pub, prv, V[i]));
	}
	//free
	
	return V;
}
paillier_ciphertext_t*** protocol::SF_P2(paillier_ciphertext_t*** S, int w, int k){
	int i = 0;
	int j = 0;

	paillier_ciphertext_t*** U = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(k));
	for( i = 0 ; i < k ; i++ )
	{
		U[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(w));
		for( j = 0 ; j < w ; j++ )
		{
			U[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(U[i][j]->c);
		}
	}
	for( i = 0 ; i < k ; i++ )
	{
		for( j = 0 ; j < w ; j++ )
		{
			if(strcmp(paillier_plaintext_to_str(paillier_dec(0,pubkey,prvkey,S[i][j])),paillier_plaintext_to_str(plain_zero)) == 0 )
			{
				U[i][j] = paillier_create_enc_one();
			}
			else
			{
				U[i][j] = paillier_create_enc_zero();
			}
		}
	}
	return U;
}
paillier_ciphertext_t* protocol::SCMC(paillier_ciphertext_t** Entire_set, paillier_ciphertext_t** K_set, int w, int k){
	printf("===SCMS start====\n");
	int i = 0;
	protocol_2_printf(Entire_set, w, "Entire_set", true);
	protocol_2_printf(K_set, k, "K_set", true);
	int idx = 0;
	paillier_ciphertext_t** V = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(w));
	V = SF_P1(Entire_set, K_set, w, k);
	protocol_2_printf(V, w, "V", false);
	idx = MAXn(V, w);
	free(V);
	return Entire_set[idx];
}
