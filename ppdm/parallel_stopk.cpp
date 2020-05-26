#include <iostream>
#include <fstream>
#include <gmp.h>
#include <string>
#include <sstream>
#include <thread>
#include <chrono>
#include <vector>
#include <iostream>
#include "protocol.h"
#include "parallel_protocol.h"
using namespace std;


int** protocol::STopk_PB(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode)
{
	return 0;
}

int** protocol::STopk_PGI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode)
{
	printf("\n=== STopk_PGI start ===\n");
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;
	time_t startTime = 0;	time_t endTime = 0;	float gap = 0.0;

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

	Parallel_GSRO_Topk(Q, Q, alpha, node, NumNode, &cnt, &NumNodeGroup, true);
///////////////////////////////////////////////////////////////
/*
	startTime = clock();
	for( i = 0 ; i < NumNode ; i++ )
	{	
		alpha[i] = DP_GSRO(Q, Q, node[i].LL, node[i].RR);
		
		if(!verify_flag)
		{	//  검증 단계에서는 수행하지 않음
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
*/
///////////////////////////////////////////////////////////////
	while(1)
	{
		cnt = 0;
		if(!verify_flag)
		{	
			startTime = clock();

			cand = Parallel_GSRO_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
///////////////////////////////////////////////////////////////
			//cand = GSRO_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
///////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////
			temp_cand = Parallel_GSRO_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
///////////////////////////////////////////////////////////////
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = clock();
			data_extract_second_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			printf("Node expansion time : %f\n", data_extract_second_time);
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

		ComputeScoreinMultithread(cand, q, SCORE, &cnt, true);
///////////////////////////////////////////////////////////////
/*
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
*/
///////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////
			MAXnMultithread2(cnt, SCORE_MINUS_MAX, SCORE, MAX, C_RAND);
/*
			for( i = 0 ; i < cnt ; i++ )
			{
				paillier_subtract(pubkey, SCORE_MINUS_MAX[i], SCORE[i], MAX);
				SCORE_MINUS_MAX[i] = SM_p1(SCORE_MINUS_MAX[i], C_RAND);
				if(Print)
				{
					gmp_printf("SCORE-MAX : %Zd \n", paillier_dec(0, pubkey, prvkey, SCORE_MINUS_MAX[i]));
				}
			}
*/
///////////////////////////////////////////////////////////////////
			V = Topk_sub(SCORE_MINUS_MAX, cnt);
///////////////////////////////////////////////////////////////////
			MAXnMultithread3(cnt, s, V, V2, SCORE, cand, MIN, RESULT);
/*
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
*/
///////////////////////////////////////////////////////////////////
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
			if(cnt < k)	 // 노드의 fanout이 k보다 적은 경우를 처리함
				kth_score = computeScore(q, RESULT[cnt-1]);
			else	// k개가 찾아졌다면, k번째 데이터의 score를 계산
				kth_score = computeScore(q, RESULT[k-1]);
			//gmp_printf("%dth score : %Zd\n", k, paillier_dec(0, pubkey, prvkey, kth_score));
///////////////////////////////////////////////////////////////////
			for( i = 0 ; i < NumNode ; i++ )
			{	
				for( j = 0 ; j < dim ; j++ )
				{
					paillier_mul(pubkey, Q[j], SM_p1(node[i].RR[j], psi[j]), SM_p1(node[i].LL[j], SBN(psi[j])));	
				}
				alpha[i] = GSCMP(kth_score , computeScore(q, Q));	// k번째 score보다 높을 가능성이 있으면 alpha=1
				
				if(Print)
				{
					gmp_printf("computeScore(q, Q) : %Zd %d line\n", paillier_dec(0, pubkey, prvkey, computeScore(q, Q)), __LINE__);
					gmp_printf("alpha : %Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
				}
			}
			endTime = clock();
			node_expansion_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
///////////////////////////////////////////////////////////////////
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
		printf("\n");
	}
	cout << "End Line" <<endl;
	return SkNNm_Bob2(RESULT, rand, k, dim);
}

int** protocol::STopk_PAI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode)
{
	return 0;
}



paillier_ciphertext_t*** protocol::Parallel_GSRO_sNodeRetrievalforTopk(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now Parallel_GSRO_sNodeRetrievalforTopk starts =====\n");
	printf("NumNode : %d thread_num : %d\n", NumNode, thread_num);
	int i=0, j=0, m=0;

	int** node_group;
	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);


	paillier_ciphertext_t*** cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(*NumNodeGroup)*FanOut);		
	for(int i = 0 ; i < *NumNodeGroup * FanOut ; i++ ) {
		cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j ++ ){
			cand[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cand[i][j]->c);
		}
	}

	// Extract Data
	int threadNum = thread_num;
	if(threadNum > NumNode)
	{
		threadNum = NumNode;
	}

	int remained = 0;

	float progress = 0.1;

	std::vector<int*> *inputNodeGroup = new std::vector<int *>[threadNum];
	std::vector<int> *inputJ = new std::vector<int>[threadNum];
	std::vector<int> *inputCnt = new std::vector<int>[threadNum];
	int currentCount = 0;
	for(int i = 0; i < *NumNodeGroup; i++)
	{
		remained = node_group[i][0];
		int *nodeGroup = node_group[i];
		for(int j = 0; j < FanOut; j++)
		{
			inputNodeGroup[currentCount].push_back(nodeGroup);
			inputJ[currentCount].push_back(j);
			inputCnt[currentCount].push_back(*cnt);
			(*cnt)++;
			currentCount = (currentCount + 1) % threadNum;
		}
	}

	std::thread *sNodeRetrievalforDataThread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		sNodeRetrievalforDataThread[i] = std::thread(sNodeRetrievalforDatainThread, std::ref(inputJ[i]), std::ref(inputNodeGroup[i]), std::ref(remained), std::ref(inputCnt[i]), std::ref(node), std::ref(data), std::ref(alpha), std::ref(cand), std::ref(*this), i);
	}

	for(int i = 0; i < threadNum; i++)
	{
		sNodeRetrievalforDataThread[i].join();
	}

	delete[] inputNodeGroup;
	delete[] inputJ;
	delete[] inputCnt;
	delete[] sNodeRetrievalforDataThread;

	return cand;
}

void protocol::ComputeScoreinMultithread(paillier_ciphertext_t*** cand, paillier_ciphertext_t** q, paillier_ciphertext_t** SCORE, int * cnt, bool type)
{
	int NumThread = thread_num;
	if(NumThread > (*cnt) )
	{
		NumThread = (*cnt);
	}
	std::vector<int> *CSInput = new std::vector<int>[NumThread];
	int count = 0;
	for(int i = 0; i < (*cnt) ; i++)
	{
		CSInput[count++].push_back(i);
		count %= NumThread;
	}
	std::thread *CSThread = new std::thread[NumThread];
	for(int i = 0; i < NumThread; i++)
	{
		CSThread[i] = std::thread(ComputeScore_Multithread, std::ref(cand), std::ref(q), std::ref(SCORE), std::ref(CSInput[i]), std::ref(*this), true);
	}
	for(int i = 0; i < NumThread; i++)
	{
		CSThread[i].join();
	}
	delete[] CSInput;
	delete[] CSThread;
}

void protocol::Parallel_GSRO_Topk(paillier_ciphertext_t** cipher_qLL, paillier_ciphertext_t** cipher_qRR, paillier_ciphertext_t** alpha, boundary* node, int NumNode, int* cnt, int* NumNodeGroup, bool type)
{
	int NumThread = thread_num;
	if(NumThread > NumNode )
	{
		NumThread = NumNode;
	}
	std::vector<int> *NumNodeInput = new std::vector<int>[NumThread];
	int count = 0;
	for(int i = 0; i < NumNode ; i++)
	{
		NumNodeInput[count++].push_back(i);
		count %= NumThread;
	}
	std::thread *GSROThread = new std::thread[NumThread];
	for(int i = 0; i < NumThread; i++)
	{
		GSROThread[i] = std::thread(DP_GSRO_MultithreadforTopk, std::ref(cipher_qLL), std::ref(cipher_qLL), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), true);
	}
	for(int i = 0; i < NumThread; i++)
	{
		GSROThread[i].join();
	}
	delete[] NumNodeInput;
	delete[] GSROThread;
}

void protocol::MAXnMultithread2(int cnt, paillier_ciphertext_t** DIST_MINUS_MIN, paillier_ciphertext_t** DIST, paillier_ciphertext_t* MIN, paillier_ciphertext_t* C_RAND)
{
	int numThread = thread_num;
	if(numThread > cnt)
	{
		numThread = cnt;
	}

	std::vector<int> *inputCnt = new std::vector<int> [numThread];
	int count = 0;
	for(int i = 0; i < cnt; i++)
	{
		inputCnt[count++].push_back(i);
		count %= numThread;
	}

	std::thread *SMSThread = new std::thread[numThread];
	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i] = std::thread(MAXnInthread2, std::ref(inputCnt[i]), std::ref(DIST_MINUS_MIN), std::ref(DIST), std::ref(MIN), std::ref(C_RAND), std::ref(*this), i);
	}

	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i].join();
	}

	delete[] inputCnt;
	delete[] SMSThread;
}

void protocol::MAXnMultithread3(int cnt, int s, paillier_ciphertext_t **V, paillier_ciphertext_t ***V2, paillier_ciphertext_t **SCORE, paillier_ciphertext_t ***cand, paillier_ciphertext_t *MIN, paillier_ciphertext_t ***Result)
{
	int numThread = thread_num;
	if(numThread > cnt)
	{
		numThread = cnt;
	}
	
	std::vector<int> *inputCnt = new std::vector<int>[numThread];
	int count = 0;
	for(int i = 0; i < cnt; i++)
	{
		inputCnt[count++].push_back(i);
		count %= numThread;
	}

	std::thread *SMSThread = new std::thread[numThread];
	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i] = std::thread(MAXnInThread3_1, i, s, std::ref(inputCnt[i]), std::ref(V), std::ref(SCORE), std::ref(cand), std::ref(MIN), std::ref(*this));
	}
	
	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i].join();
	}
	delete[] SMSThread;
	delete[] inputCnt;

	SMSThread = new std::thread[numThread];
	inputCnt = new std::vector<int>[numThread];
	int offset[numThread];
	int record_size[numThread];
	int quot = (int)(cnt/numThread);
	int remain = (int)(cnt%numThread);

	paillier_ciphertext_t*** thread_result = new paillier_ciphertext_t **[numThread];
	for(int i=0; i<numThread; i++)
	{
		record_size[i] = quot;
		thread_result[i] = new paillier_ciphertext_t*[dim];
                for(int j = 0; j < dim; j++)
                {
                        thread_result[i][j] = paillier_create_enc_zero();
                }
		if(i<remain)
		{
			record_size[i]++;
		}
	}
	offset[0] = 0;
	for(int i=1; i<numThread; i++)
	{
		offset[i] = record_size[i-1] + offset[i-1];
	}

	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i] = std::thread(MAXnInThread3_2_v2, i, offset[i], record_size[i], dim, std::ref(thread_result[i]), std::ref(V), std::ref(V2), std::ref(SCORE), std::ref(cand), std::ref(MIN), std::ref(*this));
	}
	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i].join();
	}

	for(int i=0; i<numThread; i++)
	{
		for(int j=0; j<dim; j++)
		{
			if(i == 0)
			{
				Result[s][j] = thread_result[i][j];
			}
			else
			{
				paillier_mul(pubkey, Result[s][j], thread_result[i][j], Result[s][j]);
			}
		}
	}

	delete[] SMSThread;
	delete[] inputCnt;
}