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

void protocol::kNN_SSED_SBD_inMultithread(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_distance)
{
	cout << "kNN_SSED_SBD_inMultithread" << endl;
	int threadNum = thread_num;
	if(threadNum > NumData)
	{
		threadNum = NumData;
	}

	std::vector<int> *inputIdx = new std::vector<int>[threadNum];
	int count = 0;

	for(int i = 0; i < NumData; i++)
	{
		inputIdx[count++].push_back(i);
		count %= threadNum;
	}
	
	
	std::thread *SSED_SBD_thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		SSED_SBD_thread[i] = std::thread(SSED_SBD_inThread, i, std::ref(inputIdx[i]), std::ref(query), std::ref(data), dim,	std::ref(cipher_distance), std::ref(cipher_SBD_distance), std::ref(*this));
	}

	for(int i = 0; i < threadNum; i++)
	{
		SSED_SBD_thread[i].join();
	}
	
	delete[] SSED_SBD_thread;
	delete[] inputIdx;
}

void protocol::recalculate_DISTforkNN_inMultithread(paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_distance)
{
	cout << "recalculate_DISTforkNN_inMultithread" << endl;
	int threadNum = thread_num;
	if(threadNum > NumData)
	{
		threadNum = NumData;
	}

	std::vector<int> *inputIdx = new std::vector<int>[threadNum];
	int count = 0;

	for(int i = 0; i < NumData; i++)
	{
		inputIdx[count++].push_back(i);
		count %= threadNum;
	}
	
	std::thread *recalculate_DISTforkNN_Thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		recalculate_DISTforkNN_Thread[i] = std::thread(recalculate_DISTforkNN_inThread, std::ref(inputIdx[i]), std::ref(cipher_SBD_distance), std::ref(cipher_distance), std::ref(*this));
	}

	for(int i = 0; i < threadNum; i++)
	{
		recalculate_DISTforkNN_Thread[i].join();
	}
	
	delete[] recalculate_DISTforkNN_Thread;
	delete[] inputIdx;
}


void protocol::UPDATE_SBD_SCORE_InKNN_PB_inMultithread(paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_V, int NumData)
{
	cout << "UPDATE_SBD_SCORE_InKNN_PB_inMultithread" << endl;
	int threadNum = thread_num;
	if(threadNum > NumData)
	{
		threadNum = NumData;
	}

	std::vector<int> *inputIdx = new std::vector<int>[threadNum];
	int count = 0;

	for(int i = 0; i < NumData; i++)
	{
		inputIdx[count++].push_back(i);
		count %= threadNum;
	}
	std::thread *UPDATE_SBD_SCORE_InKNN_PB_Thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		UPDATE_SBD_SCORE_InKNN_PB_Thread[i] = std::thread(UPDATE_SBD_SCORE_InKNN_PB_inThread, std::ref(inputIdx[i]), std::ref(cipher_SBD_distance), std::ref(cipher_V), std::ref(*this));
	}
	for(int i = 0; i < threadNum; i++)
	{
		UPDATE_SBD_SCORE_InKNN_PB_Thread[i].join();
	}

	delete[] UPDATE_SBD_SCORE_InKNN_PB_Thread;
	delete[] inputIdx;
}


int** protocol::SkNN_PB(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, int k, int NumData)
{
	cout << "\n=========================== SkNN_PB start ============================\n" << endl;

	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}

	int i=0, s=0, n=0, t=0, j=0;
	int rand = 5;
	n = NumData;
	
	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* cipher_binary;
	paillier_ciphertext_t* cipher_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_rand	= paillier_create_enc(rand);

	paillier_ciphertext_t** cipher_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);
	paillier_ciphertext_t*** cipher_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	paillier_ciphertext_t** cipher_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);
	paillier_ciphertext_t** cipher_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);
	paillier_ciphertext_t** cipher_V=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);


	paillier_ciphertext_t*** cipher_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	for( i = 0 ; i < size; i++ ){
		cipher_Smin[i] = cipher_zero;
	}

	for( i = 0 ; i < NumData ; i++ ){
		cipher_distance[i] 	= cipher_zero;
		cipher_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		cipher_V[i]	= cipher_zero;
	}

	for( i = 0 ; i < k ; i++ ){
		cipher_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( j = 0 ; j < dim; j ++ ){
			cipher_result[i][j] = paillier_create_enc_zero();
		}
	}

	kNN_SSED_SBD_inMultithread(data, query, cipher_SBD_distance, cipher_distance);

	for( s = 0 ; s < k ; s++ ){
		printf("\n%dth sMINn start \n", s+1);

		startTime = clock();
		cipher_Smin = Smin_n_Multithread(cipher_SBD_distance, NumData); // bit로 표현된 암호화 min 거리 추출
		endTime = clock();
		gap = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
		sMINn_first_time += gap;
		printf("%dth sMINn time : %f\n", s+1, gap);	
	
		// bit로 표현된 암호화 min 거리를 암호화 정수로 변환
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			cipher_binary = paillier_create_enc(t);
			cipher_binary = SM_p1(cipher_binary, cipher_Smin[size-j]);
			paillier_mul(pubkey, cipher_min, cipher_binary, cipher_min);
			paillier_freeciphertext(cipher_binary);
		}
		paillier_print("min dist : ", cipher_min);
		
		printf("\n== recalculate query<->data distances ===\n");
		// 이전 iteration에서 min값으로 선택된 데이터의 거리가 secure하게 MAX로 변환되었기 때문에, 
		// 질의-데이터 간 거리 계산을 모두 다시 수행
		if( s != 0 )
		{
			recalculate_DISTforkNN_inMultithread(cipher_SBD_distance, cipher_distance);
		}
		
		// 질의-데이터 거리와 min 거리와의 차를 구함 (min 데이터의 경우에만 0으로 만들기 위함)
		DuringProcessedInTOPK_PB_inMultithread(cipher_distance, cipher_mid, cipher_min, cipher_rand, NumData);
		cipher_V = SkNNm_sub(cipher_mid, NumData);
		// min 데이터 추출
		ExtractTOPKInTOPK_PB_inMultithread(data, cipher_V, cipher_result, NumData, s);
		printf("cipher_result : ");	
		for( j = 0 ; j < dim; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey,cipher_result[s][j]));
		}
		printf("\n");		

		// Data SBOR 수행
		UPDATE_SBD_SCORE_InKNN_PB_inMultithread(cipher_SBD_distance, cipher_V, NumData);
		cipher_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
	}
	
	// user(Bob)에게 결과 전송을 위해 random 값 삽입
	for(i=0;i<k;i++){
		for(j=0; j<dim; j++){
			paillier_mul(pubkey,cipher_result[i][j],cipher_result[i][j],cipher_rand);
		}
	}


	
	return SkNNm_Bob2(cipher_result, rand, k, dim);
}



int** protocol::SkNN_PGI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode)
{
	printf("\n=== SkNN_PGI start ===\n");
	int rand = 5;
	int cnt = 0;	
	int NumNodeGroup = 0;
	//Print = true;
	bool verify_flag = false;

	std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point end;
    //std::chrono::duration<float> duration_sec;
	
	float progress;
	
	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;

	paillier_ciphertext_t **alpha = new paillier_ciphertext_t *[NumNode];
	for(int i = 0; i < NumNode; i++)
	{
		alpha[i] = new paillier_ciphertext_t;
		mpz_init(alpha[i]->c);
	}
	paillier_ciphertext_t ***cand;
	paillier_ciphertext_t ***temp_cand;

	paillier_ciphertext_t ***Result = new paillier_ciphertext_t **[k];

	for(int i = 0; i < k; i++)
	{
		Result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim));
		for(int j = 0 ; j < dim; j++)
		{
			Result[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(Result[i][j]->c);
		}
	}

	paillier_ciphertext_t*		C_RAND			= paillier_create_enc(rand);
	paillier_ciphertext_t*		TMP_alpha		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	paillier_ciphertext_t*		MAX				= paillier_create_enc(pow(2, size)-1);
	paillier_ciphertext_t*		K_DIST = 0;	

	while(1)
	{
		progress = 0.1;
		startTime = clock();
		start = std::chrono::system_clock::now();

		if(!verify_flag)
		{	
			std::cout << "!!!!!!!!GSRO Start!!!!!!!!" << std::endl;
			Parallel_GSRO_kNN(q, q, alpha, node, NumNode, &cnt, &NumNodeGroup, verify_flag);
		}
		else
		{
			std::cout << "!!!!!!!!Seconde Node Expansion!!!!!!!!" << std::endl;

			int NumThread = thread_num;
			if(NumThread > NumNode)
			{
				NumThread = NumNode;
			}

			std::cout << "numThread: " << NumThread << std::endl;

			std::vector<int> *NumNodeInput = new std::vector<int>[NumThread];
			int count = 0;
			for(int i = 0; i < NumNode; i++)
			{
				NumNodeInput[count++].push_back(i);
				count %= NumThread;
			}

			std::thread *NodeExpansionThread = new std::thread[NumThread];
			for(int i = 0; i < NumThread; i++)
			{
				NodeExpansionThread[i] = std::thread(NodeExpansion_Multithread, std::ref(q), std::ref(K_DIST), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), i);
			}

			for(int i = 0; i < NumThread; i++)
			{
				NodeExpansionThread[i].join();
			}

			delete[] NumNodeInput;
			delete[] NodeExpansionThread;
			
		}
		
		
		cnt = 0;

		if(!verify_flag)
		{	
			cand = sNodeRetrievalforClassificationMultiThread(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			printf("Node Retrieval time : %f\n", data_extract_first_time);
		}else{
			temp_cand = sNodeRetrievalforClassificationMultiThread(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			printf("Node Retrieval time : %f\n", data_extract_second_time);

			if(cnt == 0)
			{	
				break;
			}

			cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for(int i = 0 ; i < cnt+k ; i++ )
			{
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
				for(int j = 0 ; j < dim ; j ++ )
				{
					cand[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
					mpz_init(cand[i][j]->c);
				}
			}

			for(int i = 0 ; i < k ; i++ )
			{
				for(int j = 0 ; j < dim ; j ++ )
				{
					cand[i][j] = Result[i][j];
				}
			}
			for(int i = 0 ; i < cnt ; i++ )
			{
				for(int j = 0 ; j < dim ; j ++ )
				{
					cand[i+k][j] = temp_cand[i][j];
				}
			}
			cnt = cnt + k;
		}

		if(Print)
		{
			std::cout << "cnt : "<<cnt<<std::endl;
			std::cout << "!!!!EXTRACT CAND LIST!!!!"<<std::endl;
			for(int i = 0 ; i < cnt ; i ++)
			{
				for(int j = 0 ; j < dim ; j ++ )
				{
					gmp_printf("%Zd\t", paillier_dec(0, pubkey, prvkey, cand[i][j]));
				}
				std::cout<<std::endl;
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

		for(int i = 0 ; i < cnt ; i ++ )
		{
			DIST_MINUS_MIN[i]	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			DIST[i]				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			V[i]				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			V2[i]				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim));		
			mpz_init(DIST_MINUS_MIN[i]->c);
			mpz_init(DIST[i]->c);
			mpz_init(V[i]->c);
			for(int j = 0 ; j < dim ; j++ ){
				V2[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
				mpz_init(V2[i][j]->c);
			}	
		}
		
		SSEDMultiThread(DIST, cand, q, cnt);

		
		for(int s = 0 ; s < k ; s ++ )
		{
			std::cout << s+1 <<" th Classification Start !!!!!!!!!!!!"<<std::endl;
			idx = SMSn(DIST, cnt);
			MIN = DIST[idx];
			SMSnMultithread2(cnt, DIST_MINUS_MIN, DIST, MIN, C_RAND);
			V = SkNNm_sub(DIST_MINUS_MIN, cnt);
			SMSnMultithread3(cnt, s, V, V2, DIST, cand, MAX, Result);
		}
		
		if(!verify_flag)
		{	
			K_DIST = SSEDm(q, Result[k-1], dim); 
		}

		if(!verify_flag)
			verify_flag = true;
		else
			break;
	}

	
	paillier_plaintext_t* A = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	int** sknn = (int**)malloc(sizeof(int*)*k);
	for(int i = 0 ; i < k ; i++ ){
		sknn[i] = (int*)malloc(sizeof(int)*(dim));
		for(int j = 0 ; j < dim ; j++ ){
			A = paillier_dec(0, pubkey, prvkey, Result[i][j]);
			sknn[i][j] = mpz_get_ui(A->m);
			printf("%d ", sknn[i][j]);
		}
		printf("\n");
	}
	return SkNNm_Bob2(Result, 0, k, dim);
}


int** protocol::SkNN_PAI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, int k, int NumData, int NumNode)
{
	return 0;
}

void protocol::Parallel_GSRO_kNN(paillier_ciphertext_t** cipher_qLL, paillier_ciphertext_t** cipher_qRR, paillier_ciphertext_t** alpha, boundary* node, int NumNode, int* cnt, int* NumNodeGroup, bool type)
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
		GSROThread[i] = std::thread(GSRO_Multithread, std::ref(cipher_qLL), std::ref(cipher_qLL), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), true);
	}
	for(int i = 0; i < NumThread; i++)
	{
		GSROThread[i].join();
	}
	delete[] NumNodeInput;
	delete[] GSROThread;
}


