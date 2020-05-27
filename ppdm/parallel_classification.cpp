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

paillier_ciphertext_t** protocol::Smin_n_Multithread(paillier_ciphertext_t*** cipher, int number)
{
	paillier_ciphertext_t*** copy_cipher = new paillier_ciphertext_t **[number];
	for(int i = 0; i < number ; i++)
	{
		copy_cipher[i] = cipher[i];
	}

	paillier_ciphertext_t **min;
	paillier_ciphertext_t **sbd = SBD(paillier_create_enc_zero());
	
	int numThread = thread_num;
	if(numThread > number)
	{
		numThread = number;
	}

	std::vector<int> *inputIdx = new std::vector<int>[numThread];
	paillier_ciphertext_t ***inputMin = new paillier_ciphertext_t **[numThread];

	int count = 0;
	for(int i = 0; i < number; i++)
	{

		inputIdx[count++].push_back(i);
		count %= numThread;
	}

	std::thread *sminThread = new std::thread[numThread];
	for(int i = 0; i < numThread; i++)
	{
		sminThread[i] = std::thread(Smin_n_InThread, i, std::ref(inputIdx[i]), std::ref(cipher), std::ref(inputMin), std::ref(*this));
	}

	for(int i = 0; i < numThread; i++)
	{
		sminThread[i].join();
	}

	min = Smin_n(inputMin, numThread);

	delete[] inputIdx;
	delete[] sminThread;
	delete[] copy_cipher;
	delete sbd;

	return min;
}

void protocol::SSEDMultiThread(paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t **q, int cnt)
{
	int NumThread = thread_num;
	if(thread_num > cnt)
	{
		NumThread = cnt;
	}

	std::vector<int> *inputCnt = new std::vector<int> [NumThread];
	int count = 0;
	for(int i = 0; i < cnt; i++)
	{
		inputCnt[count++].push_back(i);
		count %= NumThread;
	}

	std::thread *SSEDThread = new std::thread[NumThread];
	for(int i = 0; i < NumThread; i++)
	{
		SSEDThread[i] = std::thread(SSEDInThread, std::ref(inputCnt[i]), std::ref(DIST), std::ref(cand), std::ref(q), cnt, std::ref(*this));
	}

	for(int i = 0; i < NumThread; i++)
	{
		SSEDThread[i].join();
	}

	delete[] inputCnt;
	delete[] SSEDThread;
}

paillier_ciphertext_t*** protocol::sNodeRetrievalforClassificationMultiThread(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha,	int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now sNodeRetrievalforClassification starts =====\n");
	
	int **node_group;

	node_group = sRange_sub(alpha, NumNode, NumNodeGroup);

	if(*NumNodeGroup == 0)
		return 0;

	paillier_ciphertext_t*** cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(*NumNodeGroup)*FanOut);
		
	for(int i = 0 ; i < *NumNodeGroup * FanOut ; i++ ) {
		cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim+1);
		for(int j = 0 ; j < dim+1 ; j ++ ){
			cand[i][j] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
			mpz_init(cand[i][j]->c);
		}
	}

	int threadNum = thread_num;
	if(threadNum > FanOut)
	{
		threadNum = FanOut;
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

	std::thread *sNodeRetrievalforClassificationThread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		sNodeRetrievalforClassificationThread[i] = std::thread(sNodeRetrievalforClassificationinThread, std::ref(inputJ[i]), std::ref(inputNodeGroup[i]), std::ref(remained), std::ref(inputCnt[i]), std::ref(node), std::ref(data), std::ref(alpha), std::ref(cand), std::ref(*this), i);
	}

	for(int i = 0; i < threadNum; i++)
	{
		sNodeRetrievalforClassificationThread[i].join();
	}

	delete[] inputNodeGroup;
	delete[] inputJ;
	delete[] inputCnt;
	delete[] sNodeRetrievalforClassificationThread;

	return cand;
}


void protocol::SMSnMultithread2(int cnt, paillier_ciphertext_t** DIST_MINUS_MIN, paillier_ciphertext_t** DIST, paillier_ciphertext_t* MIN, paillier_ciphertext_t* C_RAND)
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
		SMSThread[i] = std::thread(SMSnInthread2, std::ref(inputCnt[i]), std::ref(DIST_MINUS_MIN), std::ref(DIST), std::ref(MIN), std::ref(C_RAND), std::ref(*this), i);
	}

	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i].join();
	}

	delete[] inputCnt;
	delete[] SMSThread;
}


void protocol::SMSnMultithread3(int cnt, int s, paillier_ciphertext_t **V, paillier_ciphertext_t ***V2, paillier_ciphertext_t **DIST, paillier_ciphertext_t ***cand, paillier_ciphertext_t *MAX, paillier_ciphertext_t ***Result)
{
	int dimension = dim;
	if( query == CLASSIFICATION) dimension = dim+1;

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
		SMSThread[i] = std::thread(SMSnInThread3_1, i, s, std::ref(inputCnt[i]), std::ref(V), std::ref(DIST), std::ref(cand), std::ref(MAX), std::ref(*this));
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
		thread_result[i] = new paillier_ciphertext_t*[dimension];
                for(int j = 0; j < dimension; j++)
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
		SMSThread[i] = std::thread(SMSnInThread3_2_v2, i, offset[i], record_size[i], dimension, std::ref(thread_result[i]), std::ref(V), std::ref(V2), std::ref(DIST), std::ref(cand),
			std::ref(MAX), std::ref(*this));
	}
	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i].join();
	}


	for(int i=0; i<numThread; i++)
	{
		for(int j=0; j<dimension; j++)
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




int** protocol::Classification_PB(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, paillier_ciphertext_t** Entire_set, boundary* node, int k, int NumData, int NumNode, int Entire_num)
{
	cout << "\n==============================  Classification_PB  START  ==============================\n" << endl;
	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}

	int n = 0, t = 0;
	int rand = 5;
	n = NumData;

	std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point end;

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

	paillier_ciphertext_t*** ciper_rabel_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	for(int i = 0 ; i < size; i++ ){
		ciper_Smin[i] = ciper_zero;
	}

	for(int i = 0 ; i < n ; i++ ){
		ciper_distance[i] 	= ciper_zero;
		ciper_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		ciper_V[i]	= ciper_zero;
		ciper_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
	}

	for(int i = 0 ; i < k ; i++ ){
		ciper_rabel_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
		for(int j = 0 ; j < dim+1 ; j++)
		{
			if (j == 0)
			{
				ciper_V2[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
				mpz_init(ciper_V2[i][j]->c);
			}
			ciper_rabel_result[i][j] = paillier_create_enc_zero();
		}
	}

	
	printf("\n=== Data SSED & SBD start ===\n");

	start = std::chrono::system_clock::now();
	
	
	int threadNum = thread_num;
	if(threadNum > n)
	{
		threadNum = n;
	}

	std::vector<int> *inputIdx = new std::vector<int>[threadNum];
	int count = 0;

	for(int i = 0; i < n; i++)
	{
		inputIdx[count++].push_back(i);
		count %= threadNum;
	}
	
	
	std::thread *SSED_SBD_thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		SSED_SBD_thread[i] = std::thread(SSED_SBD_inThread, i, std::ref(inputIdx[i]), std::ref(query), std::ref(data), dim,
			std::ref(ciper_distance), std::ref(ciper_SBD_distance), std::ref(*this));
	}

	for(int i = 0; i < threadNum; i++)
	{
		SSED_SBD_thread[i].join();
	}
	
	delete[] SSED_SBD_thread;
	delete[] inputIdx;
	
	end = std::chrono::system_clock::now();
	std::chrono::duration<float> SSED_SBD_time = end - start;
	data_SSED_SBD_time += SSED_SBD_time.count();
	//data_SSED_SBD_time += (float)(endTime-startTime)/(CLOCKS_PER_SEC);
	printf("data SSED & SBD time : %f\n", SSED_SBD_time.count());

	for(int s = 0; s < k; s++ )
	{
		start = std::chrono::system_clock::now();

		ciper_Smin = Smin_n_Multithread(ciper_SBD_distance, n);
		end = std::chrono::system_clock::now();
		std::chrono::duration<float> gaptime = end - start;
		gap = gaptime.count();
		sMINn_first_time += gap;
		printf("%dth sMINn time : %f\n", s+1, gap);
		
		for(int j = size; j > 0; j--)
		{
			t = (int)pow(2, j-1);
			ciper_binary = paillier_create_enc(t);
			ciper_binary = SM_p1(ciper_binary, ciper_Smin[size-j], 0);
			paillier_mul(pubkey, ciper_min, ciper_binary, ciper_min);
			paillier_freeciphertext(ciper_binary);

		}
		paillier_print("min dist : ", ciper_min);
		
		
		if(s != 0)
		{
			threadNum = thread_num;
			if(threadNum > n)
			{
				threadNum = n;
			}

			inputIdx = new std::vector<int>[threadNum];
			count = 0;
			for(int i = 0; i < n; i++)
			{
				inputIdx[count++].push_back(i);
				count %= threadNum;
			}

			std::thread *calcdistanceThread = new std::thread[threadNum];
			for(int i = 0; i < threadNum; i++)
			{
				calcdistanceThread[i] = std::thread(caldistanceInThread, i, std::ref(inputIdx[i]), std::ref(ciper_SBD_distance), std::ref(ciper_distance),
					std::ref(*this));
			}

			for(int i = 0; i < threadNum; i++)
			{
				calcdistanceThread[i].join();
			}

			delete[] inputIdx;
			delete[] calcdistanceThread;
		}


		
		threadNum = thread_num;
		if(threadNum > n)
		{
			threadNum = n;
		}

		inputIdx = new std::vector<int>[threadNum];
		std::thread *classipmThread = new std::thread[threadNum];
		count = 0;
		for(int i = 0; i < n; i++)
		{
			inputIdx[count++].push_back(i);
			count %= threadNum;
		}

		for(int i = 0; i < threadNum; i++)
		{
			classipmThread[i] = std::thread(classipmInThread1, i, std::ref(inputIdx[i]), std::ref(ciper_distance), std::ref(ciper_min), std::ref(ciper_rand),
				std::ref(ciper_mid), std::ref(*this));
		}

		for(int i = 0; i < threadNum; i++)
		{
			classipmThread[i].join();
		}
		
		delete[] classipmThread;
		delete[] inputIdx;
		

		ciper_V = SkNNm_sub(ciper_mid, n);



		
		for(int i = 0; i < n; i++ )
		{
			for (int j = 0 ; j < dim+1 ; j++)
			{
				ciper_V2[i][j] = SM_p1(ciper_V[i], data[i][j], 0);
				//paillier_print("ciper_V2 : ", ciper_V2[i]);
				paillier_mul(pubkey, ciper_rabel_result[s][j], ciper_V2[i][j], ciper_rabel_result[s][j]);
			}
		}

		printf("ciper_rabel_result : ");	
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper_rabel_result[s][dim]));
		
		
		start = std::chrono::system_clock::now();
		threadNum = thread_num;
		if(threadNum > n)
		{
			threadNum = n;
		}

		inputIdx = new std::vector<int>[threadNum];
		classipmThread = new std::thread[threadNum];
		count = 0;
		for(int i = 0; i < n; i++)
		{
			inputIdx[count++].push_back(i);
			count %= threadNum;
		}

		for(int i = 0; i < threadNum; i++)
		{
			classipmThread[i] = std::thread(SBOR_inThread, i, std::ref(inputIdx[i]), std::ref(ciper_V), std::ref(ciper_SBD_distance), std::ref(*this));
		}

		for(int i = 0; i < threadNum; i++)
		{
			classipmThread[i].join();
		}
		
		delete[] classipmThread;
		delete[] inputIdx;
		

		end = std::chrono::system_clock::now();
		std::chrono::duration<float> SBOR_time = end - start;
		data_SBOR_time += SBOR_time.count();

		ciper_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
	}



	paillier_ciphertext_t **ciper_rabel = new paillier_ciphertext_t*[k];
	
	paillier_plaintext_t* A = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	int** sknn = (int**)malloc(sizeof(int*)*k);
	for(int i = 0 ; i < k ; i++ ){
		sknn[i] = (int*)malloc(sizeof(int)*(dim+1));
		ciper_rabel[i] = ciper_rabel_result[i][dim];
		for(int j = 0 ; j < dim+1 ; j++ ){
			A = paillier_dec(0, pubkey, prvkey, ciper_rabel_result[i][j]);
			sknn[i][j] = mpz_get_ui(A->m);
			printf("%d ", sknn[i][j]);
		}
		printf("\n");
	}

	paillier_freeciphertext(temp_dist);	
	SCMC(Entire_set, ciper_rabel, Entire_num, k);
	return sknn;
}	


int** protocol::Classification_PGI(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, paillier_ciphertext_t** Entire_set, boundary* node, int k, int NumData, int NumNode, int Entire_num)
{
	printf("\n=== Classification_PGI start ===\n");

	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}


	int rand = 5;
	int cnt = 0;	 
	int NumNodeGroup = 0;
	bool verify_flag = false;

	std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point end;
	
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

	paillier_ciphertext_t **Rabel_Result = new paillier_ciphertext_t *[k];
	paillier_ciphertext_t ***Result = new paillier_ciphertext_t **[k];

	for(int i = 0; i < k; i++)
	{
		Rabel_Result[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));

		Result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
		for(int j = 0 ; j < dim + 1; j++)
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
				int NumThread = thread_num;
				if(NumThread > NumNode)
				{
					NumThread = NumNode;
				}
				std::vector<int> *NumNodeInput = new std::vector<int>[NumThread];
				int count = 0;
				for(int i = 0; i < NumNode; i++)
				{
					NumNodeInput[count++].push_back(i);
					count %= NumThread;
				}

				std::thread *GSROThread = new std::thread[NumThread];
				for(int i = 0; i < NumThread; i++)
				{
					GSROThread[i] = std::thread(DP_GSRO_Multithread, std::ref(query), std::ref(query), std::ref(NumNodeInput[i]),
						std::ref(node), std::ref(alpha), std::ref(*this), i);
				}

				for(int i = 0; i < NumThread; i++)
				{
					GSROThread[i].join();
				}

				delete[] NumNodeInput;
				delete[] GSROThread;
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
				NodeExpansionThread[i] = std::thread(NodeExpansion_Multithread, std::ref(query), std::ref(K_DIST), std::ref(NumNodeInput[i]),
					std::ref(node), std::ref(alpha), std::ref(*this), i);
			}

			for(int i = 0; i < NumThread; i++)
			{
				NodeExpansionThread[i].join();
			}

			delete[] NumNodeInput;
			delete[] NodeExpansionThread;
			
			if(Print)
			{
				for(int i = 0; i < NumNode; i++)
				{
					std::cout<< i+1 <<"Node Expansion alpha : ";
					gmp_printf("%Zd\n", paillier_dec(0, pubkey, prvkey, alpha[i]));
				}
			}
		}
		
		if(!verify_flag)	
		{
			//endTime = clock();
			//node_SRO_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			end = std::chrono::system_clock::now();
			std::chrono::duration<float> SRO_time = end - start;
			node_SRO_time = SRO_time.count();
			printf("GSRO time : %f\n", node_SRO_time);
		}
		else 
		{
			//endTime = clock();
			//node_expansion_time = (float)(endTime-startTime)/(CLOCKS_PER_SEC);
			end = std::chrono::system_clock::now();
			std::chrono::duration<float> SRO_time = end - start;
			node_expansion_time = SRO_time.count();
			printf("Node expansion time : %f\n", node_expansion_time);
		}
		
		cnt = 0;

		if(!verify_flag)
		{	
			cand = sNodeRetrievalforClassificationMultiThread(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			printf("Node Retrieval time : %f\n", data_extract_first_time);
		}else{
			start = std::chrono::system_clock::now();
			std::cout << "second GSRO_sNodeRetrievalforkNN start" << std::endl;
			temp_cand = sNodeRetrievalforClassificationMultiThread(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;

			end = std::chrono::system_clock::now();
			std::chrono::duration<float> SRO_time = end - start;
			data_extract_second_time = SRO_time.count();
			printf("Node Retrieval time : %f\n", data_extract_second_time);

			if(cnt == 0)
			{	
				break;
			}

			cand = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(cnt+k));
			for(int i = 0 ; i < cnt+k ; i++ )
			{
				cand[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim+1);
				for(int j = 0 ; j < dim+1 ; j ++ )
				{
					cand[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
					mpz_init(cand[i][j]->c);
				}
			}

			for(int i = 0 ; i < k ; i++ )
			{
				for(int j = 0 ; j < dim+1 ; j ++ )
				{
					cand[i][j] = Result[i][j];
				}
			}
			for(int i = 0 ; i < cnt ; i++ )
			{
				for(int j = 0 ; j < dim+1 ; j ++ )
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
		for(int i = 0 ; i < cnt ; i ++ )
		{
			DIST_MINUS_MIN[i]	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			DIST[i]				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			V[i]				= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			V2[i]				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
			
			mpz_init(DIST_MINUS_MIN[i]->c);
			mpz_init(DIST[i]->c);
			mpz_init(V[i]->c);
			
			for(int j = 0 ; j < dim+1 ; j++ ){
				V2[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
				mpz_init(V2[i][j]->c);
			}
			
		}
		
		SSEDMultiThread(DIST, cand, query, cnt);
	
		for(int s = 0 ; s < k ; s ++ )
		{
			std::cout << s+1 <<" th Classification Start !!!!!!!!!!!!"<<std::endl;
			
			idx = SMSn(DIST, cnt);			
			MIN = DIST[idx];
			SMSnMultithread2(cnt, DIST_MINUS_MIN, DIST, MIN, C_RAND);

			V = SkNNm_sub(DIST_MINUS_MIN, cnt);

			SMSnMultithread3(cnt, s, V, V2, DIST, cand, MAX, Result);
		}
		if(Print)
		{
			for(int i = 0 ; i < k ; i++ )
			{
				std::cout << "MIDDLE RESULT : ";
				for(int j = 0 ; j < dim+1 ; j++ )
				{
					gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, Result[i][j]));
				}
				std::cout << std::endl;
			}
		}
		if(!verify_flag)
		{	
			K_DIST = DP_SSED(query, Result[k-1], dim); 
		}

		if(!verify_flag)
			verify_flag = true;
		else
			break;
	}
	for(int i = 0 ; i < k ; i++ ){
		*Rabel_Result[i] = *Result[i][dim];
	}
	
	paillier_plaintext_t* A = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	int** sknn = (int**)malloc(sizeof(int*)*k);
	for(int i = 0 ; i < k ; i++ ){
		sknn[i] = (int*)malloc(sizeof(int)*(dim+1));
		for(int j = 0 ; j < dim+1 ; j++ ){
			A = paillier_dec(0, pubkey, prvkey, Result[i][j]);
			sknn[i][j] = mpz_get_ui(A->m);
			printf("%d ", sknn[i][j]);
		}
		printf("\n");
	}

	return sknn;
}

int** protocol::Classification_PAI(paillier_ciphertext_t*** data, paillier_ciphertext_t** query, paillier_ciphertext_t** Entire_set, boundary* node, int k, int NumData, int NumNode, int Entire_num)
{
	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}

}
