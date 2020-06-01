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
	printf("\n===== Now sNodeRetrievalforClassificationMultiThread starts =====\n");
	
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
		SMSThread[i] = std::thread(SMSnInThread3_2_v2, i, offset[i], record_size[i], dimension, std::ref(thread_result[i]), std::ref(V), std::ref(V2), std::ref(DIST), std::ref(cand), std::ref(MAX), std::ref(*this));
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

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* cipher_binary;
	paillier_ciphertext_t* cipher_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_dist	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_rand	= paillier_create_enc(rand);
	paillier_ciphertext_t* temp_dist = paillier_create_enc_zero();

	paillier_ciphertext_t** cipher_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);
	paillier_ciphertext_t*** cipher_SBD_distance = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	paillier_ciphertext_t** cipher_mid = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);
	paillier_ciphertext_t** cipher_Smin = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);
	paillier_ciphertext_t** cipher_V=(paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);
	paillier_ciphertext_t*** cipher_V2 = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);

	paillier_ciphertext_t*** cipher_rabel_result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	for(int i = 0 ; i < size; i++ ){
		cipher_Smin[i] = cipher_zero;
	}

	for(int i = 0 ; i < NumData ; i++ ){
		cipher_distance[i] 	= cipher_zero;
		cipher_SBD_distance[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*size);
		cipher_V[i]	= cipher_zero;
		cipher_V2[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
	}

	for(int i = 0 ; i < k ; i++ ){
		cipher_rabel_result[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim+1));
		for(int j = 0 ; j < dim+1 ; j++)
		{
			if (j == 0)
			{
				cipher_V2[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
				mpz_init(cipher_V2[i][j]->c);
			}
			cipher_rabel_result[i][j] = paillier_create_enc_zero();
		}
	}

	
	printf("\n=== Data SSED & SBD start ===\n");



	startTime = std::chrono::system_clock::now(); // startTime check			
	
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
		SSED_SBD_thread[i] = std::thread(SSED_SBD_inThread, i, std::ref(inputIdx[i]), std::ref(query), std::ref(data), dim,
			std::ref(cipher_distance), std::ref(cipher_SBD_distance), std::ref(*this));
	}

	for(int i = 0; i < threadNum; i++)
	{
		SSED_SBD_thread[i].join();
	}
	
	delete[] SSED_SBD_thread;
	delete[] inputIdx;
	
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["SSED&SBD"] = time_variable.find("SSED&SBD")->second + duration_sec.count();


	for(int s = 0; s < k; s++ )
	{
		startTime = std::chrono::system_clock::now(); // startTime check			

		cipher_Smin = Smin_n_Multithread(cipher_SBD_distance, NumData);

		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["SMIN"] = time_variable.find("SMIN")->second + duration_sec.count();
		
		for(int j = size; j > 0; j--)
		{
			t = (int)pow(2, j-1);
			cipher_binary = paillier_create_enc(t);
			cipher_binary = SM_p1(cipher_binary, cipher_Smin[size-j], 0);
			paillier_mul(pubkey, cipher_min, cipher_binary, cipher_min);
			paillier_freeciphertext(cipher_binary);

		}
		paillier_print("min dist : ", cipher_min);
		

		startTime = std::chrono::system_clock::now(); // startTime check			
		if(s != 0)
		{
			threadNum = thread_num;
			if(threadNum > NumData)
			{
				threadNum = NumData;
			}

			inputIdx = new std::vector<int>[threadNum];
			count = 0;
			for(int i = 0; i < NumData; i++)
			{
				inputIdx[count++].push_back(i);
				count %= threadNum;
			}

			std::thread *calcdistanceThread = new std::thread[threadNum];
			for(int i = 0; i < threadNum; i++)
			{
				calcdistanceThread[i] = std::thread(caldistanceInThread, i, std::ref(inputIdx[i]), std::ref(cipher_SBD_distance), std::ref(cipher_distance),
					std::ref(*this));
			}

			for(int i = 0; i < threadNum; i++)
			{
				calcdistanceThread[i].join();
			}

			delete[] inputIdx;
			delete[] calcdistanceThread;
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["Bi to De"] = time_variable.find("Bi to De")->second + duration_sec.count();


		startTime = std::chrono::system_clock::now(); // startTime check					
		threadNum = thread_num;
		if(threadNum > NumData)
		{
			threadNum = NumData;
		}

		inputIdx = new std::vector<int>[threadNum];
		std::thread *classipmThread = new std::thread[threadNum];
		count = 0;
		for(int i = 0; i < NumData; i++)
		{
			inputIdx[count++].push_back(i);
			count %= threadNum;
		}

		for(int i = 0; i < threadNum; i++)
		{
			classipmThread[i] = std::thread(classipmInThread1, i, std::ref(inputIdx[i]), std::ref(cipher_distance), std::ref(cipher_min), std::ref(cipher_rand),
				std::ref(cipher_mid), std::ref(*this));
		}

		for(int i = 0; i < threadNum; i++)
		{
			classipmThread[i].join();
		}
		
		delete[] classipmThread;
		delete[] inputIdx;

		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["subtract"] = time_variable.find("subtract")->second + duration_sec.count();


		cipher_V = SkNNm_sub(cipher_mid, NumData);



		startTime = std::chrono::system_clock::now(); // startTime check					
		ExtractKNNInCLASSIFICATION_PB_inMultithread(data, cipher_V, cipher_rabel_result, NumData, s);
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["kNN"] = time_variable.find("kNN")->second + duration_sec.count();


		printf("cipher_rabel_result : ");	
		gmp_printf("%Zd \n", paillier_dec(0, pubkey, prvkey, cipher_rabel_result[s][dim]));
		
		
		startTime = std::chrono::system_clock::now(); // startTime check					
		threadNum = thread_num;
		if(threadNum > NumData)
		{
			threadNum = NumData;
		}

		inputIdx = new std::vector<int>[threadNum];
		classipmThread = new std::thread[threadNum];
		count = 0;
		for(int i = 0; i < NumData; i++)
		{
			inputIdx[count++].push_back(i);
			count %= threadNum;
		}

		for(int i = 0; i < threadNum; i++)
		{
			classipmThread[i] = std::thread(SBOR_inThread, i, std::ref(inputIdx[i]), std::ref(cipher_V), std::ref(cipher_SBD_distance), std::ref(*this));
		}

		for(int i = 0; i < threadNum; i++)
		{
			classipmThread[i].join();
		}
		
		delete[] classipmThread;
		delete[] inputIdx;

		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["update"] = time_variable.find("update")->second + duration_sec.count();

		cipher_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
	}



	paillier_ciphertext_t **cipher_rabel = new paillier_ciphertext_t*[k];
	
	paillier_plaintext_t* A = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	int** sknn = (int**)malloc(sizeof(int*)*k);
	for(int i = 0 ; i < k ; i++ ){
		sknn[i] = (int*)malloc(sizeof(int)*(dim+1));
		cipher_rabel[i] = cipher_rabel_result[i][dim];
		for(int j = 0 ; j < dim+1 ; j++ ){
			A = paillier_dec(0, pubkey, prvkey, cipher_rabel_result[i][j]);
			sknn[i][j] = mpz_get_ui(A->m);
			printf("%d ", sknn[i][j]);
		}
		printf("\n");
	}

	paillier_freeciphertext(temp_dist);	


	startTime = std::chrono::system_clock::now(); // startTime check					
	SCMC(Entire_set, cipher_rabel, Entire_num, k);
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["classify"] = time_variable.find("classify")->second + duration_sec.count();


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

		startTime = std::chrono::system_clock::now(); // startTime check					

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
					GSROThread[i] = std::thread(GSRO_Multithread, std::ref(query), std::ref(query), std::ref(NumNodeInput[i]),
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
		{	
			startTime = std::chrono::system_clock::now(); // startTime check						
			cand = sNodeRetrievalforClassificationMultiThread(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["nodeRetrieval"] = time_variable.find("nodeRetrieval")->second + duration_sec.count();
		}else{
			startTime = std::chrono::system_clock::now(); // startTime check						
			std::cout << "second GSRO_sNodeRetrievalforkNN start" << std::endl;
			temp_cand = sNodeRetrievalforClassificationMultiThread(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["exnodeRetrieval"] = time_variable.find("exnodeRetrieval")->second + duration_sec.count();

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
	
		startTime = std::chrono::system_clock::now(); // startTime check						
		SSEDMultiThread(DIST, cand, query, cnt);
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["SSED"] = time_variable.find("SSED")->second + duration_sec.count();

		for(int s = 0 ; s < k ; s ++ )
		{
			std::cout << s+1 <<" th Classification Start !!!!!!!!!!!!"<<std::endl;
			startTime = std::chrono::system_clock::now(); // startTime check						
			idx = SMSn(DIST, cnt);
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["SMIN"] = time_variable.find("SMIN")->second + duration_sec.count();
			MIN = DIST[idx];

			startTime = std::chrono::system_clock::now(); // startTime check						
			SMSnMultithread2(cnt, DIST_MINUS_MIN, DIST, MIN, C_RAND);
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["subtract"] = time_variable.find("subtract")->second + duration_sec.count();

			V = SkNNm_sub(DIST_MINUS_MIN, cnt);

			startTime = std::chrono::system_clock::now(); // startTime check						
			SMSnMultithread3(cnt, s, V, V2, DIST, cand, MAX, Result);
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["kNN"] = time_variable.find("kNN")->second + duration_sec.count();
		}

		if(!verify_flag)
		{	
			K_DIST = SSEDm(query, Result[k-1], dim); 
		}

		if(!verify_flag)
			verify_flag = true;
		else
			break;
	}
	for(int i = 0 ; i < k ; i++ ){
		*Rabel_Result[i] = *Result[i][dim];
	}
	

	paillier_ciphertext_t**	cipher_rabel = new paillier_ciphertext_t*[k];
	paillier_plaintext_t* A = (paillier_plaintext_t*)malloc(sizeof(paillier_plaintext_t));
	int** sknn = (int**)malloc(sizeof(int*)*k);
	for(int i = 0 ; i < k ; i++ ){
		cipher_rabel[i] = new paillier_ciphertext_t;
		mpz_init(cipher_rabel[i]->c);
		sknn[i] = (int*)malloc(sizeof(int)*(dim+1));
		for(int j = 0 ; j < dim+1 ; j++ ){
			A = paillier_dec(0, pubkey, prvkey, Result[i][j]);
			sknn[i][j] = mpz_get_ui(A->m);
			printf("%d ", sknn[i][j]);
			if ( j == dim)
			{
				cipher_rabel[i] = Result[i][j];
			}
		}
		printf("\n");
	}

	startTime = std::chrono::system_clock::now(); // startTime check					
	SCMC(Entire_set, cipher_rabel, Entire_num, k);
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["classify"] = time_variable.find("classify")->second + duration_sec.count();



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


void protocol::ExtractKNNInCLASSIFICATION_PB_inMultithread(paillier_ciphertext_t*** data, paillier_ciphertext_t** cipher_V, paillier_ciphertext_t*** cipher_result, int NumData, int serialKNN)
{
	cout << "ExtractKNNInCLASSIFICATION_PB_inMultithread" << endl;
	int threadNum = thread_num;
	if(threadNum > NumData)
	{
		threadNum = NumData;
	}

	int dimension = dim + 1;

	paillier_ciphertext_t*** cipher_thread_result = new paillier_ciphertext_t**[threadNum];
	for( int i = 0 ; i < threadNum ; i++ )
	{
		cipher_thread_result[i] = new paillier_ciphertext_t*[dimension];
		for( int j = 0 ; j < dimension ; j++ )
		{
			cipher_thread_result[i][j] = paillier_create_enc(0);
		}
	}

	std::vector<int> *inputIdx = new std::vector<int>[threadNum];
	int count = 0;

	for(int i = 0; i < NumData; i++)
	{
		inputIdx[count++].push_back(i);
		count %= threadNum;
	}
	std::thread *ExtractKNNInCLASSIFICATION_PB_Thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		ExtractKNNInCLASSIFICATION_PB_Thread[i] = std::thread(ExtractKNNInCLASSIFICATION_PB_inThread, std::ref(inputIdx[i]), std::ref(data), std::ref(cipher_V), std::ref(cipher_thread_result[i]), std::ref(*this), dimension);
	}
	for(int i = 0; i < threadNum; i++)
	{
		ExtractKNNInCLASSIFICATION_PB_Thread[i].join();
	}


	for( int i = 0 ; i < threadNum ; i++ )
	{
		for( int j = 0 ; j < dimension ; j++ )
		{
			paillier_mul(pubkey, cipher_result[serialKNN][j], cipher_result[serialKNN][j], cipher_thread_result[i][j]);
		}
	}

	delete[] cipher_thread_result;
	delete[] ExtractKNNInCLASSIFICATION_PB_Thread;
	delete[] inputIdx;
}