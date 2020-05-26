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
			//std::cout << "inputCnt = " << *cnt << std::endl;
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
	/*
	SMSThread = new std::thread[numThread];
	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i] = std::thread(SMSnInThread3_2, i, s, std::ref(inputCnt[i]), std::ref(V), std::ref(V2), std::ref(DIST), std::ref(cand), std::ref(MAX), std::ref(Result),
		std::ref(*this));
	}

	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i].join();
	}
	*/
	delete[] inputCnt;
	/*
	if(numThread > dim + 1)
	{
		numThread = dim + 1;
	}
	*/
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
		thread_result[i] = new paillier_ciphertext_t*[dim+1];
                for(int j = 0; j < dim+1; j++)
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
	/*
	count = 0;
	for(int i = 0; i < dim + 1; i++)
	{
		inputCnt[count++].push_back(i);
		count %= numThread;
	}
	*/
	//for(int i = 0; i < cnt; i++)
	//{
	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i] = std::thread(SMSnInThread3_2_v2, i, offset[i], record_size[i], dim, std::ref(thread_result[i]), std::ref(V), std::ref(V2), std::ref(DIST), std::ref(cand),
			std::ref(MAX), std::ref(*this));
	}
	for(int i = 0; i < numThread; i++)
	{
		SMSThread[i].join();
	}

	//}
	for(int i=0; i<numThread; i++)
	{
		for(int j=0; j<dim+1; j++)
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