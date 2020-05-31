#include <iostream>
#include <fstream>
#include <gmp.h>
#include <string>
#include <sstream>
#include <thread>
#include <chrono>
#include <mutex>
#include <vector>
#include <iostream>
#include "protocol.h"
#include "parallel_protocol.h"

using namespace std;

// PARALLEL + BASIC (RANGE)
int** protocol::sRange_PB(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num)
{
	printf("\n=====now sRange_PB start=====\n");
	
	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}

	int i=0, j=0, m=0;	
	int rand = 5;
	int NumNodeGroup = 0;

	paillier_ciphertext_t*** cipher_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** cipher_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);

	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	
	cout << "\n=====================SBD QUERY=====================\n" << endl;
	startTime = std::chrono::system_clock::now(); // startTime check
	for(i=0; i<dim; i++) {
		cipher_qLL_bit[i] = SBD_for_SRO(q.LL[i], 0);		// query LL bound 변환
		cipher_qRR_bit[i] = SBD_for_SRO(q.RR[i], 1);		// query RR bound 변환
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["SBD_Q"] = time_variable.find("SBD_Q")->second + duration_sec.count();


	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumData);

	cout << "\n=====================PARALLEL SBD DATA and SRO=====================\n" << endl;

	startTime = std::chrono::system_clock::now(); // startTime check
	int threadNum = thread_num;
	if(threadNum > NumData)
	{
		threadNum = NumData;
	}

	std::vector<int> *inputIdx = new std::vector<int>[threadNum];
	int count = 0;

	for(int i = 0; i < NumData ; i++)
	{
		inputIdx[count++].push_back(i);
		count %= threadNum;
	}
	
	
	std::thread *SRO_SBD_Thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		SRO_SBD_Thread[i] = std::thread(SRO_SBD_inThread, std::ref(inputIdx[i]), std::ref(cipher_qLL_bit), std::ref(cipher_qRR_bit), std::ref(data), std::ref(alpha), std::ref(*this));
	}

	for(int i = 0; i < threadNum; i++)
	{
		SRO_SBD_Thread[i].join();
	}
	
	delete[] SRO_SBD_Thread;
	delete[] inputIdx;

	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["SBD_D + SRO_D"] = time_variable.find("SBD_D + SRO_D")->second + duration_sec.count();


	cout << "\n=======================EXTRACT RESULT=====================\n" << endl;


	paillier_ciphertext_t*** result = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData); //할당방법 변경해야하고

	//값 확인 && If αi = 1, E(t’i)를 result에 삽입
	for(i=0; i<NumData; i++){
		if(strcmp( paillier_plaintext_to_str( paillier_dec(0, pub, prv, alpha[i])), paillier_plaintext_to_str(plain_one)) == 0) 
		{
			result[(*result_num)]=data[i]; //result에 삽입
			(*result_num)++;
			printf("\n");
		}	
	}

	// user(Bob)에게 결과 전송을 위해 random 값 삽입
	for( i = 0 ; i<(*result_num) ; i++ )
	{
		for(j=0; j<dim; j++) 
		{	
				paillier_mul(pubkey,result[i][j],result[i][j],cipher_rand);
		}
	}
	return FsRange_Bob(result, rand, (*result_num), dim);
}



// PARALLEL + GARBLED + INDEX (RANGE)
int** protocol::sRange_PGI(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num)
{
	printf("\n===== Now sRange_PGI starts =====\n");
	int i=0, j=0, m=0;
	int rand = 5;
	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}
	int NumNodeGroup = 0;

	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	paillier_print("rand : ",cipher_rand);


	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수

	paillier_ciphertext_t*** cand ;


	cand = Parallel_GSRO_sNodeRetrievalforRange(data, q.LL, q.RR, node, NumData, NumNode, &cnt, &NumNodeGroup, 1);
	totalNumOfRetrievedNodes += NumNodeGroup;



	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);


	startTime = std::chrono::system_clock::now(); // startTime check
	Parallel_GSRO_inMultithread(cnt, cand, q, alpha, cipher_rand);
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["SRO_D"] = time_variable.find("SRO_D")->second + duration_sec.count();
	
	
	cout<<"sRangeG NumNodeGroup : " << NumNodeGroup<<endl;
	cout<<"sRangeG cnt : "<<cnt<<endl;
	paillier_ciphertext_t*** result;
	result = sRange_result(alpha, cand, cnt, NumNodeGroup, result_num);
	
	cout<<"sRange_result end"<<endl;

	return FsRange_Bob(result, rand, (*result_num), dim);
}


// PARALLEL + AS_CMP + INDEX (RANGE)
int** protocol::sRange_PAI(paillier_ciphertext_t*** data, boundary q, boundary* node, int NumData, int NumNode, int* result_num)
{

	printf("\n===== Now sRange_PAI starts =====\n");
	int i=0, j=0, m=0;
	int rand = 5;

	int NumNodeGroup = 0;

	paillier_ciphertext_t* cipher_rand = paillier_create_enc(rand);
	paillier_print("rand : ",cipher_rand);


	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수

	paillier_ciphertext_t*** cand ;

	cand = Parallel_GSRO_sNodeRetrievalforRange(data, q.LL, q.RR, node, NumData, NumNode, &cnt, &NumNodeGroup, 1);
	totalNumOfRetrievedNodes += NumNodeGroup;


	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*cnt);

	Parallel_GSRO_inMultithread(cnt, cand, q, alpha, cipher_rand);
	
	cout<<"sRangeG NumNodeGroup : " << NumNodeGroup<<endl;
	cout<<"sRangeG cnt : "<<cnt<<endl;
	paillier_ciphertext_t*** result;
	result = sRange_result(alpha, cand, cnt, NumNodeGroup, result_num);
	
	return FsRange_Bob(result, rand, (*result_num), dim);
}


paillier_ciphertext_t*** protocol::Parallel_GSRO_sNodeRetrievalforRange(paillier_ciphertext_t*** data, paillier_ciphertext_t** cipher_qLL, paillier_ciphertext_t** cipher_qRR, boundary* node, int NumData, int NumNode, int* cnt, int* NumNodeGroup, int type)
{
	printf("\n===== Now Parallel_GSRO_sNodeRetrievalfor %d starts =====\n", (int)query);
	printf("NumNode : %d thread_num : %d\n", NumNode, thread_num);
	int i=0, j=0, m=0;


	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	for(i=0; i<NumNode; i++){
		alpha[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(alpha[i]->c);
	}

	startTime = std::chrono::system_clock::now(); // startTime check
	// Check Overlap
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
		GSROThread[i] = std::thread(GSRO_Multithread, std::ref(cipher_qLL), std::ref(cipher_qRR), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), type);
	}

	for(int i = 0; i < NumThread; i++)
	{
		GSROThread[i].join();
	}

	delete[] NumNodeInput;
	delete[] GSROThread;

	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["nodeRetrievalSRO"] = time_variable.find("nodeRetrievalSRO")->second + duration_sec.count();


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

	startTime = std::chrono::system_clock::now(); // startTime check
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

	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["nodeRetrieval"] = time_variable.find("nodeRetrieval")->second + duration_sec.count();

	return cand;
}

paillier_ciphertext_t* protocol::GSRO_inThread(paillier_ciphertext_t** qLL, paillier_ciphertext_t** qRR,	paillier_ciphertext_t** nodeLL,	paillier_ciphertext_t** nodeRR, int idx)
{
	//cout << "\nGSRO start" <<endl;
	int i = 0;
	int random = 0;
	int cnt = 0 ;
	int * Flag_array = (int *)malloc(sizeof(int)*(dim*2));
	int * iR1_array = (int *)malloc(sizeof(int)*(dim*2));
	int * iR2_array = (int *)malloc(sizeof(int)*(dim*2));
	paillier_ciphertext_t** F1_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));
	paillier_ciphertext_t** F2_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));

	paillier_ciphertext_t** R1_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));
	paillier_ciphertext_t** R2_array = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim*2));
	
	for( i = 0 ; i < dim*2 ; i++ )
	{
		F1_array[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		F2_array[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		R1_array[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		R2_array[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(F1_array[i]->c);
		mpz_init(F2_array[i]->c);
		mpz_init(R1_array[i]->c);
		mpz_init(R2_array[i]->c);
		iR1_array[i] = 0;
		iR2_array[i] = 0;
		R1_array[i] = paillier_create_enc(iR1_array[i]);
		R2_array[i] = paillier_create_enc(iR2_array[i]);
		//	gmp_printf("%d : %Zd\t", i, F1_array[i]);		gmp_printf("%Zd\t", F2_array[i]);		gmp_printf("%Zd\t", R1_array[i]);		gmp_printf("%Zd\t", R2_array[i]);		cout << endl;
	}
	srand((unsigned int)time(NULL));
	for( i = 0 ; i < dim ; i ++ )
	{
		//random = rand()%2;
		//gmp_printf("%d : %Zd\t", random, paillier_dec(0, pub, prv, qLL[i]));		gmp_printf("%Zd\t", paillier_dec(0, pub, prv, nodeRR[i]));		cout << endl;
		if (random==0)
		{
			Flag_array[i] = 0;
			paillier_mul(pub, F1_array[i], qLL[cnt], R1_array[i]);
			paillier_mul(pub, F2_array[i], nodeRR[cnt++], R2_array[i]);
		}
		else
		{
			Flag_array[i] = 1;
			paillier_mul(pub, F1_array[i], nodeRR[cnt], R2_array[i]);
			paillier_mul(pub, F2_array[i], qLL[cnt++], R1_array[i]);
		}
	}

	for( i = dim, cnt = 0 ; i < dim*2 ; i++ )
	{
		//random = rand()%2;
		//gmp_printf("%d : %Zd\t", random, paillier_dec(0, pub, prv, qRR[cnt]));		gmp_printf("%Zd\t", paillier_dec(0, pub, prv, nodeLL[cnt]));		cout << endl;
		if( random == 0 )
		{
			Flag_array[i] = 0;
			paillier_mul(pub, F1_array[i], nodeLL[cnt], R1_array[i]);
			paillier_mul(pub, F2_array[i], qRR[cnt++], R2_array[i]);
		}
		else
		{
			Flag_array[i] = 1;
			paillier_mul(pub, F1_array[i], qRR[cnt], R2_array[i]);
			paillier_mul(pub, F2_array[i], nodeLL[cnt++], R1_array[i]);
		}
	}

	mtx.lock();
	F1_array = GSRO_sub(F1_array, F2_array, iR1_array, iR2_array);
	mtx.unlock();

	for( i = 0 ; i < dim*2 ; i ++){
		F1_array[0] = SM_p1(F1_array[0], F1_array[i]);
	}

	return F1_array[0];
}

paillier_ciphertext_t*** protocol::Parallel_GSRO_inMultithread(int cnt, paillier_ciphertext_t*** cand, boundary q, paillier_ciphertext_t** alpha, paillier_ciphertext_t* cipher_rand)
{	
	int NumThread = thread_num;
	if(NumThread > cnt)
	{
		NumThread = cnt;
	}
	std::vector<int> *NumNodeInput = new std::vector<int>[NumThread];
	int count = 0;
	for(int i = 0; i < cnt; i++)
	{
		NumNodeInput[count++].push_back(i);
		count %= NumThread;
	}

	std::thread *GSROThread = new std::thread[NumThread];
	for(int i = 0; i < NumThread; i++)
	{
		if (query == RANGE)
		{
			GSROThread[i] = std::thread(GSRO_dataquery_Multithread, std::ref(cand), std::ref(q), std::ref(NumNodeInput[i]), std::ref(alpha), std::ref(*this), i, std::ref(cipher_rand));
		}
		else if (query != RANGE)
		{
			cout << "THIE QUERY IS NOT RANGE" << endl;
			//GSROThread[i] = std::thread(GSRO_Multithread, std::ref(q), std::ref(q), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), i);
		}		
	}

	for(int i = 0; i < NumThread; i++)
	{
		GSROThread[i].join();
	}

	delete[] NumNodeInput;
	delete[] GSROThread;
}
