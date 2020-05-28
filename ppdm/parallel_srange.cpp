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
	
	time_t startTime = 0;
	time_t endTime = 0;
	float gap = 0.0;
	int rand = 5;
	int NumNodeGroup = 0;

	paillier_ciphertext_t*** ciper_qLL_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);
	paillier_ciphertext_t*** ciper_qRR_bit = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*dim);

	paillier_ciphertext_t* ciper_rand = paillier_create_enc(rand);
	
	cout << "\n=====================SBD QUERY=====================\n" << endl;
	for(i=0; i<dim; i++) {
		ciper_qLL_bit[i] = SBD_for_SRO(q.LL[i], 0);		// query LL bound 변환
		ciper_qRR_bit[i] = SBD_for_SRO(q.RR[i], 1);		// query RR bound 변환
	}


	paillier_ciphertext_t** alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*NumData);

	cout << "\n=====================PARALLEL SBD DATA and SRO=====================\n" << endl;
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
		SRO_SBD_Thread[i] = std::thread(SRO_SBD_inThread, std::ref(inputIdx[i]), std::ref(ciper_qLL_bit), std::ref(ciper_qRR_bit), std::ref(data), std::ref(alpha), std::ref(*this));
	}

	for(int i = 0; i < threadNum; i++)
	{
		SRO_SBD_Thread[i].join();
	}
	
	delete[] SRO_SBD_Thread;
	delete[] inputIdx;

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
				paillier_mul(pubkey,result[i][j],result[i][j],ciper_rand);
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

	Parallel_GSRO_inMultithread(cnt, cand, q, alpha, cipher_rand);
	
	
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
		GSROThread[i] = std::thread(DP_GSRO_Multithread, std::ref(cipher_qLL), std::ref(cipher_qRR), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), type);
	}

	for(int i = 0; i < NumThread; i++)
	{
		GSROThread[i].join();
	}

	delete[] NumNodeInput;
	delete[] GSROThread;


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

paillier_ciphertext_t* protocol::DP_GSRO_inThread(paillier_ciphertext_t** qLL, paillier_ciphertext_t** qRR,	paillier_ciphertext_t** nodeLL,	paillier_ciphertext_t** nodeRR, int idx)
{
	int *Flag_array = new int[dim * 2];
	int *R1 = new int[dim * 2];
	int *R2 = new int[dim * 2];
	for(int i = 0; i < dim * 2; i++)
	{
		Flag_array[i] = 0;
		R1[i] = 5;
		R2[i] = 5;
	}
	
	paillier_ciphertext_t **F1_array = new paillier_ciphertext_t *[dim * 2];
	paillier_ciphertext_t **F2_array;

	paillier_plaintext_t *shift = paillier_plaintext_from_ui(0);
	paillier_plaintext_t *p1_array = paillier_plaintext_from_ui(0);
	paillier_plaintext_t *p2_array = paillier_plaintext_from_ui(0);
	paillier_plaintext_t **R1_array = new paillier_plaintext_t *[dim * 2];
	paillier_plaintext_t **R2_array = new paillier_plaintext_t *[dim * 2];
	paillier_plaintext_t *tmp = paillier_plaintext_from_ui(0);
	
	paillier_ciphertext_t *c1_array;
	paillier_ciphertext_t *c2_array;
	paillier_ciphertext_t *ciper_tmp = new paillier_ciphertext_t;
	mpz_init(ciper_tmp->c);

	paillier_ciphertext_t **temp_qLL = new paillier_ciphertext_t *[dim];
	paillier_ciphertext_t **temp_qRR = new paillier_ciphertext_t *[dim];
	paillier_ciphertext_t **temp_nodeLL = new paillier_ciphertext_t *[dim];
	paillier_ciphertext_t **temp_nodeRR = new paillier_ciphertext_t *[dim];

	for(int i = 0; i <dim; i++)
	{
		temp_qLL[i] = new paillier_ciphertext_t;
		temp_qRR[i] = new paillier_ciphertext_t;
		temp_nodeLL[i] = new paillier_ciphertext_t;
		temp_nodeRR[i] = new paillier_ciphertext_t;
		mpz_init(temp_qLL[i]->c);
		mpz_init(temp_qRR[i]->c);
		mpz_init(temp_nodeLL[i]->c);
		mpz_init(temp_nodeRR[i]->c);
	}

	// after shifting and add random value under Array P1
	for(int i = 0; i < dim * 2; i++)
	{
		mpz_ui_pow_ui(shift->m, 2, DP * i);
		R1_array[i] = paillier_plaintext_from_ui(R1[i]);
		mpz_mul(tmp->m, R1_array[i]->m, shift->m);
		mpz_add(p1_array->m, p1_array->m, tmp->m);
	}

	// after shifting and add random value under Array P2
	for(int i = 0; i < dim * 2; i++)
	{
		mpz_ui_pow_ui(shift->m, 2, DP * i);
		R2_array[i] = paillier_plaintext_from_ui(R2[i]);
		mpz_mul(tmp->m, R2_array[i]->m, shift->m);
		mpz_add(p2_array->m, p2_array->m, tmp->m);
	}
	
	// enc Array P1, enc Array P2
	c1_array = paillier_enc(0, pubkey, p1_array, paillier_get_rand_devurandom);
	c2_array = paillier_enc(0, pubkey, p2_array, paillier_get_rand_devurandom);

	int random = 0 ;
	int cnt = 0 ;	

	// 2qLL 2nodeLL
	for(int i = 0; i < dim; i++)
	{
		mpz_ui_pow_ui(shift->m,2,1);
		paillier_exp(pubkey, temp_qLL[i], qLL[i], shift);
		paillier_exp(pubkey, temp_nodeLL[i], nodeLL[i], shift);
	}
	// 2qRR+1 2nodeRR+1
	for(int i = 0; i < dim; i++){
		paillier_exp(pubkey, temp_nodeRR[i], nodeRR[i], shift); 
		paillier_mul(pubkey, temp_nodeRR[i], temp_nodeRR[i], ciper_one); 
		
		paillier_exp(pubkey, temp_qRR[i], qRR[i], shift);
		paillier_mul(pubkey, temp_qRR[i], temp_qRR[i], ciper_one); 
	}

	// add to Array QueryLL, nodeRR value
	for(int i = 0; i < dim; i++)
	{
		random = rand() % 2;
		mpz_ui_pow_ui(shift->m, 2, DP * i);
		if (random == 0)
		{
			Flag_array[i] = 0;
			paillier_exp(pubkey, ciper_tmp, temp_qLL[cnt], shift); //multi
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
			paillier_exp(pubkey, ciper_tmp, temp_nodeRR[cnt++], shift); //multi
			paillier_mul(pubkey, c2_array, c2_array, ciper_tmp); //E_add			
		}
		else
		{
			Flag_array[i] = 1;
			paillier_exp(pubkey, ciper_tmp, temp_nodeRR[cnt], shift); //multi
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
			paillier_exp(pubkey, ciper_tmp, temp_qLL[cnt++], shift); //multi
			paillier_mul(pubkey, c2_array, c2_array, ciper_tmp); //E_add
		}
	}
	// add to Array QueryRR, nodeLL value	
	for(int i = dim, cnt = 0; i < dim * 2; i++)
	{
		random = rand() % 2;
		mpz_ui_pow_ui(shift->m, 2, DP * i);
		if(random == 0)
		{
			Flag_array[i] = 0;
			paillier_exp(pubkey, ciper_tmp, temp_nodeLL[cnt], shift); //multi
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
			paillier_exp(pubkey, ciper_tmp, temp_qRR[cnt++], shift); //multi
			paillier_mul(pubkey, c2_array, c2_array, ciper_tmp); //E_add			
		}
		else
		{
			Flag_array[i] = 1;
			paillier_exp(pubkey, ciper_tmp, temp_qRR[cnt], shift); //multi
			paillier_mul(pubkey, c1_array, c1_array, ciper_tmp); //E_add
			paillier_exp(pubkey, ciper_tmp, temp_nodeLL[cnt++], shift); //multi
			paillier_mul(pubkey, c2_array, c2_array, ciper_tmp); //E_add			
		}
	}

	mtx.lock();
	F2_array = unDP_GSRO(c1_array, c2_array, R1, R2);
	mtx.unlock();
	
	for(int i = 0; i < dim * 2; i++)
	{
		F1_array[i] = paillier_create_enc(Flag_array[i]);
		F1_array[i] = SBXOR(F1_array[i], F2_array[(dim * 2) - i - 1]);
		
	}

	for(int i = 0; i < dim * 2; i++)
	{
		F1_array[0] = SM_p1(F1_array[0], F1_array[i]);
	}

	delete[] Flag_array;
	delete[] R1;
	delete[] R2;
	delete shift;
	delete p1_array;
	delete p2_array;
	delete tmp;
	delete ciper_tmp;
	for(int i = 0; i < dim * 2; i++)
	{
		delete[] R1_array[i];
		delete[] R2_array[i];
		delete[] F2_array[i];
	}
	delete[] R1_array;
	delete[] R2_array;
	delete[] F2_array;

	for(int i = 0; i < dim; i++)
	{
		delete[] temp_qLL[i];
		delete[] temp_qRR[i];
		delete[] temp_nodeLL[i];
		delete[] temp_nodeRR[i];
	}
	delete[] temp_qLL;
	delete[] temp_qRR;
	delete[] temp_nodeLL;
	delete[] temp_nodeRR;

	//mtx.unlock();
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
			GSROThread[i] = std::thread(DP_GSRO_dataquery_Multithread, std::ref(cand), std::ref(q), std::ref(NumNodeInput[i]), std::ref(alpha), std::ref(*this), i, std::ref(cipher_rand));
		}
		else if (query != RANGE)
		{
			cout << "THIE QUERY IS NOT RANGE" << endl;
			//GSROThread[i] = std::thread(DP_GSRO_Multithread, std::ref(q), std::ref(q), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), i);
		}		
	}

	for(int i = 0; i < NumThread; i++)
	{
		GSROThread[i].join();
	}

	delete[] NumNodeInput;
	delete[] GSROThread;
}
