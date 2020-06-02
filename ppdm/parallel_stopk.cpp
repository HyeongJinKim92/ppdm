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
	cout << "\n ======================= STopk_PB start ===============================\n " << endl;

	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}


	int i = 0;
	int s = 0;
	int t = 0;
	int j = 0;
	int rand = 5;
	
	paillier_ciphertext_t * hint = paillier_create_enc_zero();

	paillier_plaintext_t * pt = paillier_plaintext_from_ui(0);

	paillier_ciphertext_t* cipher_binary;
	paillier_ciphertext_t* cipher_min	= paillier_create_enc_zero();
	paillier_ciphertext_t* cipher_dist	= paillier_create_enc_zero();
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


	printf("\n=== PARALLEL Data compute & SBD start ===\n");
	startTime = std::chrono::system_clock::now(); // startTime check	
	TOPK_CS_SBD_inMultithread(data, q, cipher_SBD_distance, cipher_distance);
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["CS&SBD"] = time_variable.find("CS&SBD")->second + duration_sec.count();





	for( s = 0 ; s < k ; s++ )
	{
		printf("\n%dth sMAX start \n", s+1);
		startTime = std::chrono::system_clock::now(); // startTime check	
		cipher_Smin = Smax_n_Multithread(cipher_SBD_distance, NumData); // bit로 표현된 암호화 min 거리 추출 parallel
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["SMAX"] = time_variable.find("SMAX")->second + duration_sec.count();


		startTime = std::chrono::system_clock::now(); // startTime check		
		// bit로 표현된 암호화 max 거리를 암호화 정수로 변환
		for( j = size ; j > 0 ; j-- ){
			t = (int)pow(2, j-1);
			cipher_binary = paillier_create_enc(t);
			cipher_binary = SM_p1(cipher_binary, cipher_Smin[size-j]);
			paillier_mul(pubkey, cipher_min, cipher_binary, cipher_min);
			paillier_freeciphertext(cipher_binary);
		}
		paillier_print("max dist : ", cipher_min);
		
		printf("\n== recalculate query<->data distances ===\n");
		// 이전 iteration에서 min값으로 선택된 데이터의 거리가 secure하게 MAX로 변환되었기 때문에, 
		// 질의-데이터 간 거리 계산을 모두 다시 수행
		if( s != 0 ){
			recalculate_SCORE_inMultithread(cipher_SBD_distance, cipher_distance, NumData);
		}
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["Bi to De"] = time_variable.find("Bi to De")->second + duration_sec.count();

		// 질의-데이터 거리와 min 거리와의 차를 구함 (min 데이터의 경우에만 0으로 만들기 위함)
		startTime = std::chrono::system_clock::now(); // startTime check		
		DuringProcessedInTOPK_PB_inMultithread(cipher_distance, cipher_mid, cipher_min, cipher_rand, NumData);
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["subtract"] = time_variable.find("subtract")->second + duration_sec.count();

		cipher_V = Topk_sub(cipher_mid, NumData);

		// min 데이터 추출
		startTime = std::chrono::system_clock::now(); // startTime check		
		ExtractTOPKInTOPK_PB_inMultithread(data, cipher_V, cipher_result, NumData, s);
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["Topk"] = time_variable.find("Topk")->second + duration_sec.count();

		startTime = std::chrono::system_clock::now(); // startTime check		
		UPDATE_SBD_SCORE_InTOPK_PB_inMultithread(cipher_SBD_distance, cipher_V, NumData);
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["update"] = time_variable.find("update")->second + duration_sec.count();

		cipher_min = paillier_enc(0, pubkey, pt, paillier_get_rand_devurandom);	
	}
	// user(Bob)에게 결과 전송을 위해 random 값 삽입
	for(i=0;i<k;i++){
		gmp_printf("%d : %Zd : \n", i, paillier_dec(0, pubkey, prvkey, computeScore(q, cipher_result[i])));
		for(j=0; j<dim; j++){
			paillier_mul(pubkey, cipher_result[i][j], cipher_result[i][j], cipher_rand);
		}
	}
	return SkNNm_Bob2(cipher_result, rand, k, dim);
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


	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}


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
	
	startTime = std::chrono::system_clock::now(); // startTime check	
	hint = q[dim];

	for( i = 0 ; i < dim ; i++ )
	{
		paillier_mul(pubkey, coeff[i], q[i], hint);
	}
	for( i = 0 ; i < dim ; i++ )
	{
		psi[i] = GSCMP(hint, coeff[i]);
	}
	for( i = 0 ; i < dim ; i++ )
	{	
		Q[i] = SM_p1(max_val[i], psi[i]);
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["Query"] = time_variable.find("Query")->second + duration_sec.count();


	startTime = std::chrono::system_clock::now(); // startTime check	

	Parallel_GSRO_Topk(Q, Q, alpha, node, NumNode, &cnt, &NumNodeGroup, true);
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["nodeRetrievalSRO"] = time_variable.find("nodeRetrievalSRO")->second + duration_sec.count();

	while(1)
	{
		cnt = 0;
		if(!verify_flag)
		{	
			startTime = std::chrono::system_clock::now(); // startTime check	
			cand = PARALLEL_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			cout << "NumNodeGroup : " << NumNodeGroup << endl;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["nodeRetrieval"] = time_variable.find("nodeRetrieval")->second + duration_sec.count();
		}
		else
		{
			startTime = std::chrono::system_clock::now(); // startTime check	
			temp_cand = PARALLEL_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			cout << "Expand NumNodeGroup : " << NumNodeGroup << endl;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["exnodeRetrieval"] = time_variable.find("exnodeRetrieval")->second + duration_sec.count();

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

		startTime = std::chrono::system_clock::now(); // startTime check	
		ComputeScoreinMultithread(cand, q, SCORE, &cnt, true);
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["CS"] = time_variable.find("CS")->second + duration_sec.count();

		for( s = 0 ; s < k ; s++ )
		{
			cout << s+1<<"dth MAXn_Topk start"<<endl;
			startTime = std::chrono::system_clock::now(); // startTime check	
			MAX_idx = MAXn(SCORE, cnt);
			MAX = SCORE[MAX_idx];
			if(!verify_flag) {
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["SMAX"] = time_variable.find("SMAX")->second + duration_sec.count();
			}
			else {
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["exSMAX"] = time_variable.find("exSMAX")->second + duration_sec.count();
			}
			startTime = std::chrono::system_clock::now(); // startTime check	
			MAXnMultithread2(cnt, SCORE_MINUS_MAX, SCORE, MAX, C_RAND);
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["subtract"] = time_variable.find("subtract")->second + duration_sec.count();

			V = Topk_sub(SCORE_MINUS_MAX, cnt);

			startTime = std::chrono::system_clock::now(); // startTime check	
			MAXnMultithread3(cnt, s, V, V2, SCORE, cand, MIN, RESULT);
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["Topk&update"] = time_variable.find("Topk&update")->second + duration_sec.count();
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
			if(cnt < k)	 // 노드의 fanout이 k보다 적은 경우를 처리함
				kth_score = computeScore(q, RESULT[cnt-1]);
			else	// k개가 찾아졌다면, k번째 데이터의 score를 계산
				kth_score = computeScore(q, RESULT[k-1]);

			startTime = std::chrono::system_clock::now(); // startTime check	
			int threadNum = thread_num;
			if(threadNum > NumNode)
			{
				threadNum = NumNode;
			}
			std::vector<int> *inputIdx = new std::vector<int>[threadNum];
			int count = 0;

			for(int i = 0; i < NumNode ; i++)
			{
				inputIdx[count++].push_back(i);
				count %= threadNum;
			}					
			std::thread *GSCMPforTopk_Thread = new std::thread[threadNum];
			for(int i = 0; i < threadNum; i++)
			{
				GSCMPforTopk_Thread[i] = std::thread(GSCMPforTopk_inThread, std::ref(inputIdx[i]), std::ref(node), std::ref(psi), std::ref(kth_score), std::ref(q), std::ref(alpha), std::ref(*this));
			}
			for(int i = 0; i < threadNum; i++)
			{
				GSCMPforTopk_Thread[i].join();
			}
			delete[] GSCMPforTopk_Thread;
			delete[] inputIdx;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["excheck"] = time_variable.find("excheck")->second + duration_sec.count();
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
		gmp_printf(" %Zd \n", paillier_dec(0, pubkey, prvkey, computeScore(q, RESULT[i])));
	
		for( j = 0 ; j < dim ; j++ )
		{
			paillier_mul(pubkey, RESULT[i][j], RESULT[i][j], C_RAND);
		}
	}
	cout << "End Line" <<endl;
	return SkNNm_Bob2(RESULT, rand, k, dim);
}

int** protocol::STopk_PAI(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, boundary* node, paillier_ciphertext_t** max_val, int NumData, int NumNode)
{
	printf("\n=== STopk_PAI start ===\n");
	int i=0, s=0, n=0, t=0, j=0, m=0;
	int rand = 5;
	int cnt = 0;	 // 질의 영역과 겹치는 노드 내에 존재하는 총 데이터의 수
	int NumNodeGroup = 0;
	Print = false;
	bool verify_flag = false;


	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}


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
	
	startTime = std::chrono::system_clock::now(); // startTime check	
	hint = q[dim];

	for( i = 0 ; i < dim ; i++ )
	{
		paillier_mul(pubkey, coeff[i], q[i], hint);
	}
	for( i = 0 ; i < dim ; i++ )
	{
		psi[i] = AS_CMP_MIN_BOOL(hint, coeff[i]);
	}
	for( i = 0 ; i < dim ; i++ )
	{	
		Q[i] = SM_p1(max_val[i], psi[i]);
	}
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["Query"] = time_variable.find("Query")->second + duration_sec.count();


	startTime = std::chrono::system_clock::now(); // startTime check	
	PARALLEL_AS_CMP_SRO_Topk(Q, Q, alpha, node, NumNode, &cnt, &NumNodeGroup, true);
	//Parallel_GSRO_Topk(Q, Q, alpha, node, NumNode, &cnt, &NumNodeGroup, true);
	endTime = std::chrono::system_clock::now(); // endTime check
	duration_sec = endTime - startTime;  // calculate duration 
	time_variable["nodeRetrievalSRO"] = time_variable.find("nodeRetrievalSRO")->second + duration_sec.count();

	while(1)
	{
		cnt = 0;
		if(!verify_flag)
		{	
			startTime = std::chrono::system_clock::now(); // startTime check	
			cand = PARALLEL_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			cout << "NumNodeGroup : " << NumNodeGroup << endl;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["nodeRetrieval"] = time_variable.find("nodeRetrieval")->second + duration_sec.count();
		}
		else
		{
			startTime = std::chrono::system_clock::now(); // startTime check	
			temp_cand = PARALLEL_sNodeRetrievalforTopk(data, node, alpha, NumData, NumNode, &cnt, &NumNodeGroup);
			totalNumOfRetrievedNodes += NumNodeGroup;
			cout << "Expand NumNodeGroup : " << NumNodeGroup << endl;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["exnodeRetrieval"] = time_variable.find("exnodeRetrieval")->second + duration_sec.count();

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

		startTime = std::chrono::system_clock::now(); // startTime check	
		ComputeScoreinMultithread(cand, q, SCORE, &cnt, true);
		endTime = std::chrono::system_clock::now(); // endTime check
		duration_sec = endTime - startTime;  // calculate duration 
		time_variable["CS"] = time_variable.find("CS")->second + duration_sec.count();

		for( s = 0 ; s < k ; s++ )
		{
			cout << s+1<<"dth MAXn_Topk start"<<endl;
			startTime = std::chrono::system_clock::now(); // startTime check	

			MAX = AS_CMP_MAXn(SCORE, cnt);
			//gmp_printf("Max : %Zd \n", paillier_dec(0, pubkey, prvkey, MAX));

			if(!verify_flag) {
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["SMAX"] = time_variable.find("SMAX")->second + duration_sec.count();
			}
			else {
				endTime = std::chrono::system_clock::now(); // endTime check
				duration_sec = endTime - startTime;  // calculate duration 
				time_variable["exSMAX"] = time_variable.find("exSMAX")->second + duration_sec.count();
			}
			startTime = std::chrono::system_clock::now(); // startTime check	
			MAXnMultithread2(cnt, SCORE_MINUS_MAX, SCORE, MAX, C_RAND);
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["subtract"] = time_variable.find("subtract")->second + duration_sec.count();

			V = Topk_sub(SCORE_MINUS_MAX, cnt);

			startTime = std::chrono::system_clock::now(); // startTime check	
			MAXnMultithread3(cnt, s, V, V2, SCORE, cand, MIN, RESULT);
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["Topk&update"] = time_variable.find("Topk&update")->second + duration_sec.count();
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
			if(cnt < k)	 // 노드의 fanout이 k보다 적은 경우를 처리함
				kth_score = computeScore(q, RESULT[cnt-1]);
			else	// k개가 찾아졌다면, k번째 데이터의 score를 계산
				kth_score = computeScore(q, RESULT[k-1]);

			startTime = std::chrono::system_clock::now(); // startTime check	
			int threadNum = thread_num;
			if(threadNum > NumNode)
			{
				threadNum = NumNode;
			}
			std::vector<int> *inputIdx = new std::vector<int>[threadNum];
			int count = 0;

			for(int i = 0; i < NumNode ; i++)
			{
				inputIdx[count++].push_back(i);
				count %= threadNum;
			}					
			std::thread *GSCMPforTopk_Thread = new std::thread[threadNum];
			for(int i = 0; i < threadNum; i++)
			{
				GSCMPforTopk_Thread[i] = std::thread(GSCMPforTopk_inThread, std::ref(inputIdx[i]), std::ref(node), std::ref(psi), std::ref(kth_score), std::ref(q), std::ref(alpha), std::ref(*this));
			}
			for(int i = 0; i < threadNum; i++)
			{
				GSCMPforTopk_Thread[i].join();
			}
			delete[] GSCMPforTopk_Thread;
			delete[] inputIdx;
			endTime = std::chrono::system_clock::now(); // endTime check
			duration_sec = endTime - startTime;  // calculate duration 
			time_variable["excheck"] = time_variable.find("excheck")->second + duration_sec.count();
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
		gmp_printf(" %Zd \n", paillier_dec(0, pubkey, prvkey, computeScore(q, RESULT[i])));
	
		for( j = 0 ; j < dim ; j++ )
		{
			paillier_mul(pubkey, RESULT[i][j], RESULT[i][j], C_RAND);
		}
	}
	cout << "End Line" <<endl;
	return SkNNm_Bob2(RESULT, rand, k, dim);
}



void protocol::TOPK_CS_SBD_inMultithread(paillier_ciphertext_t*** data, paillier_ciphertext_t** q, paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_distance)
{
	cout << "TOPK_CS_SBD_inMultithread" << endl;
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
	
	
	std::thread *CS_SBD_thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		CS_SBD_thread[i] = std::thread(CS_SBD_inThread, std::ref(inputIdx[i]), std::ref(q), std::ref(data), std::ref(cipher_distance), std::ref(cipher_SBD_distance), std::ref(*this));
	}

	for(int i = 0; i < threadNum; i++)
	{
		CS_SBD_thread[i].join();
	}
	
	delete[] CS_SBD_thread;
	delete[] inputIdx;
}

paillier_ciphertext_t** protocol::Smax_n_Multithread(paillier_ciphertext_t*** cipher, int number)
{
	cout << "Smax_n_Multithread" << endl;
	paillier_ciphertext_t*** copy_cipher = new paillier_ciphertext_t **[number];
	for(int i = 0; i < number ; i++)
	{
		copy_cipher[i] = cipher[i];
	}

	paillier_ciphertext_t **min;
	
	int numThread = thread_num;
	if(numThread > number)
	{
		numThread = number;
	}

	std::vector<int> *inputIdx = new std::vector<int>[numThread];
	paillier_ciphertext_t ***inputMax = new paillier_ciphertext_t **[numThread];

	int count = 0;
	for(int i = 0; i < number; i++)
	{

		inputIdx[count++].push_back(i);
		count %= numThread;
	}

	std::thread *smaxThread = new std::thread[numThread];
	for(int i = 0; i < numThread; i++)
	{
		smaxThread[i] = std::thread(Smax_n_InThread, i, std::ref(inputIdx[i]), std::ref(cipher), std::ref(inputMax), std::ref(*this));
	}

	for(int i = 0; i < numThread; i++)
	{
		smaxThread[i].join();
	}

	min = Smax_n(inputMax, numThread);

	delete[] inputIdx;
	delete[] smaxThread;
	delete[] copy_cipher;

	return min;
}

void protocol::recalculate_SCORE_inMultithread(paillier_ciphertext_t*** cipher_SBD_dist, paillier_ciphertext_t** cipher_distance, int NumData)
{
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
	std::thread *recalculate_dist_Thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		recalculate_dist_Thread[i] = std::thread(recalculate_dist_inThread, std::ref(inputIdx[i]), std::ref(cipher_SBD_dist), std::ref(cipher_distance), std::ref(*this));
	}
	for(int i = 0; i < threadNum; i++)
	{
		recalculate_dist_Thread[i].join();
	}
	delete[] recalculate_dist_Thread;
	delete[] inputIdx;
}

void protocol::DuringProcessedInTOPK_PB_inMultithread(paillier_ciphertext_t** cipher_distance, paillier_ciphertext_t** cipher_mid, paillier_ciphertext_t* cipher_min, paillier_ciphertext_t* cipher_rand, int NumData)
{
	cout << "DuringProcessedInTOPK_PB_inMultithread" << endl;
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
	std::thread *DuringProcessedInTOPK_PB_Thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		DuringProcessedInTOPK_PB_Thread[i] = std::thread(DuringProcessedInTOPK_PB_inThread, std::ref(inputIdx[i]), std::ref(cipher_distance), std::ref(cipher_mid), std::ref(cipher_min), std::ref(cipher_rand), std::ref(*this));
	}
	for(int i = 0; i < threadNum; i++)
	{
		DuringProcessedInTOPK_PB_Thread[i].join();
	}
	delete[] DuringProcessedInTOPK_PB_Thread;
	delete[] inputIdx;
}

void protocol::ExtractTOPKInTOPK_PB_inMultithread(paillier_ciphertext_t*** data, paillier_ciphertext_t** cipher_V, paillier_ciphertext_t*** cipher_result, int NumData, int serialTOPK)
{
	cout << "ExtractTOPKInTOPK_PB_inMultithread" << endl;
	int threadNum = thread_num;
	if(threadNum > NumData)
	{
		threadNum = NumData;
	}

	paillier_ciphertext_t*** cipher_thread_result = new paillier_ciphertext_t**[threadNum];
	for( int i = 0 ; i < threadNum ; i++ )
	{
		cipher_thread_result[i] = new paillier_ciphertext_t*[dim];
		for( int j = 0 ; j < dim ; j++ )
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
	std::thread *ExtractTOPKInTOPK_PB_Thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		ExtractTOPKInTOPK_PB_Thread[i] = std::thread(ExtractTOPKInTOPK_PB_inThread, std::ref(inputIdx[i]), std::ref(data), std::ref(cipher_V), std::ref(cipher_thread_result[i]), std::ref(*this));
	}
	for(int i = 0; i < threadNum; i++)
	{
		ExtractTOPKInTOPK_PB_Thread[i].join();
	}


	for( int i = 0 ; i < threadNum ; i++ )
	{
		for( int j = 0 ; j < dim ; j++ )
		{
			paillier_mul(pubkey, cipher_result[serialTOPK][j], cipher_result[serialTOPK][j], cipher_thread_result[i][j]);
		}
	}

	delete[] cipher_thread_result;
	delete[] ExtractTOPKInTOPK_PB_Thread;
	delete[] inputIdx;
}

void protocol::UPDATE_SBD_SCORE_InTOPK_PB_inMultithread(paillier_ciphertext_t*** cipher_SBD_distance, paillier_ciphertext_t** cipher_V, int NumData)
{
	cout << "UPDATE_SBD_SCORE_InTOPK_PB_inMultithread" << endl;
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
	std::thread *UPDATE_SBD_SCORE_InTOPK_PB_Thread = new std::thread[threadNum];
	for(int i = 0; i < threadNum; i++)
	{
		UPDATE_SBD_SCORE_InTOPK_PB_Thread[i] = std::thread(UPDATE_SBD_SCORE_InTOPK_PB_inThread, std::ref(inputIdx[i]), std::ref(cipher_SBD_distance), std::ref(cipher_V), std::ref(*this));
	}
	for(int i = 0; i < threadNum; i++)
	{
		UPDATE_SBD_SCORE_InTOPK_PB_Thread[i].join();
	}

	delete[] UPDATE_SBD_SCORE_InTOPK_PB_Thread;
	delete[] inputIdx;
}

paillier_ciphertext_t*** protocol::PARALLEL_sNodeRetrievalforTopk(paillier_ciphertext_t*** data, boundary* node, paillier_ciphertext_t** alpha, int NumData, int NumNode, int* cnt, int* NumNodeGroup)
{
	printf("\n===== Now PARALLEL_sNodeRetrievalforTopk starts =====\n");
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
		GSROThread[i] = std::thread(GSRO_MultithreadforTopk, std::ref(cipher_qLL), std::ref(cipher_qLL), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), true);
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

//PARALLEL_AS_CMP_SRO
void protocol::PARALLEL_AS_CMP_SRO_Topk(paillier_ciphertext_t** cipher_qLL, paillier_ciphertext_t** cipher_qRR, paillier_ciphertext_t** alpha, boundary* node, int NumNode, int* cnt, int* NumNodeGroup, bool type)
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
		GSROThread[i] = std::thread(GSRO_MultithreadforTopk, std::ref(cipher_qLL), std::ref(cipher_qLL), std::ref(NumNodeInput[i]), std::ref(node), std::ref(alpha), std::ref(*this), true);
	}
	for(int i = 0; i < NumThread; i++)
	{
		GSROThread[i].join();
	}
	delete[] NumNodeInput;
	delete[] GSROThread;
}
