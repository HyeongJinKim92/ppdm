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

// KMEANS_PB CALCULATE THE NEW CLUSTER

void protocol::ComputeNEWCLUSTER_forKMEANS_PBinMultithread(paillier_ciphertext_t*** data, paillier_ciphertext_t*** formerSumCluster, paillier_ciphertext_t*** NewSumCluster, paillier_ciphertext_t** NewSumCntCluster)
{
	cout << "ComputeNEWCLUSTER_forKMEANS_PB" << endl;
	int NumThread = thread_num;
	if(NumThread > NumData )
	{
		NumThread = NumData;
	}


	paillier_ciphertext_t**** tmp_NewSumCluster = new paillier_ciphertext_t***[NumThread];
	paillier_ciphertext_t*** tmp_NewSumCntCluster = new paillier_ciphertext_t**[NumThread];
	for ( int i = 0 ; i < NumThread ; i++ )
	{
		tmp_NewSumCluster[i] = new paillier_ciphertext_t**[k];
		tmp_NewSumCntCluster[i] = new paillier_ciphertext_t*[k];
		for ( int j = 0 ; j < k ; j++ )
		{
			tmp_NewSumCluster[i][j] = new paillier_ciphertext_t*[dim];
			tmp_NewSumCntCluster[i][j] = paillier_create_enc(0);
			for ( int l = 0 ; l < dim ; l++ )
			{
				tmp_NewSumCluster[i][j][l] = paillier_create_enc(0);
			}
		}
	}


	std::vector<int> *Input = new std::vector<int>[NumThread];
	int count = 0;
	for(int i = 0; i < NumData ; i++)
	{
		Input[count++].push_back(i);
		count %= NumThread;
	}
	std::thread *ComputeNEWCLUSTER_forKMEANS_PBThread = new std::thread[NumThread];
	for(int i = 0; i < NumThread; i++)
	{
		ComputeNEWCLUSTER_forKMEANS_PBThread[i] = std::thread(ComputeNEWCLUSTER_forKMEANS_PBinThread, std::ref(Input[i]), std::ref(data), std::ref(formerSumCluster),  std::ref(tmp_NewSumCluster[i]), std::ref(tmp_NewSumCntCluster[i]), std::ref(*this));
	}
	for(int i = 0; i < NumThread; i++)
	{
		ComputeNEWCLUSTER_forKMEANS_PBThread[i].join();
	}

	for ( int i = 0 ; i < NumThread ; i++ )
	{
		for ( int j = 0 ; j < k ; j++ )
		{
			paillier_mul(pub, NewSumCntCluster[j], NewSumCntCluster[j], tmp_NewSumCntCluster[i][j]);
			for ( int l = 0 ; l < dim ; l++ )
			{
				paillier_mul(pub, NewSumCluster[j][l], NewSumCluster[j][l], tmp_NewSumCluster[i][j][l]);
			}
		}
	}

	delete[] tmp_NewSumCluster;
	delete[] tmp_NewSumCntCluster;
	delete[] Input;
	delete[] ComputeNEWCLUSTER_forKMEANS_PBThread;
}


void protocol::ComputeNEWCLUSTER_forKMEANS_PGIinMultithread(paillier_ciphertext_t*** data, paillier_ciphertext_t*** formerSumCluster, paillier_ciphertext_t*** NewSumCluster, paillier_ciphertext_t** NewSumCntCluster)
{
	cout << "ComputeNEWCLUSTER_forKMEANS_PGIinMultithread" << endl;
	int NumThread = thread_num;
	if(NumThread > NumData )
	{
		NumThread = NumData;
	}


	paillier_ciphertext_t**** tmp_NewSumCluster = new paillier_ciphertext_t***[NumThread];
	paillier_ciphertext_t*** tmp_NewSumCntCluster = new paillier_ciphertext_t**[NumThread];
	for ( int i = 0 ; i < NumThread ; i++ )
	{
		tmp_NewSumCluster[i] = new paillier_ciphertext_t**[k];
		tmp_NewSumCntCluster[i] = new paillier_ciphertext_t*[k];
		for ( int j = 0 ; j < k ; j++ )
		{
			tmp_NewSumCluster[i][j] = new paillier_ciphertext_t*[dim];
			tmp_NewSumCntCluster[i][j] = paillier_create_enc(0);
			for ( int l = 0 ; l < dim ; l++ )
			{
				tmp_NewSumCluster[i][j][l] = paillier_create_enc(0);
			}
		}
	}


	std::vector<int> *Input = new std::vector<int>[NumThread];
	int count = 0;
	for(int i = 0; i < NumData ; i++)
	{
		Input[count++].push_back(i);
		count %= NumThread;
	}
	std::thread *ComputeNEWCLUSTER_forKMEANS_PBThread = new std::thread[NumThread];
	for(int i = 0; i < NumThread; i++)
	{
		ComputeNEWCLUSTER_forKMEANS_PBThread[i] = std::thread(ComputeNEWCLUSTER_forKMEANS_PBinThread, std::ref(Input[i]), std::ref(data), std::ref(formerSumCluster),  std::ref(tmp_NewSumCluster[i]), std::ref(tmp_NewSumCntCluster[i]), std::ref(*this));
	}
	for(int i = 0; i < NumThread; i++)
	{
		ComputeNEWCLUSTER_forKMEANS_PBThread[i].join();
	}

	for ( int i = 0 ; i < NumThread ; i++ )
	{
		for ( int j = 0 ; j < k ; j++ )
		{
			paillier_mul(pub, NewSumCntCluster[j], NewSumCntCluster[j], tmp_NewSumCntCluster[i][j]);
			for ( int l = 0 ; l < dim ; l++ )
			{
				paillier_mul(pub, NewSumCluster[j][l], NewSumCluster[j][l], tmp_NewSumCluster[i][j][l]);
			}
		}
	}

	delete[] tmp_NewSumCluster;
	delete[] tmp_NewSumCntCluster;
	delete[] Input;
	delete[] ComputeNEWCLUSTER_forKMEANS_PBThread;
}



paillier_ciphertext_t*** protocol::Parallel_Clustering_m_Grid(paillier_ciphertext_t*** cipher, int NumData, int k, int b, paillier_ciphertext_t*** former_Center)
{
	printf("Parallel_Clustering_m_Grid start\n");
	
	paillier_ciphertext_t*** origin_data				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	for(int i = 0 ; i < NumData ; i++ ){
		origin_data[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j++){
			origin_data[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(origin_data[i][j]->c);
			origin_data[i][j] = cipher[i][j];
		}
	}

	paillier_ciphertext_t*** Distance_center_data		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	paillier_ciphertext_t** minDist						= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumData);

	paillier_ciphertext_t*** formerSumCluster			= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t** formerSumCntCluster			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	paillier_ciphertext_t*** NewSumCluster				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t** NewSumCntCluster			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	

	paillier_ciphertext_t** data_arr					= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	paillier_ciphertext_t** idx_arr;
	paillier_ciphertext_t* beta							= paillier_create_enc(b);
	paillier_ciphertext_t* f							= paillier_create_enc(1);
	paillier_ciphertext_t** f_arr						= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	paillier_plaintext_t* gamma;


	for(int i = 0 ; i < k ; i++ ){
		f_arr[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		NewSumCntCluster[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		formerSumCntCluster[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		data_arr[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		
		mpz_init(f_arr[i]->c);
		mpz_init(data_arr[i]->c);
		mpz_init(NewSumCntCluster[i]->c);
		mpz_init(formerSumCntCluster[i]->c);
		
		NewSumCntCluster[i] = paillier_create_enc_zero();
		formerSumCntCluster[i] = paillier_create_enc(1);
		NewSumCluster[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		formerSumCluster[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j++){
			NewSumCluster[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			formerSumCluster[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(NewSumCluster[i][j]->c);
			mpz_init(formerSumCluster[i][j]->c);
		}
	}

	for(int i = 0 ; i < NumData ; i++){
		minDist[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(minDist[i]->c);
		Distance_center_data[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
		for( int j = 0 ; j < k ; j++){
			Distance_center_data[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		}
	}

	srand(time(NULL));
	
	for(int i = 0 ; i < k ; i++ ){
		for(int j = 0 ; j < dim ; j++ ){
			formerSumCluster[i][j] = former_Center[i][j];
		}
	}
	
	for( int i = 0 ; i < k ; i++ ){
		for( int j = 0 ; j < dim ; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, formerSumCluster[i][j]));
		}
		printf("\n");
	}
	
	
	for(int i = 0 ; i < 1000 ; i++)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ROUND : %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", i);
		int t;
		for(int i = 0 ; i < k ; i++ ){
			NewSumCntCluster[i] = paillier_create_enc(0);
			for(int j = 0 ; j < dim ; j++){
				NewSumCluster[i][j] = paillier_create_enc(0);
			}
		}
		
		ComputeNEWCLUSTER_forKMEANS_PGIinMultithread(cipher, formerSumCluster, NewSumCluster, NewSumCntCluster);
		//Parallel_Compute_Cluster(NumData, cipher, formerSumCluster, NewSumCluster, NewSumCntCluster);
		////////////////////////////
/*
		for(int i = 0 ; i < NumData ; i++){

			for( int j = 0 ; j < k ; j++){
				Distance_center_data[i][j] = SSEDm(cipher[i], formerSumCluster[j], dim);
				//Distance_center_data[i][j] = DP_SSED(cipher[i], formerSumCluster[j], dim);
				//Distance_center_data[i][j] = OP_SSED(cipher[i], formerSumCluster[j], formerSumCntCluster, i, k);
				//gmp_printf("dist : %Zd\n", paillier_dec(0, pubkey, prvkey, Distance_center_data[i][j]));
			}
			
			idx_arr = Smin_bool(Distance_center_data[i], k);
			
			
			for( int j = 0 ; j < k ; j++ ){
				paillier_mul(pubkey, NewSumCntCluster[j], idx_arr[j], NewSumCntCluster[j]);
				for( int e = 0 ; e < dim ; e++ ){
					data_arr[e] = SM_p1(idx_arr[j], origin_data[i][e]);
					paillier_mul(pubkey, NewSumCluster[j][e], NewSumCluster[j][e], data_arr[e]);
				}
			}
		}
*/
		////////////////////////////
		printf("save point 1\n");
		
		printf("new center cnt\n");
		for(int i = 0 ; i < k ; i++ ){
			gmp_printf("%Zd\t", paillier_dec(0, pub, prv, NewSumCntCluster[i]));
		}
		printf("\n");


		unsigned int convert_int;
		unsigned int convert_cnt_int;
	
		paillier_plaintext_t* convert_center;
		paillier_plaintext_t* convert_cnt;

		for(int i = 0 ; i < k ; i++ ){
			convert_cnt = paillier_dec(0, pubkey, prvkey, NewSumCntCluster[i]);
			convert_cnt_int = mpz_get_ui(convert_cnt->m);
			for(int j = 0 ; j < dim ; j++ ){
				convert_center = paillier_dec(0, pubkey, prvkey, NewSumCluster[i][j]);
				convert_int = mpz_get_ui(convert_center->m);
				convert_int = convert_int/convert_cnt_int;
				NewSumCluster[i][j] = paillier_create_enc(convert_int);
			}
			NewSumCntCluster[i] = paillier_create_enc(1);
		}
		gamma = SETC(formerSumCluster, formerSumCntCluster, NewSumCluster, NewSumCntCluster);
		
		if( mpz_get_ui(gamma->m) == 1 ){
			break;
		}else{
			for(int i = 0 ; i < k ; i ++){
				for( int j = 0 ; j < dim ; j++ ){
					formerSumCluster[i][j] = NewSumCluster[i][j];
				}
				formerSumCntCluster[i] = NewSumCntCluster[i];
			}
		}
	}
/*
	for(int i = 0 ; i < k ; i++ ){
		printf("Sum : ");
		for(int j = 0 ; j < dim ; j++){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, NewSumCluster[i][j]));
		}
		printf("\t\t");
	}
	printf("\n");
*/
	return NewSumCluster;
}



paillier_ciphertext_t*** protocol::Clustering_PB(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k)
{
	printf("!!!!!!!!!!!!!!!!!!!!!!!!Clustering_PB start!!!!!!!!!!!!!!!!!!!!!!!!\n");
	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}

	
	paillier_ciphertext_t*** origin_data				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	for(int i = 0 ; i < NumData ; i++ ){
		origin_data[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j++){
			origin_data[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(origin_data[i][j]->c);
			origin_data[i][j] = ciper[i][j];
		}
	}

	paillier_ciphertext_t*** formerSumCluster			= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t** formerSumCntCluster			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	paillier_ciphertext_t*** NewSumCluster				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t** NewSumCntCluster			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	

	paillier_ciphertext_t** data_arr					= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	paillier_ciphertext_t** idx_arr;
	paillier_ciphertext_t* beta							= paillier_create_enc(k);
	paillier_ciphertext_t* f							= paillier_create_enc(1);
	paillier_ciphertext_t** f_arr						= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	paillier_plaintext_t* gamma;


	for(int i = 0 ; i < k ; i++ ){
		f_arr[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		NewSumCntCluster[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		formerSumCntCluster[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		data_arr[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		
		mpz_init(f_arr[i]->c);
		mpz_init(data_arr[i]->c);
		mpz_init(NewSumCntCluster[i]->c);
		mpz_init(formerSumCntCluster[i]->c);
		
		NewSumCntCluster[i] = paillier_create_enc_zero();
		formerSumCntCluster[i] = paillier_create_enc(1);
		NewSumCluster[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		formerSumCluster[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j++){
			NewSumCluster[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			formerSumCluster[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(NewSumCluster[i][j]->c);
			mpz_init(formerSumCluster[i][j]->c);
		}
	}


	srand(time(NULL));

	for(int i = 0 ; i < k ; i++ ){
		for(int j = 0 ; j < dim ; j++ ){
			formerSumCluster[i][j] = origin_data[rand()%NumData][j];
		}
	}
	
	for(int i = 0 ; i < 100 ; i++)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ROUND : %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", i);
		int t;
		for(int i = 0 ; i < k ; i++ ){
			NewSumCntCluster[i] = paillier_create_enc(0);
			for(int j = 0 ; j < dim ; j++){
				NewSumCluster[i][j] = paillier_create_enc(0);
			}
		}
		/////////////////////////////////

		ComputeNEWCLUSTER_forKMEANS_PBinMultithread(origin_data, formerSumCluster, NewSumCluster, NewSumCntCluster);

/*
		for(int i = 0 ; i < NumData ; i++){
			for( int j = 0 ; j < k ; j++){
				Distance_center_data[i][j] = SSEDm(origin_data[i], formerSumCluster[j], dim);
			}
			
			idx_arr = Smin_bool(Distance_center_data[i], k);
			
			for( int j = 0 ; j < k ; j++ ){
				paillier_mul(pubkey, NewSumCntCluster[j], idx_arr[j], NewSumCntCluster[j]);
				for( int e = 0 ; e < dim ; e++ ){
					data_arr[e] = SM_p1(idx_arr[j], origin_data[i][e]);
					paillier_mul(pubkey, NewSumCluster[j][e], NewSumCluster[j][e], data_arr[e]);
				}
			}
		}
*/
		/////////////////////////////////

		printf("new center cnt\n");
		for(int i = 0 ; i < k ; i++ ){
			gmp_printf("%Zd\t", paillier_dec(0, pub, prv, NewSumCntCluster[i]));
		}
		printf("\n");


		unsigned int convert_int;
		unsigned int convert_cnt_int;
	
		paillier_plaintext_t* convert_center;
		paillier_plaintext_t* convert_cnt;

		for(int i = 0 ; i < k ; i++ ){
			convert_cnt = paillier_dec(0, pubkey, prvkey, NewSumCntCluster[i]);
			convert_cnt_int = mpz_get_ui(convert_cnt->m);
			for(int j = 0 ; j < dim ; j++ ){
				convert_center = paillier_dec(0, pubkey, prvkey, NewSumCluster[i][j]);
				convert_int = mpz_get_ui(convert_center->m);
				convert_int = convert_int/convert_cnt_int;
				NewSumCluster[i][j] = paillier_create_enc(convert_int);
			}
			NewSumCntCluster[i] = paillier_create_enc(1);
		}

		gamma = SETC(formerSumCluster, formerSumCntCluster, NewSumCluster, NewSumCntCluster);
		
		if( mpz_get_ui(gamma->m) == 1 ){
			break;
		}else{
			for(int i = 0 ; i < k ; i ++){
				for( int j = 0 ; j < dim ; j++ ){
					formerSumCluster[i][j] = NewSumCluster[i][j];
				}
				formerSumCntCluster[i] = NewSumCntCluster[i];
			}
		}
		for( int i = 0 ; i < k ; i++ ){
			for( int j = 0 ; j < dim ; j++ ){
				gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, formerSumCluster[i][j]));
			}
			printf("\n");
		}
	}
	return NewSumCluster;
}






paillier_ciphertext_t*** protocol::Clustering_PGI(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k)
{
	printf("Start Clustering_Grid_preprocessing\n");

	if( thread_num < 1 ) 
	{
		cout << "THREAD NUM IS ERROR" << endl;
		exit(1);
	}

	int i = 0, j = 0, m = 0;
											
	paillier_ciphertext_t*** formerCenter = PreClustering(ciper, node, NumNode, NumData, k);
	

	for( i = 0 ; i < k ; i++ ){
		for( j = 0 ; j < dim ; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, formerCenter[i][j]));
		}
		printf("\n");
	}

	paillier_ciphertext_t*** NewCluster = Parallel_Clustering_m_Grid(ciper, NumData, k, 100, formerCenter);
	

	printf("End Clustering_Grid\n");
	return NewCluster;

}
paillier_ciphertext_t*** protocol::Clustering_PAI(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k)
{
	return 0;	
}


