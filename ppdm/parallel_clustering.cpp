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

void protocol::Parallel_Compute_Cluster(int NumData, paillier_ciphertext_t*** cipher, paillier_ciphertext_t*** former_Center, paillier_ciphertext_t*** NewSumCluster, paillier_ciphertext_t** NewSumCntCluster)
{
	cout << "Parallel_Compute_Cluster" << endl;
	int NumThread = thread_num;
	if(NumThread > NumData )
	{
		NumThread = NumData;
	}
	std::vector<int> *NumNodeInput = new std::vector<int>[NumThread];
	int count = 0;
	for(int i = 0; i < NumData ; i++)
	{
		NumNodeInput[count++].push_back(i);
		count %= NumThread;
	}
	std::thread *ComputeCluster = new std::thread[NumThread];
	for(int i = 0; i < NumThread; i++)
	{
		ComputeCluster[i] = std::thread(Comp_Cluster_inThread, std::ref(cipher), std::ref(former_Center), std::ref(NumNodeInput[i]), std::ref(NewSumCluster), std::ref(NewSumCntCluster), std::ref(*this), i);
	}
	for(int i = 0; i < NumThread; i++)
	{
		ComputeCluster[i].join();
	}
	delete[] NumNodeInput;
	delete[] ComputeCluster;
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

		////////////////////////////
		Parallel_Compute_Cluster(NumData, cipher, formerSumCluster, NewSumCluster, NewSumCntCluster);
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
	return 0;
}
paillier_ciphertext_t*** protocol::Clustering_PGI(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k)
{
	printf("Start Clustering_Grid_preprocessing\n");
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


