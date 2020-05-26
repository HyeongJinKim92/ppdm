#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <iostream>
#include <time.h>
#include "protocol.h"



paillier_ciphertext_t*** protocol::Clustering_m(paillier_ciphertext_t*** ciper, int NumData, int k, int b){
	printf("!!!!!!!!!!!!!!!!!!!!!!!!Clustering_m start!!!!!!!!!!!!!!!!!!!!!!!!\n");

	
	paillier_ciphertext_t*** origin_data				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	for(int i = 0 ; i < NumData ; i++ ){
		origin_data[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j++){
			origin_data[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(origin_data[i][j]->c);
			origin_data[i][j] = ciper[i][j];
		}
	}
	/*
	for(int i = 0 ; i < NumData ; i++ ){
		for(int j = 0 ; j < dim ; j++){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, ciper[i][j]));
		}
		printf("\n");
	}
	*/
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
		/*
		printf("former CNT\n");
		for(int i = 0 ; i < k; i++){
			gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, formerSumCntCluster[i]));
		}
		printf("\n");

		printf("former Center\n");
		for(int i = 0 ; i < k; i++){
			for(int j = 0 ; j < dim ; j++){
				gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, formerSumCluster[i][j]));
			}
			printf("\n");
		}
		printf("\n");
		*/
		
		/*
		printf("New CNT : ");
		for(int i = 0 ; i < k; i++){
			gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, NewSumCntCluster[i]));
		}
		printf("\n");

		printf("New Center\n");
		for(int i = 0 ; i < k; i++){
			for(int j = 0 ; j < dim ; j++){
				gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, NewSumCluster[i][j]));
			}
			printf("\n");
		}
		printf("\n");
		*/
		
		for(int i = 0 ; i < NumData ; i++){
						
			for( int j = 0 ; j < k ; j++){
				Distance_center_data[i][j] = SSEDm(origin_data[i], formerSumCluster[j], dim);
				//Distance_center_data[i][j] = OP_SSED(ciper[i], formerSumCluster[j], formerSumCntCluster, i, k);
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

paillier_ciphertext_t* protocol::OP_SSED(paillier_ciphertext_t** data, paillier_ciphertext_t** Cluster_Center, paillier_ciphertext_t** Center_cnt, int tmp_k, int k){
	//printf("Start op_ssed");
	paillier_ciphertext_t* b = paillier_create_enc(1);
	paillier_ciphertext_t* all_b = paillier_create_enc(1);

	paillier_ciphertext_t** a  = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim));
	paillier_ciphertext_t** a_ = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*(dim));
	
	paillier_ciphertext_t* result_dist = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));

	for(int i = 0 ; i < dim ; i++){
		a[i]  = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		a_[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	}
	
	for(int i = 0 ; i < k ; i++ ){
		all_b = SM_p1(Center_cnt[i], all_b);
		
		if( i != tmp_k ){
			b = SM_p1(b, Center_cnt[i]);
		}
	}
	
	for(int i = 0 ; i < dim ; i++ ){
		a[i] = SM_p1(data[i], all_b);
		a_[i] = SM_p1(Cluster_Center[i], b);
	}

	result_dist = SSEDm(a, a_, dim);
	
	for(int i = 0 ; i < dim ; i++){
		free(a[i]);
		free(a_[i]);
	}
	free(a);
	free(a_);
	free(b);
	free(all_b);
	
	//printf("End op_ssed\n");

	return result_dist;
}

paillier_plaintext_t* protocol::SETC(paillier_ciphertext_t*** former, paillier_ciphertext_t** formerCnt, paillier_ciphertext_t*** New, paillier_ciphertext_t** NewCnt){
	paillier_ciphertext_t** tau = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	paillier_ciphertext_t** V = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	paillier_ciphertext_t** Z = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	paillier_ciphertext_t* Y = paillier_create_enc(1);
	
	paillier_ciphertext_t*** G = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t*** G_ = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);

	paillier_ciphertext_t** H = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	paillier_ciphertext_t** H_ = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	paillier_ciphertext_t* L = paillier_create_enc(1);
	paillier_ciphertext_t* R = paillier_create_enc(1);

	paillier_plaintext_t* beta = paillier_plaintext_from_ui(10);

	for(int i = 0 ; i < k ; i++ ){
		G[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		G_[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	}

	for(int i = 0 ; i < k ; i++ ){
		tau[i] = SM_p1(formerCnt[i], NewCnt[i]);
	}
	
	for(int i = 0 ; i < k ; i++ ){
		V[i] = paillier_create_enc(1);
		for(int j = 0 ; j < k ; j++ ){
			if(i != j){
				V[i] = SM_p1(V[i], tau[j]);
			}
		}
		Z[i] = SM_p1(V[i], V[i]);
	}
	
	Y = SM_p1(V[0], tau[0]);
	Y = SM_p1(Y, Y);

	for(int i = 0 ; i < k ; i++ ){
		for(int j = 0 ; j < dim ; j++ ){
			G[i][j] = SM_p1(former[i][j], NewCnt[i]);
			G_[i][j] = SM_p1(New[i][j], formerCnt[i]);
		}
	}

	for(int i = 0 ; i < k ; i++ ){
		H[i] = SSEDm(G[i], G_[i], dim);
		H_[i] = SM_p1(H[i], Z[i]);
	}

	for(int i = 0 ; i < k ; i++ ){
		L = SM_p1(H_[i], L);
	}
	paillier_exp(pubkey, R, Y, beta);
	
	gmp_printf("L : %Zd\n", paillier_dec(0, pubkey, prvkey, L));
	gmp_printf("R : %Zd\n", paillier_dec(0, pubkey, prvkey, R));
	gmp_printf("SC : %Zd\n", paillier_dec(0, pubkey, prvkey, SC(L,R)));
	return paillier_dec(0, pubkey, prvkey, SC(L,R));
}

paillier_ciphertext_t*** protocol::Clustering_Grid(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k){
	printf("!!!!!!!!!!!!!!!!!!!!!!Start Clustering_Grid!!!!!!!!!!!!!!!!!!!!!!\n");
	
	paillier_ciphertext_t*** origin_data				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	for(int i = 0 ; i < NumData ; i++ ){
		origin_data[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j++){
			origin_data[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(origin_data[i][j]->c);
			origin_data[i][j] = ciper[i][j];
		}
	}

	
	paillier_ciphertext_t*** formercenter				= PreClustering(ciper, node, NumNode, NumData, k);
	

	paillier_ciphertext_t** formercenter_cnt			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);

	paillier_ciphertext_t*** newcenter					= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*k);
	paillier_ciphertext_t** newcenter_cnt				= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	
	paillier_ciphertext_t** data_arr			= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
	paillier_ciphertext_t** idx_arr;
	


	paillier_ciphertext_t*		temp_sum		= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	paillier_ciphertext_t***	GSRO_alpha		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumNode);
	paillier_ciphertext_t**		Sum_GSRO_alpha	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*NumNode);
	
	paillier_ciphertext_t**		tmp_sum_center	= (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	
	paillier_ciphertext_t*		filter_signal	= (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	
	
	for(int i = 0 ; i < k ; i++ ){
		formercenter_cnt[i] = paillier_create_enc(1);
		data_arr[i] = paillier_create_enc(0);
	}

	paillier_ciphertext_t* tmp_max_dist_inSRO;
	paillier_ciphertext_t* temp_coord1 ;
	paillier_ciphertext_t* temp_coord2 ;
	paillier_ciphertext_t* temp_coord3 ;
	paillier_ciphertext_t** psi = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*3);
	
	paillier_ciphertext_t** GSCMP_alpha = (paillier_ciphertext_t**) malloc(sizeof(paillier_ciphertext_t*)*k);
	paillier_ciphertext_t** shortestPoint = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	for(int i=0; i<dim; i++){
		shortestPoint[i] = (paillier_ciphertext_t*) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(shortestPoint[i]->c);
	}

	paillier_ciphertext_t* temp_nodedist = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
	mpz_init(temp_nodedist->c);
	
	for( int i = 0 ; i < k ; i++ ){
		for( int j = 0 ; j < dim ; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, formercenter[i][j]));
		}
		printf("\n");
	}
	paillier_plaintext_t* gamma;
	int cnt = 0;
	while(1){
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ROUND : %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", cnt++);
		for(int i = 0 ; i < k ; i++){
			newcenter[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
			newcenter_cnt[i] = paillier_create_enc(0);
			for(int j = 0 ; j < dim ; j++ ){
				newcenter[i][j] = paillier_create_enc(0);
			}
		}

		for(int i = 0 ; i < NumNode ; i++ ){
			GSRO_alpha[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
			Sum_GSRO_alpha[i] = paillier_create_enc(0);
		}


		for(int i = 0 ; i < k ; i++){
			for(int j = 0 ; j < NumNode ; j++){
				GSRO_alpha[j][i] = DP_GSRO(formercenter[i], formercenter[i], node[j].LL, node[j].RR);
				paillier_mul(pubkey, Sum_GSRO_alpha[j], Sum_GSRO_alpha[j], GSRO_alpha[j][i]);	
			}
		}
		for(int i = 0 ; i < NumNode ; i++ ){
			//printf("NumNode : %d\n", i);
			if( strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, Sum_GSRO_alpha[i])), paillier_plaintext_to_str(plain_one)) == 0){
				//printf("0000 Box in Point\n");
				for(int j = 0 ; j < k ; j++){
					for(int x = 0 ; x < dim ; x++){
						if( j == 0 ){
							tmp_sum_center[x] = SM_p1(GSRO_alpha[i][j], formercenter[j][x]);
						}else{
							paillier_mul(pubkey, tmp_sum_center[x], tmp_sum_center[x], SM_p1(GSRO_alpha[i][j], formercenter[j][x]));
						}
					}
				}
				tmp_max_dist_inSRO = boundary_dist(node[i], tmp_sum_center, NumNode);
				
				for(int m=0; m<k; m++){
					for(int j=0; j<dim; j++) {	
						psi[0] = SC(formercenter[m][j], node[i].LL[j]);		// 프사이 계산 (1이면 q의 좌표가 노드의 LL bound보다 크고, 0이면 q의 좌표가 작음)
						psi[1] = SC(formercenter[m][j], node[i].RR[j]);	// 프사이 계산 (1이면 q의 좌표가 노드의 UR bound보다 작고, 0이면 q의 좌표가 큼)
						psi[2] = SBXOR(psi[0], psi[1]);	// psi[2]가 1이면, 해당 차원에서의 질의 좌표가 노드의 LL bound 와 UR bound 사이에 속함을 의미
						temp_coord1 = SM_p1(psi[2], formercenter[m][j]);		// psi[2]가 1이면, 질의에서의 수선의 발이 최단거리 이므로, 질의의 좌표를 살려야 함
						temp_coord2 = SM_p1(psi[0], node[i].LL[j]);		// psi[0]이 0이면, LL bound의 좌표를 살림
						temp_coord3 = SM_p1(SBN(psi[0]), node[i].RR[j]);		// psi[0]이 1이면, RR bound의 좌표를 살림
						paillier_mul(pubkey, temp_coord3, temp_coord2, temp_coord3);
						temp_coord3 = SM_p1(SBN(psi[2]), temp_coord3);
						paillier_mul(pubkey, shortestPoint[j], temp_coord1, temp_coord3);
					}
					temp_nodedist = SSEDm(formercenter[m], shortestPoint, dim);
					GSCMP_alpha[m] = SC(tmp_max_dist_inSRO, temp_nodedist);
					paillier_mul(pubkey, GSCMP_alpha[m], GSCMP_alpha[m], GSRO_alpha[i][m]);
					if(m==0){
						filter_signal = GSCMP_alpha[m];
					}else{
						filter_signal = SM_p1(GSCMP_alpha[m],filter_signal);
					}
					//gmp_printf("nodedist :%Zd %Zd\t GSCMP_alpha : %Zd\n", paillier_dec(0, pub, prv, tmp_max_dist_inSRO), paillier_dec(0, pub, prv, temp_nodedist), paillier_dec(0, pub, prv, GSCMP_alpha[m]));
				}
				//gmp_printf("FILTERING SIGNAL : %Zd\n", paillier_dec(0, pub, prv, filter_signal));
				if( strcmp( paillier_plaintext_to_str( paillier_dec(0, pubkey, prvkey, filter_signal)), paillier_plaintext_to_str(plain_one)) == 0){
					printf("0000 1111 start FILTERING\n");
					idx_arr = GSRO_alpha[i];
					/*
					printf("idx_arr : ");
					for( int a = 0 ; a < k ; a++ ){
						gmp_printf("%Zd\t", paillier_dec(0, pub, prv, idx_arr[a]));
					}
					printf("\n");
					*/
					for( int l = 0 ; l < node[i].NumData ; l++){
						for( int j = 0 ; j < k ; j++ ){
							paillier_mul(pubkey, newcenter_cnt[j], idx_arr[j], newcenter_cnt[j]);
							//printf("data_arr : %d \t", j);
							for( int e = 0 ; e < dim ; e++ ){
								data_arr[e] = SM_p1(idx_arr[j], origin_data[node[i].indata_id[l]][e]);
								//gmp_printf("%Zd\t", paillier_dec(0, pub, prv, data_arr[e]));
								paillier_mul(pubkey, newcenter[j][e], newcenter[j][e], data_arr[e]);	
							}
							//printf("\n");
						}
					}
				}else{
					//printf("0000 2222 INCLUDE No FILTERING\n");
					paillier_ciphertext_t*** Distance_center_data		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*node[i].NumData);

					for( int j = 0 ; j < node[i].NumData ; j++){
						Distance_center_data[j] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
						for( int m = 0 ; m < k ; m++ ){
							Distance_center_data[j][m] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
							mpz_init(Distance_center_data[j][m]->c);
							Distance_center_data[j][m] = SSEDm(origin_data[node[i].indata_id[j]], formercenter[m], dim);
						}
						idx_arr = Smin_bool(Distance_center_data[j], k);
						
						for( int a = 0 ; a < k ; a++ ){
							paillier_mul(pubkey, newcenter_cnt[a], idx_arr[a], newcenter_cnt[a]);
							//printf("data_arr : %d \t", j);
							for( int e = 0 ; e < dim ; e++ ){
								data_arr[e] = SM_p1(idx_arr[a], origin_data[node[i].indata_id[j]][e]);
								//gmp_printf("%Zd\t", paillier_dec(0, pub, prv, data_arr[e]));
								paillier_mul(pubkey, newcenter[a][e], newcenter[a][e], data_arr[e]);
							}
							//printf("\n");
						}

					}
				}
			}else{
				//printf("3333 POINT is not included\n");
				paillier_ciphertext_t*** Distance_center_data		= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*node[i].NumData);
				for( int j = 0 ; j < node[i].NumData ; j++){
					Distance_center_data[j] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*k);
					for( int m = 0 ; m < k ; m++ ){
						Distance_center_data[j][m] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
						mpz_init(Distance_center_data[j][m]->c);
						Distance_center_data[j][m] = SSEDm( origin_data[node[i].indata_id[j]], formercenter[m], dim);
					}
					idx_arr = Smin_bool(Distance_center_data[j], k);
					
					for( int a = 0 ; a < k ; a++ ){
						paillier_mul(pubkey, newcenter_cnt[a], idx_arr[a], newcenter_cnt[a]);
						//printf("data_arr : %d \t", j);
						for( int e = 0 ; e < dim ; e++ ){
							data_arr[e] = SM_p1(idx_arr[a], origin_data[node[i].indata_id[j]][e]);
							//gmp_printf("%Zd\t", paillier_dec(0, pub, prv, data_arr[e]));
							paillier_mul(pubkey, newcenter[a][e], newcenter[a][e], data_arr[e]);
						}
						//printf("\n");
					}
				}
			}
		}
		/*
		printf("save point 1\n");
		
		printf("new center cnt\n");
		for(int i = 0 ; i < k ; i++ ){
			gmp_printf("%Zd\t", paillier_dec(0, pub, prv, newcenter_cnt[i]));
		}
		printf("\n");
		*/
		unsigned int convert_int;
		unsigned int convert_cnt_int;
		paillier_plaintext_t* convert_center;
		paillier_plaintext_t* convert_cnt;

		for(int i = 0 ; i < k ; i++ ){
			convert_cnt = paillier_dec(0, pubkey, prvkey, newcenter_cnt[i]);
			convert_cnt_int = mpz_get_ui(convert_cnt->m);
			//gmp_printf("convert_cnt : %Zd %d\n", convert_cnt->m, convert_cnt_int);
			for(int j = 0 ; j < dim ; j++ ){
				convert_center = paillier_dec(0, pubkey, prvkey, newcenter[i][j]);
				//gmp_printf("convert_center : %Zd\n", convert_center->m);
				convert_int = mpz_get_ui(convert_center->m);
				//gmp_printf("convert_center : %Zd %d\n", convert_center->m, convert_int);
				if(convert_cnt == 0){
					convert_int = 0;
				}
				convert_int = convert_int/convert_cnt_int;
				newcenter[i][j] = paillier_create_enc((int)convert_int);
			}
			newcenter_cnt[i] = paillier_create_enc(1);
		}
		
		gamma = SETC(formercenter, formercenter_cnt, newcenter, newcenter_cnt);
		
		if( mpz_get_ui(gamma->m) == 1 ){
			break;
		}else{
			for(int i = 0 ; i < k ; i ++){
				for( int j = 0 ; j < dim ; j++ ){
					formercenter[i][j] = newcenter[i][j];
				}
				formercenter_cnt[i] = newcenter_cnt[i];
			}
		}
	}

	printf("End Clustering_Grid\n");
	return newcenter;
}

paillier_ciphertext_t*** protocol::PreClustering(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k){
	printf("Start PreClustering\n");
	/*
	for(int i = 0 ; i < NumNode ; i++){
		printf("node Cnt : %d\n", node[i].NumData);
	}
	*/
	int Sample_Cnt = 0;

	paillier_ciphertext_t*** Sample = extract_sample(ciper, node, NumNode, NumData, &Sample_Cnt);
	/*
	printf("Sample_Cnt : %d\n", Sample_Cnt);
	printf("Sample Data\n");
	for(int i = 0 ; i < Sample_Cnt ; i++ ){
		for(int j = 0 ; j < dim ; j++){
			gmp_printf("%Zd ", paillier_dec(0, pub, prv, Sample[i][j]));
		}
		printf("\n");
	}
	*/
	paillier_ciphertext_t*** Pre_Cluster = Clustering_m(Sample, Sample_Cnt, k, 100);
/*
	for(int i = 0 ; i < k ; i++ ){
		printf("in PreClustering Sum : ");
		for(int j = 0 ; j < dim ; j++){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, Pre_Cluster[i][j]));
		}
		printf("\t\t");
	}
	printf("\n");
*/
	
	return Pre_Cluster;
}


paillier_ciphertext_t*** protocol::extract_sample(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int* Sample_Cnt){
	printf("Start extract_sample\n");
	int partition = 10;
	int* NodeCnt = (int*)malloc(sizeof(int)*NumNode);

	for(int i = 0; i < NumNode ; i++){
		NodeCnt[i] = node[i].NumData/partition;
		//printf("NodeCnt[%d]  = %d / %d = %d\n", i, node[i].NumData, partition, NodeCnt[i]);
		*Sample_Cnt = (*Sample_Cnt)+NodeCnt[i];
	}
	printf("Sample_Cnt : %d\n", *Sample_Cnt);

	//Sample init
	paillier_ciphertext_t*** Sample = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*(*Sample_Cnt));
	for(int i = 0 ; i < (*Sample_Cnt) ; i++){
		Sample[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j++){
			Sample[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(Sample[i][j]->c);
		}
	}

	int Sample_idx = 0;
	int data_idx = 0;
	for(int i = 0 ; i < NumNode ; i++){
		data_idx = 0;
		while( NodeCnt[i] > 0 ){
			//gmp_printf("%d ", node[i].indata_id[data_idx]);
			for(int j = 0 ; j < dim ; j++ ){
				Sample[Sample_idx][j] = ciper[node[i].indata_id[data_idx]][j];
			}
			data_idx++;
			Sample_idx++;
			NodeCnt[i]--;
		}
		//printf("\n");
	}
	/*
	printf("Sample Data\n");
	for(int i = 0 ; i < (*Sample_Cnt) ; i++ ){
		for(int j = 0 ; j < dim ; j++){
			gmp_printf("%Zd ", paillier_dec(0, pub, prv, Sample[i][j]));
		}
		printf("\n");
	}
	*/
	free(NodeCnt);

	printf("End extract_sample\n");
	return Sample;
}

paillier_ciphertext_t*** protocol::Clustering_m_Grid(paillier_ciphertext_t*** ciper, int NumData, int k, int b, paillier_ciphertext_t*** former_Center){
	printf("Clustering_m_Grid start\n");
	
	paillier_ciphertext_t*** origin_data				= (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	for(int i = 0 ; i < NumData ; i++ ){
		origin_data[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for(int j = 0 ; j < dim ; j++){
			origin_data[i][j] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(origin_data[i][j]->c);
			origin_data[i][j] = ciper[i][j];
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
		/*
		printf("former CNT\n");
		for(int i = 0 ; i < k; i++){
			gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, formerSumCntCluster[i]));
		}
		printf("\n");

		printf("former Center\n");
		for(int i = 0 ; i < k; i++){
			for(int j = 0 ; j < dim ; j++){
				gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, formerSumCluster[i][j]));
			}
			printf("\n");
		}
		printf("\n");
		*/
		
		/*
		printf("New CNT : ");
		for(int i = 0 ; i < k; i++){
			gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, NewSumCntCluster[i]));
		}
		printf("\n");

		printf("New Center\n");
		for(int i = 0 ; i < k; i++){
			for(int j = 0 ; j < dim ; j++){
				gmp_printf(" %Zd ", paillier_dec(0, pubkey, prvkey, NewSumCluster[i][j]));
			}
			printf("\n");
		}
		printf("\n");
		*/
		
		for(int i = 0 ; i < NumData ; i++){
						
			for( int j = 0 ; j < k ; j++){
				Distance_center_data[i][j] = SSEDm(ciper[i], formerSumCluster[j], dim);
				//Distance_center_data[i][j] = DP_SSED(ciper[i], formerSumCluster[j], dim);
				//Distance_center_data[i][j] = OP_SSED(ciper[i], formerSumCluster[j], formerSumCntCluster, i, k);
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

paillier_ciphertext_t* protocol::boundary_dist(boundary node, paillier_ciphertext_t** Center, int NumNode){
	//printf("boundary_dist start\n");
	
	
	int vertex_cnt = pow(2, dim);
	paillier_ciphertext_t*** vertex = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*vertex_cnt);
	paillier_ciphertext_t** vertex_distance = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*vertex_cnt);

	for( int i = 0 ; i < vertex_cnt ; i++ ){
		vertex_distance[i] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
		mpz_init(vertex_distance[i]->c);
		vertex[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		for( int m = 0 ; m < dim ; m++ ){
			vertex[i][m] = (paillier_ciphertext_t*)malloc(sizeof(paillier_ciphertext_t));
			mpz_init(vertex[i][m]->c);
		}
	}

	int* c = (int *)malloc(sizeof(int)*dim);
	for( int i = 0 ; i < dim ; i++ ){
		c[i] = 0;
	}

	bool a = true;
	
	for( int x = 0 ; x < dim ; x++){
		for( int m = 0 ; m < vertex_cnt ; m++ ){
			if(	c[x] == vertex_cnt/pow(2,(x+1)) ){
				a = !a;
				c[x] = 0;
			}
			c[x]++;
				if( a ){
				vertex[m][x] = node.LL[x];
			}else{
				vertex[m][x] = node.RR[x];
			}
		}
	}
	/*
	printf("node boundary\n");
	for(int i = 0 ; i < dim ; i++){
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, node.LL[i]));
	}
	printf("\n");
	for(int i = 0 ; i < dim ; i++){
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, node.RR[i]));
	}
	printf("\n");
	printf("vertex\n");
	for(int i = 0 ; i < vertex_cnt ; i++){
		for(int j = 0 ; j < dim ; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, vertex[i][j]));
		}
		printf("\n");
	}
	*/

	/*
	printf("vertex\n");
	for(int i = 0 ; i < vertex_cnt ; i++){
		for(int j = 0 ; j < dim ; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, vertex[i][j]));
		}
		printf("\n");
	}
	for(int i = 0 ; i < dim ; i++){
		gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, Center[i]));
	}
	printf("\n");
	printf("vertex_distance\n");
	*/
	for(int i = 0 ; i < vertex_cnt ; i++ ){
		vertex_distance[i] = SSEDm(Center, vertex[i], dim);
		//gmp_printf("DIST %Zd\n", paillier_dec(0, pubkey, prvkey, vertex_distance[i]));
	}
	
	int temp = MAXn(vertex_distance, vertex_cnt);

	//gmp_printf("max : %Zd\n", paillier_dec(0, pubkey, prvkey, vertex_distance[temp]));

	//printf("boundary_dist end\n");
	return vertex_distance[temp];
}


paillier_ciphertext_t*** protocol::Clustering_Grid_preprocessing(paillier_ciphertext_t*** ciper, boundary* node, int NumNode, int NumData, int k){
	printf("Start Clustering_Grid_preprocessing\n");
	int i = 0, j = 0, m = 0;
											
	paillier_ciphertext_t*** formerCenter = PreClustering(ciper, node, NumNode, NumData, k);
	

	for( i = 0 ; i < k ; i++ ){
		for( j = 0 ; j < dim ; j++ ){
			gmp_printf("%Zd ", paillier_dec(0, pubkey, prvkey, formerCenter[i][j]));
		}
		printf("\n");
	}

	paillier_ciphertext_t*** NewCluster = Clustering_m_Grid(ciper, NumData, k, 100, formerCenter);
	

	printf("End Clustering_Grid\n");
	return NewCluster;
}