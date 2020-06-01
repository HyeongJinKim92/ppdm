#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <sstream>
#include <fstream>
#include <map>

#include "protocol.h"
#include "setting.h"

#define TOKEN " "
using namespace std;

boundary* setting:: KdInfo_read(char* kdFilename, int dim ,int* NumNode){	
	boundary* node;

	paillier_plaintext_t * plain;
	FILE *pFile = NULL;
	char *token;
	char strTemp[1000000];
	char *pStr = NULL;
	char file_name[128]; 
	int i=0;
	int j=0;

	pFile = fopen(kdFilename, "r" );

	printf("kd file name : %s\n",kdFilename);
	
	//전체 노드 수 읽는
	if( pFile != NULL ){
		pStr = fgets( strTemp, sizeof(strTemp), pFile );		
		*NumNode=atoi(pStr);
	}
	
	printf("NumNode : %d\n",*NumNode);

	node = (boundary*)malloc(sizeof(boundary)*(*NumNode));
	for(i=0; i<*NumNode; i++) {
		node[i].LL =  (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
		node[i].RR =  (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	}
	
	// node setting
	i=0;
	if( pFile != NULL ){   	//파일에서 kd-tree 정보 읽기		
		 while( !feof( pFile ) ){				
            pStr = fgets( strTemp, sizeof(strTemp), pFile );
			token=strtok(pStr ,TOKEN);
			j=0;
			while(token != NULL && i< *NumNode ){				
				if(j==0){	//노드 id
					node[i].node_id=atoi(token);
					//printf("\nID : %d\n", node[i].node_id);
				}
				else if(j>0 && j<=dim){ //LL
					plain = paillier_plaintext_from_ui(atoi(token));
					node[i].LL[j-1]=paillier_enc(0, pub, plain, paillier_get_rand_devurandom);
					//printf("LL : %s\n", token);
				}else if(j>dim && j<=dim*2){ //RR
					plain = paillier_plaintext_from_ui(atoi(token));
					node[i].RR[j-dim-1]=paillier_enc(0, pub, plain, paillier_get_rand_devurandom);
					//printf("RR : %s\n", token);
				}else if(j>dim*2 && j<=dim*2+1){ // FanOut 
					node[i].NumData=atoi(token);
					node[i].indata_id=(int*)malloc(sizeof(int)*atoi(token));	
					//printf("FanOut : %d\n", node[i].NumData);
				}else if(j>dim*2+1 && j<=dim*2+1+node[i].NumData){
					//printf("i : %d , j : %d , indata : %d\n",i,j,atoi(token));
					node[i].indata_id[j-(2*dim)-2]=atoi(token);
					//printf("i : %d , j : %d , indata : %d\n",i,j,atoi(token));
					//printf("mcmcmc\n");
					//printf("%d ", node[i].indata_id[j-(2*dim)-2]);
					//printf(" j : %d \n",j);
				}				
				token = strtok(NULL, TOKEN);
				j++;
			 }
			 i++;
		}
    } //파일 읽기 끝 
	fclose(pFile);

	return node; 
}

paillier_ciphertext_t*** setting::InputData_read(char* inputFilename,int dim, int NumData){
	//input 데이터 파일 읽기
	FILE *pFile = NULL;
	char *token;
	char strTemp[255];
	char *pStr = NULL;

	pFile = fopen(inputFilename, "r" );

	printf("input file name : %s\n",inputFilename);

	int i=0;
	int j=0;

	paillier_plaintext_t * plain;

	paillier_ciphertext_t*** cipher = (paillier_ciphertext_t***)malloc(sizeof(paillier_ciphertext_t**)*NumData);
	for( i = 0 ; i < NumData ; i++ ){
		cipher[i] = (paillier_ciphertext_t**)malloc(sizeof(paillier_ciphertext_t*)*dim);
	}
	i=0;
	if( pFile != NULL ){   	
		 while( !feof( pFile ) ){	
			j=0;
			pStr = fgets( strTemp, sizeof(strTemp), pFile );
			token=strtok(pStr ,TOKEN);
			while(token != NULL && token !=" " ){
				//printf("%s\n",token);
				plain = paillier_plaintext_from_ui(atoi(token));
				cipher[i][j] = paillier_enc(0, pub, plain, paillier_get_rand_devurandom);
				token = strtok(NULL, TOKEN);
				j++;
			 }
			 i++;
		}
	} //파일 읽기 끝 
	
	fclose(pFile);
	return cipher;
	
}

void setting::TimeResult_write_cipher(float time,paillier_ciphertext_t*** result,int result_num,int NumData,char* app, int bitsize,int dim){
	FILE *fp_output;
	char file_name[128]; 
	sprintf(file_name, "output/%s/d%d_L%d.txt",app ,NumData,bitsize); 
	fp_output=fopen(file_name,"a");

	fprintf(fp_output,"time : %f\n", time);	 // 시간 기록

	//result file write
	for(int i=0;i<result_num;i++){
		for(int j=0; j<dim; j++) {				
				paillier_plaintext_t* output_for_file=paillier_dec(0, pub, prv, result[i][j]);
				gmp_fprintf(fp_output,"%Zd  ",output_for_file);
			}
		fprintf(fp_output,"\n");
	}	
	fprintf(fp_output,"\n\n");
	fclose(fp_output);	
}

void setting::TimeResult_write_map(char* query, int** result, int result_num, char* app, protocol proto)
{

	std::time_t rawtime;
	std::tm* timeinfo;
	char buffer [80];
	
	std::time(&rawtime);
	timeinfo = std::localtime(&rawtime);

	std::strftime(buffer,80,"%Y-%m-%d-%H-%M-%S",timeinfo);
	std::puts(buffer);


	FILE *fp_output;
	char file_name[128]; 
	sprintf(file_name, "output/%s/d%d_m%d_L%d_k%d_h%d_K%d.txt", app, proto.NumData, proto.dim, proto.size, proto.k, proto.tree_level, proto.modul); 
	fp_output=fopen(file_name,"a");
	cout.precision(6);
  

	for(auto it = proto.time_variable.begin() ; it != proto.time_variable.end() ; it++ ){
		cout.width(20);
		cout <<  it->first << " ";
	}
	cout << endl;
	for(auto it = proto.time_variable.begin() ; it != proto.time_variable.end() ; it++ ){
		cout.width(20);
		cout <<  it->second << " ";
	}
	cout << endl;
	for(auto it = proto.time_variable.begin() ; it != proto.time_variable.end() ; it++ ){
		cout.width(20);
		cout <<  (it->second/proto.total_time)*100 << " ";
	}
	cout << endl;

	fprintf(fp_output,"time : %s\n", buffer);
	fprintf(fp_output,"query : %s\n", query);
	fprintf(fp_output,"total time : %.6f\n", proto.total_time);
	for(auto it = proto.time_variable.begin() ; it != proto.time_variable.end() ; it++ ){
		fprintf(fp_output,"%20s", it->first);
	}
	fprintf(fp_output,"\n");
	for(auto it = proto.time_variable.begin() ; it != proto.time_variable.end() ; it++ ){
		fprintf(fp_output,"%20.6f", it->second);
	}
	fprintf(fp_output,"\n");
	for(auto it = proto.time_variable.begin() ; it != proto.time_variable.end() ; it++ ){
		fprintf(fp_output,"%20.6f", (it->second/proto.total_time)*100);
	}
	fprintf(fp_output,"\n");

	fprintf(fp_output,"# retrieved node : %d\n", proto.totalNumOfRetrievedNodes); // 탐색한 총 노드의 개수 기록

	//result file write
	for( int i = 0 ; i < result_num ; i++ ){
		fprintf(fp_output,"Result(%d) :\t", i);
		for( int j = 0 ; j < proto.dim ; j++ ) {
				fprintf(fp_output,"%d\t", result[i][j]);
			}
		fprintf(fp_output,"\n");
	}	
	fprintf(fp_output,"\n\n");
	fclose(fp_output);	
}


void setting::TimeResult_write_int(char* query, int** result, int result_num, char* app, protocol proto){
	//protocol proto;
	FILE *fp_output;
	char file_name[128]; 
	sprintf(file_name, "output/%s/d%d_m%d_L%d_k%d_h%d_K%d.txt", app, proto.NumData, proto.dim, proto.size, proto.k, proto.tree_level, proto.modul); 
	fp_output=fopen(file_name,"a");

	printf("query_Processing_time : %.0f  (%.0f%%)\n", proto.query_Processing_time, (proto.query_Processing_time/proto.total_time)*100);
	printf("node_SBD : %.0f  (%.0f%%)\n", proto.node_Processing_time, (proto.node_Processing_time/proto.total_time)*100);
	printf("node_SRO : %.0f  (%.0f%%) \t node_expansion(verify) : %.0f (%.0f%%)\n", proto.node_SRO_time, (proto.node_SRO_time/proto.total_time)*100, proto.node_expansion_time, (proto.node_expansion_time/proto.total_time)*100);
	printf("node_retrieval : %.0f  (%.0f%%) \t node_retrieval(verify) : %.0f (%.0f%%)\n", proto.data_extract_first_time, (proto.data_extract_first_time/proto.total_time)*100, proto.data_extract_second_time, (proto.data_extract_second_time/proto.total_time)*100);

	printf("data_Processing : %.0f (%.0f%%) \t data_SSED_SBD : %.0f  (%.0f%%) \t data_SBOR : %.0f (%.0f%%)\n", proto.data_Processing_time, (proto.data_Processing_time/proto.total_time)*100, proto.data_SSED_SBD_time, (proto.data_SSED_SBD_time/proto.total_time)*100, proto.data_SBOR_time, (proto.data_SBOR_time/proto.total_time)*100);
	printf("sMINn : %.0f  (%.0f%%) \t sMINn (verify): %.0f (%.0f%%)\n", proto.sMINn_first_time, (proto.sMINn_first_time/proto.total_time)*100, proto.sMINn_second_time, (proto.sMINn_second_time/proto.total_time)*100);
	printf("data_SPE(range) : %.0f  (%.0f%%) \n", proto.data_SPE_time, (proto.data_SPE_time/proto.total_time)*100);


	fprintf(fp_output,"query : %s\n", query);
	fprintf(fp_output,"total time : %.0f\n", proto.total_time);
	fprintf(fp_output,"query_Processing_time : %.0f  (%.0f%%)\n", proto.query_Processing_time, (proto.query_Processing_time/proto.total_time)*100);
	fprintf(fp_output,"node_SBD : %.0f  (%.0f%%)\n", proto.node_Processing_time, (proto.node_Processing_time/proto.total_time)*100);
	fprintf(fp_output,"node_SRO : %.0f  (%.0f%%) \t node_expansion(verify) : %.0f (%.0f%%)\n", proto.node_SRO_time, (proto.node_SRO_time/proto.total_time)*100, proto.node_expansion_time, (proto.node_expansion_time/proto.total_time)*100);
	fprintf(fp_output,"node_retrieval : %.0f  (%.0f%%) \t node_retrieval(verify) : %.0f (%.0f%%)\n",proto.data_extract_first_time, (proto.data_extract_first_time/proto.total_time)*100, proto.data_extract_second_time, (proto.data_extract_second_time/proto.total_time)*100);
	printf("data_Processing : %.0f (%.0f%%) \t data_SSED_SBD : %.0f  (%.0f%%) \t data_SBOR : %.0f (%.0f%%)\n", proto.data_Processing_time, (proto.data_Processing_time/proto.total_time)*100, proto.data_SSED_SBD_time, (proto.data_SSED_SBD_time/proto.total_time)*100, proto.data_SBOR_time, (proto.data_SBOR_time/proto.total_time)*100);
	fprintf(fp_output,"sMINn : %.0f  (%.0f%%) \t sMINn (verify): %.0f (%.0f%%)\n", proto.sMINn_first_time, (proto.sMINn_first_time/proto.total_time)*100, proto.sMINn_second_time, (proto.sMINn_second_time/proto.total_time)*100);
	fprintf(fp_output,"data_SPE(range) : %.0f  (%.0f%%) \n", proto.data_SPE_time, (proto.data_SPE_time/proto.total_time)*100);

	fprintf(fp_output,"# retrieved node : %d\n", proto.totalNumOfRetrievedNodes); // 탐색한 총 노드의 개수 기록

	//result file write
	for(int i=0;i<result_num;i++){
		for(int j=0; j<proto.dim; j++) {				
				fprintf(fp_output,"%d = %d\t", i, result[i][j]);
			}
		fprintf(fp_output,"\n");
	}	
	fprintf(fp_output,"\n\n");
	fclose(fp_output);	
}

void setting::TimeResult_write_int_forRange(float time,int** result,int k,int result_num,int NumData,char* app, int bitsize,int dim, int tree_level, int KEYSize, protocol proto){

	FILE *fp_output;
	char file_name[128]; 
	sprintf(file_name, "output/%s/d%d_m%d_L%d_k%d_h%d_M%d.txt",app ,NumData,dim,bitsize,k,tree_level, KEYSize); 
	fp_output=fopen(file_name,"a");

	printf("node_SBD : %.0f  (%.0f%%)\n", proto.node_Processing_time, (proto.node_Processing_time/proto.total_time)*100);
	printf("node_SRO : %.0f  (%.0f%%) \t node_expansion(verify) : %.0f (%.0f%%)\n", proto.node_SRO_time, (proto.node_SRO_time/proto.total_time)*100, proto.node_expansion_time, (proto.node_expansion_time/proto.total_time)*100);
	printf("node_retrieval : %.0f  (%.0f%%) \t node_retrieval(verify) : %.0f (%.0f%%)\n", proto.data_extract_first_time, (proto.data_extract_first_time/proto.total_time)*100, proto.data_extract_second_time, (proto.data_extract_second_time/proto.total_time)*100);
	printf("data_SSED_SBD : %.0f  (%.0f%%) \t data_SBOR : %.0f (%.0f%%)\n", proto.data_SSED_SBD_time, (proto.data_SSED_SBD_time/proto.total_time)*100, proto.data_SBOR_time, (proto.data_SBOR_time/proto.total_time)*100);
	printf("sMINn : %.0f  (%.0f%%) \t sMINn (verify): %.0f (%.0f%%)\n", proto.sMINn_first_time, (proto.sMINn_first_time/proto.total_time)*100, proto.sMINn_second_time, (proto.sMINn_second_time/proto.total_time)*100);
	printf("data_SPE(range) : %.0f  (%.0f%%) \n", proto.data_SPE_time, (proto.data_SPE_time/proto.total_time)*100);

	fprintf(fp_output,"total time : %.0f\n", proto.total_time);
	fprintf(fp_output,"total time : %.0f\n", proto.total_time);
	fprintf(fp_output,"node_SBD : %.0f  (%.0f%%)\n", proto.node_Processing_time, (proto.node_Processing_time/proto.total_time)*100);
	fprintf(fp_output,"node_SRO : %.0f  (%.0f%%) \t node_expansion(verify) : %.0f (%.0f%%)\n", proto.node_SRO_time, (proto.node_SRO_time/proto.total_time)*100, proto.node_expansion_time, (proto.node_expansion_time/proto.total_time)*100);
	fprintf(fp_output,"node_retrieval : %.0f  (%.0f%%) \t node_retrieval(verify) : %.0f (%.0f%%)\n",proto.data_extract_first_time, (proto.data_extract_first_time/proto.total_time)*100, proto.data_extract_second_time, (proto.data_extract_second_time/proto.total_time)*100);
	fprintf(fp_output,"data_SSED_SBD : %.0f  (%.0f%%) \t data_SBOR : %.0f (%.0f%%)\n", proto.data_SSED_SBD_time, (proto.data_SSED_SBD_time/proto.total_time)*100, proto.data_SBOR_time, (proto.data_SBOR_time/proto.total_time)*100);
	fprintf(fp_output,"sMINn : %.0f  (%.0f%%) \t sMINn (verify): %.0f (%.0f%%)\n", proto.sMINn_first_time, (proto.sMINn_first_time/proto.total_time)*100, proto.sMINn_second_time, (proto.sMINn_second_time/proto.total_time)*100);
	fprintf(fp_output,"data_SPE(range) : %.0f  (%.0f%%) \n", proto.data_SPE_time, (proto.data_SPE_time/proto.total_time)*100);

	fprintf(fp_output,"# retrieved node : %d\n", proto.totalNumOfRetrievedNodes); // 탐색한 총 노드의 개수 기록

	//printf("result num : %d\n",result_num);

	//result file write
	for(int i=0;i<result_num;i++){
		for(int j=0; j<dim; j++) {				
				fprintf(fp_output,"%d\t", result[i][j]);
			}
		fprintf(fp_output,"\n");
	}	
	fprintf(fp_output,"\n\n");
	fclose(fp_output);

}