#ifndef _SETTING_H_
#define	 _SETTING_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <sstream>
#include <fstream>

#include "protocol.h"

class setting
{
	public :
		boundary* KdInfo_read(char* kdFilename, int dim ,int* NumNode);
		paillier_ciphertext_t*** InputData_read(char* inputFilename,int dim, int NumData);
		void TimeResult_write_ciper(float time,paillier_ciphertext_t*** result,int result_num,int NumData,char* app, int bitsize,int dim);
		void TimeResult_write_int(char* query, int** result, int result_num, char* app, protocol proto);
		void TimeResult_write_map(char* query, int** result, int result_num, char* app, protocol proto);
		void TimeResult_write_int_forRange(float time,int** result,int k,int result_num,int NumData,char* app, int bitsize,int dim, int tree_level,int KEYSize , protocol proto);
};
#endif