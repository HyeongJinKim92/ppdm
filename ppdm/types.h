#ifndef _TYPES_H_
#define	 _TYPES_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <sstream>
#include <fstream>

enum QUERY
{
	RANGE,
	TOPK,
	KNN,
	CLASSIFICATION,
	KMEANS,
	ASSOCATION,
	TEST
};

enum RANGE_TYPE
{
	RANGE_B,
	RANGE_I,
	RANGE_GI,
	RANGE_PB,
	RANGE_PGI,
	RANGE_PAI
	//
};

enum TOPK_TYPE
{
	TOPK_B,
	TOPK_I,
	TOPK_GI,
	TOPK_PB,
	TOPK_PGI,
	TOPK_PAI
	//
};

enum KNN_TYPE
{
	KNN_B,
	KNN_I,
	KNN_GI,
	KNN_PB,
	KNN_PGI,
	KNN_PAI
	//
};

enum CLASSIFICATION_TYPE
{
	CLASSIFICATION_B,
	CLASSIFICATION_I,
	CLASSIFICATION_GI,
	CLASSIFICATION_PB,
	CLASSIFICATION_PGI,
	CLASSIFICATION_PAI
	//
};

enum KMEANS_TYPE
{
	KMEANS_B,
	KMEANS_I,
	KMEANS_GI,
	KMEANS_PB,
	KMEANS_PGI,
	KMEANS_PAI
	//
};

enum ASSOCIATION_TYPE
{
	//
};



struct parsed_query
{
	QUERY query;
	RANGE_TYPE rquery_t;
	TOPK_TYPE tquery_t;
	KNN_TYPE kquery_t;
	CLASSIFICATION_TYPE cquery_t;
	KMEANS_TYPE mquery_t;
	int tree_level;
	int data_n;
	int dim_n;
	int garbled;
	int thread_n;
	int key_s;
	int bit_s;
	int k;
	int hint;
	int * q;
	int * score;
	int * qLb;
	int * qUb;				
};



#endif