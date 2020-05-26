#ifndef _PPDM_H_
#define	 _PPDM_H_

#include <iostream>
#include <fstream>
#include <gmp.h>
#include <string>
#include <sstream>
#include <thread>
#include <mutex>
#include <chrono>
#include <vector>

#include "protocol.h"
#include "setting.h"
#include "types.h"
#include "setting.h"
#include "util/config.h"
#include "circuit/circuit.h"



using namespace std;

void ppdm_main(char* argv[]);

int range_main(parsed_query*);
int topk_main(parsed_query*);
int knn_main(parsed_query*);
int classification_main(parsed_query*);
int kmeans_main(parsed_query*);


parsed_query * parsing_query(char* argv[]);

void print_query();

#endif