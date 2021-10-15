#ifndef BLOCK_ASSEMBLER_H
#define BLOCK_ASSEMBLER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <bits/stdc++.h> 
#include "stringfuncs.h"
#include "structdefs.h"
#include "read_processor.h"

class block_assembler : public read_processor{

public:
std::vector<contig> s; //stores all the blocks
std::vector<contig> u,v; //u and v are contiguous blocks
std::vector<long long int> vecsizes; //sizes of all the blocks
int no_of_rounds;  //no of rounds for the assembly
std::multimap <std::string,int> pmap; //maps
int klen; //current overlap of map
long long int maxlenU,maxlenV;
double th = 0.2;
	void no_of_contigs();

	void push_block(std::vector<contig> temp);

	void push_block_size(long long int size);

	void set_no_of_rounds(int NO_OF_CLUSTERS);

	void setmap(std::string seq, int i);

	void getmap( int curr_klen);

	void infomap();

	void freemap();

	void find(int i,std::vector<int> &markU,std::vector<int> &markV);

	void run_match(std::vector<int> &markU,std::vector<int> &markV);

	void match_kernel();

	void getU(int iter);

	void getV(int iter);

	void block_assembly();

	long long int get_max_length_u();

	long long int get_max_length_v();

}; 

#endif // BLOCK_ASSEMBLER_H