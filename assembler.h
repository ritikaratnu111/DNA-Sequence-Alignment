#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <bits/stdc++.h> 
#include "stringfuncs.h"
#include "structdefs.h"
#include "read_processor.h"

class assembler : public read_processor{
protected:
	std::multimap <std::string,int> pmap;
	int klen;
	
	struct jnInfo{ 
	   int idx;
	   int len;
	   float nOv; 
	}; 

	std::vector <jnInfo> cjn; 
	int maxLen;
	float maxNOv;
	float thLen = 0.7;
	float thNov = 0.8;

public:

	void getmap(int curr_klen);

	void setmap(std::string seq, int i);

	void infomap();

	void freemap();	

	void displaymap();

	void run_match();

	void find (std::string seq, std::list <node> & sList,int i, std::vector <int> &mark, std::vector<contig> &vnew);


}; 

#endif // ASSEMBLER_H