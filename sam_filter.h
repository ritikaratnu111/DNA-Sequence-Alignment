#ifndef SAM_FILTER_H
#define SAM_FILTER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <bits/stdc++.h> 
#include <ctime> 
#include <math.h>   

#include "stringfuncs.h"

class sam_filter{
 
public:

struct record{
	std::string QNAME;
	int FLAGS;
	std::string RNAME;
	int POS;
	int MAPQ;
	std::string CIGAR;
	std::string RNEXT;
	int PNEXT;
	int TLEN;
	std::string SEQ;
};
std::vector<record> v;

friend std::istream& operator>>(std::istream& is, sam_filter::record& r);

void filter_sam(std::string filename);

void store_new_sam(std::string filename);

};

#endif // SAM_FILTER_H