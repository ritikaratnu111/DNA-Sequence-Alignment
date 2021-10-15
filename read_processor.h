#ifndef READ_PROCESSOR_H
#define READ_PROCESSOR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <bits/stdc++.h> 
#include "stringfuncs.h"
#include "structdefs.h"


class read_processor{

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

	std::vector<contig> v;
	std::vector<std::string> junctionPts;  //dont care

	int no_of_headers;
	int noOfReads = 0;

	long int pmin,pmax;
	double radius;
	double rover;

	friend std::istream& operator>>(std::istream& is, record& r);

	void header_count(int h);

	void find_min_max_pos(std::string filename);

	void set_radius(int NO_OF_CLUSTERS);

	void get_data_in_radius(std::string filename, int i);

	int no_of_unextended_reads(int read_length);

	int max_extension_length();
	
	void freereads();


};

#endif // READ_PROCESSOR_H