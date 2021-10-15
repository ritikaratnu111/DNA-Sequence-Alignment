#ifndef FASTA_CONVERTER_H
#define FASTA_CONVERTER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <bits/stdc++.h> 
#include "stringfuncs.h"
#include "structdefs.h"
#include "read_processor.h"
#include "block_assembler.h"

class fasta_converter{
public:

	std::vector<contig> final_contigs;	
	std::vector<std::string> v;
	std::vector<std::string> vprefix;
	int olen ;

	void set_olen(int ov);

	void vec_size();

	void push_contigs(std::vector<contig> temp);	

	void getv();

	void remove_redundancy();

	void write_to_file(std::string filename);
};

#endif // FASTA_CONVERTER_H
