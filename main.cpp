#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <bits/stdc++.h> 
#include <ctime> 
#include <math.h>   

#include "stringfuncs.h"
#include "structdefs.h"
#include "sam_filter.h"
#include "read_processor.h"
#include "assembler.h"
#include "block_assembler.h"
#include "fasta_converter.h"
#include "sam_filter.h"

int main(){

	int READ_LENGTH = 150;
	int NO_OF_CLUSTERS = 2048;
	long int genome_length;

	// double rover = r * 0.1;
	long long int total_reads,unext_reads,max_extension_length;
	float percent_reads_unext;

	clock_t begin,end;

	sam_filter s;
	assembler m1;
	block_assembler b1;
	fasta_converter g;

	std::string sam_original = "/home/ritika/thesis_project/assembly/reference/BCep_ef.sam" ;
	std::string isam = "/home/ritika/thesis_project/assembly/reference/BCep_ef_no_repeats.sam" ;
	std::string oname = "/home/ritika/thesis_project/assembly/assembly_output/n2048_ef.fa";


	// filtered sam with no repeats has already been created
	// s.filter_sam(sam_original);
	// s.store_new_sam(isam);

	m1.header_count(0);
	m1.find_min_max_pos(isam);
	m1.set_radius(NO_OF_CLUSTERS);
	
	long long int size = 0;

	std::cout << "Cluster Assembly::\n";

	for (int i = 0; i < NO_OF_CLUSTERS; i++){

		m1.get_data_in_radius(isam,i);
		total_reads = m1.noOfReads ;
		begin = clock();
		std::cout << "data collected\n";
		for(int kmerLen = READ_LENGTH - 1 ; kmerLen >= 13; kmerLen--){
			m1.getmap(kmerLen);
			m1.run_match();
			m1.freemap();
			//unext_reads = m1.no_of_unextended_reads(READ_LENGTH);
			//percent_reads_unext = ( float(unext_reads) / total_reads ) * 100;
			std::cout <<"kmerlen: "<< kmerLen << "\t" << m1.v.size() << "\n";
		}
		end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		std::cout << "Time: " << elapsed_secs << "\n";
		size += m1.v.size(); 
		b1.push_block(m1.v);
		b1.push_block_size(size);
		std::cout << i << "\t";
		b1.no_of_contigs();
		// m1.write_fasta("myassem.fa");

		m1.freereads();	
	}


	std::cout << "Block Assembly::\n";

	int no_of_rounds = log2(NO_OF_CLUSTERS);

	// std::cout << "Block assembly disabled\n";
	std::cout << "No of rounds: " << no_of_rounds << "\n";

	for(int round = 0; round < no_of_rounds; round++){
			std::cout << "Round " << round << "\t";
			b1.block_assembly();					
	}

 	g.push_contigs(b1.s);
	g.getv();
	// g.set_olen(150);
	// g.remove_redundancy();

	g.write_to_file(oname);
	g.vec_size();

	return 0;
}

