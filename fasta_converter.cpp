#include "fasta_converter.h"

void fasta_converter::set_olen(int ov){
	olen = ov;
}

void fasta_converter::vec_size(){
	std::cout << "No of contigs:" << v.size() << "\n";	
}

void fasta_converter :: push_contigs(std::vector<contig> temp){

	final_contigs.insert(std::end(final_contigs), std::begin(temp), std::end(temp));
}

void fasta_converter::getv(){

	for(int i = 0; i<final_contigs.size(); i++){
		v.push_back(final_contigs[i].seq);
	}		
}

void fasta_converter::remove_redundancy(){
	
	std::vector<std::string> tempvec;
	int olen = 5;
	std::string prefix;
	int maxlen, maxIdx;

	sort(v.begin(), v.end());  // sort the vector	
	std::vector<int> flag(v.size(),0);

	int i = 0;
	int j;

	while(i < v.size()){
		prefix = v[i].substr(0,olen-1);
		maxlen = v[i].length();
		maxIdx = i;
		j = i+1;
		
		while( v[j].substr(0,olen-1) == prefix ){
			if(v[j].length() > maxlen){
				maxlen = v[j].length();
				maxIdx = j;
			}				
			j++;
			if(j == v.size() - 1)
				break;					
		}
		flag[maxIdx] = 1;
		i = j;
		if(i == v.size() - 1)
			break;	
	}


	for(i = 0; i< v.size(); i++){
		if(flag[i] == 1)
			tempvec.push_back(v[i]);
	}

	v.clear();
	v = tempvec;
	tempvec.clear();

}

void fasta_converter::write_to_file(std::string filename){
	std::ofstream outfile;
	outfile.open(filename.c_str(), std::ios_base::app);
	

	for(long long int i = 0; i < v.size(); i++){
		outfile << ">contig" << i << " len=" << v[i].length() << "\n";
		outfile << v[i] << "\n";
	}

	outfile.close();
}

