#include "sam_filter.h"

std::istream& operator>>(std::istream& is, sam_filter::record& r){
	is >> r.QNAME >> r.FLAGS >> r.RNAME >> r.POS >> r.MAPQ >> r.CIGAR >> r.RNEXT >> r.PNEXT >> r.TLEN >> r.SEQ ;
	return is;
}

void sam_filter :: filter_sam(std::string filename){
	std::ifstream infile(filename.c_str()); 
	record temp;
	std::string str,read;
	std::vector<std::string> tempseq;

	if(infile){
		int no_of_headers = 6;

		while(no_of_headers > 0){ //ignore headers
			getline(infile,str);
			no_of_headers--;
		}	

		int count = 0;
		while(getline(infile,str) ){
			count++;
			std::cout << count << "\n";
			std::istringstream iss (str);
			iss >> temp; 
			if (temp.SEQ != "*"){
		
				if(temp.FLAGS == 99){       // 1st in pair, mapped as it is
					read = temp.SEQ;
				}
				else if(temp.FLAGS == 83){  // 1st in pair, RCed
					read = rev_complement(temp.SEQ);
				}
				else if(temp.FLAGS == 163){ // 2nd in pair, mapped as it is
					read = temp.SEQ;
				}
				else if(temp.FLAGS == 147){ // 2nd in pair, RCed
					read = rev_complement(temp.SEQ);
				}
				else{                 
					read = temp.SEQ;
				}	
			}
			if( std::find(tempseq.begin(), tempseq.end(),read) == tempseq.end() ){   //only insert unique entries
				v.push_back(temp);
				tempseq.push_back(read);
			}	
		}
	}
	else{
		std::cout <<  "Unable to open file\n";
	}
}

void sam_filter :: store_new_sam(std::string filename){
	std::ofstream outfile;
	outfile.open(filename.c_str(), std::ios_base::app);
	record r;

	for(long long int i = 0; i < v.size(); i++){

		r = v[i];
        outfile << r.QNAME << "\t" << r.FLAGS << "\t" << r.RNAME << "\t" << r.POS << "\t" << r.MAPQ << "\t" << r.CIGAR << "\t" << r.RNEXT << "\t" << r.PNEXT << "\t" << r.TLEN << "\t" << r.SEQ << "\n" ;
	}

	outfile.close();		
}
