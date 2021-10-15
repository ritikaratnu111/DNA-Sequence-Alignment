#include "read_processor.h"


std::istream& operator>>(std::istream& is, read_processor::record& r){
	is >> r.QNAME >> r.FLAGS >> r.RNAME >> r.POS >> r.MAPQ >> r.CIGAR >> r.RNEXT >> r.PNEXT >> r.TLEN >> r.SEQ ;
	return is;
}

void read_processor::header_count(int h){
	no_of_headers = h;
}

void read_processor::find_min_max_pos(std::string filename){
	std::ifstream infile(filename.c_str());
	std::string str;
	record temp;
	long int min = 0;
	long int max = 0;
	int no_of_records = 0;
	int negPOS = 0;
	if (infile){

		while(no_of_headers > 0){ //ignore headers
			getline(infile,str);
			no_of_headers--;
		}	
		while(getline(infile,str)){
			std::istringstream iss (str);
			iss >> temp;

			if(temp.POS < min)
				min = temp.POS ;
			else if(temp.POS > max)
				max = temp.POS ;

			no_of_records++;

			if(temp.POS < 0)
				negPOS++;

		}
		// std::cout << "\nMin: " << min << "\nMax: " << max << "\nNo of records: "<< no_of_records << "\nNo of records with negPOS: " << negPOS << "\n";
		pmin = min;
		pmax = max;
		std::cout << "Min: " << pmin << "\tMax: "<< pmax << "\n";
	}
}

void read_processor::set_radius(int NO_OF_CLUSTERS){
	std::cout << "Min: " << pmin << "\tMax: "<< pmax << "\tNo of clusters: " << NO_OF_CLUSTERS << "\n";
	double r = (pmax - pmin)/double(NO_OF_CLUSTERS) ;
	radius = r;
	rover = r*0.03;
	std::cout << "Radius: " << radius << "\n";
}


void read_processor::get_data_in_radius(std::string filename, int i){
	std::ifstream infile(filename.c_str()); 
	std::string str;
	record temp;
	std::string SNo;
	std::string name;
	std::string read;
	std::list<node> list1;
	double pos;
	std::vector<std::string> tempseq;

	if(infile){

		while(no_of_headers > 0){ //ignore headers
			getline(infile,str);
			no_of_headers--;
		}	

		int count = 0;
		while(getline(infile,str)){
			count++;
			// std::cout << count << "\n";
			std::istringstream iss (str);
			iss >> temp;  

			if (temp.SEQ != "*"){
				pos = temp.POS;
				if (pos >= (pmin + i*radius - rover) && pos <= ( pmin + (i+1)*radius + rover ) ){
					name = temp.QNAME;
					SNo = name.substr(name.find('.') + 1, name.length());

					if(temp.FLAGS == 99){       // 1st in pair, mapped as it is

						read = temp.SEQ;
						SNo += ".1";
					}

					else if(temp.FLAGS == 83){  // 1st in pair, RCed

						read = rev_complement(temp.SEQ);
						SNo+= ".1'";
					}

					else if(temp.FLAGS == 163){ // 2nd in pair, mapped as it is

						read = temp.SEQ;
						SNo += ".2";
					}

					else if(temp.FLAGS == 147){ // 2nd in pair, RCed

						read = rev_complement(temp.SEQ);
						SNo+= ".2'";
					}
					else{                 

						read = temp.SEQ;
					}	
					
					list1.push_back({SNo,0,0,0,0});
					v.push_back({read,list1});
					tempseq.push_back(read);

				}		
				// else{
				// 	std::cout<< "No!\n";
				// }			
			}


	
			list1.clear();

		}
		
	}


	else{
		std::cout << "Could not open file.\n\n";
	}
	std::cout << "No of reads: " << v.size() << "\n";
	tempseq.clear();
	noOfReads = v.size();
	// std::cout << "After removing duplicates: " << v.size() << "\n";
}


int read_processor::no_of_unextended_reads(int read_length){
	
	int curr_seq_len;
	int unextReads = 0;
	std::string read;

    for (int i =0; i< v.size(); i++){
    	read = v[i].seq ;
    	curr_seq_len = read.length();
    	if(curr_seq_len == read_length )
    		unextReads++;
	}

	return unextReads;

}


int read_processor::max_extension_length(){
	int maxReadLength = 0;
    int maxIndex;
    std::string read;
    for (int i =0; i< v.size(); i++){
    	read = v[i].seq;
    	if (read.length() > maxReadLength){
     		maxReadLength = read.length();    		
    		maxIndex = i;
    	}
	}
	return maxReadLength;
}


void read_processor::freereads(){
	v.clear();
	noOfReads = 0;
	junctionPts.clear();
}

/*
void read_processor::displaylinks(){
	std::list<node> l1;
	std::list <node> :: iterator it,itNxt;

	float ll = 0.8 * (float)max_extension_length();
	float ul = 1.2 * (float)max_extension_length();


	for(int i = 0; i< v.size(); i++){

		// if(v[i].seq.length() > ll && v[i].seq.length() < ul){
			l1 = v[i].l;

		    for(it = l1.begin(); it != l1.end(); ++it){
		    	itNxt = it;
		    	itNxt++;
		    	// if(itNxt != l1.end())
		    	// std::cout << "L\t" << it->SNo << "\t+\t" << itNxt->SNo << "\t+\t" << it->overlap << "M\n";
		    	if(itNxt != l1.end()){	

			    	// if(it == l1.begin())
			    	// 	std::cout << "L\t" << "startOfContig" << "\t+\t" << it->SNo << "\t+\t" << "0" << "M\n";
			    	
			    	std::cout << "L\t" << it->SNo << "\t+\t" << itNxt->SNo << "\t+\t" << it->overlap << "M\n";	    	
			    					
			    }

			    // else
			    // 	std::cout << "L\t" << it->SNo << "\t+\t" << "endOfContig" << "\t+\t" << "0" << "M\n";	    	
		    } 
		    l1.clear();	
		// }
	}		
}



int read_processor::checklist(std::list<node> l1){
	std::list <node> :: iterator it1,it2;
	std::string repeat;

    for(it1 = l1.begin(); it1 != l1.end(); ++it1){
    	repeat = it1->SNo;
    	it2 = it1;
    	++it2;

    	while(it2!= l1.end()){
    		if(repeat == it2->SNo){
    			return 1;
    		}
    		++it2;
    	}

    }

    return 0;
}


void read_processor::updatelinks(){
	std::list<node> l1;
	
	
	contig temp;

	for(auto i = v.begin(); i != v.end(); ++i){

		temp = *i;
		l1 = temp.l;
		int flag = checklist(l1);
		
		if(flag == 1){
			v.erase(i);
			i--;
		}
	}	
}

void read_processor::remove_duplicates(){

	std::string seq;
	std::vector<contig> temp;
	std::vector<int> flag(v.size(),0);
	int i,j;


	for(i = 0; i< v.size(); i++){
		seq  = v[i].seq;
		for( j = i+1; j < v.size(); j++){
			if(v[j].seq == seq){
				flag[j] = 1;
			}
		}
	}
	std::cout << "flags set\n";

	for(int i = 0; i< v.size(); i++){
		if (flag[i] == 0 ){
			temp.push_back(v[i]);
		}
	}
	v.clear();
	v = temp;
	temp.clear();
}
void read_processor::getnodeinfo(){

	std::list<node> l1;
	std::list <node> :: iterator it;
	int s;
	float ov = 0;
	for(int i = 0; i< v.size(); i++){
		
		l1 = v[i].l;
		s = l1.size() -  1;
		for(it = l1.begin(); it != l1.end(); ++it){
			it -> nfnodes = s;
			s--;
		}

		--it;
		ov = 0;
		
		while(it != l1.begin()){
			--it;
			ov += it-> overlap;
			it -> nOverlap = ov/ it->nfnodes;
		}

		v[i].l = l1;

	}		
}


void read_processor::displaycontigs(){

	std::list<node> l1;
	std::list <node> :: iterator it;
	
	for(int i = 0; i< v.size(); i++){
		
		if(v[i].seq.length() > 150){
			
			std::cout << v[i].seq ;

			l1 = v[i].l;
			for(it = l1.begin(); it != l1.end(); ++it){
				std::cout << "\t" << it->SNo  ;
			}
			std::cout << "\n";
		}
	}		
}

void read_processor::displaynodes(){
	std::list<node> l1;

	for(int i = 0; i< v.size(); i++){
		l1 = v[i].l;

		std::cout << "S\t" << l1.front().SNo << "\t" << "*\n";
			
	}		
}
void read_processor::reducelinks(){
	std::list<node> l1,l2;
	for(int i = 0; i< v.size(); i++){
		l1 = v[i].l;
		l2.push_back(l1.back());
		v[i].l = l2;
		l2.clear();
		l1.clear();
	}		
}


*/