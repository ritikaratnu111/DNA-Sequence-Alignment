#include "block_assembler.h"


void block_assembler :: no_of_contigs(){
	std::cout << s.size() << "\n";
}

void block_assembler :: push_block(std::vector<contig> temp){

	s.insert(std::end(s), std::begin(temp), std::end(temp));
}

void block_assembler :: push_block_size(long long int size){
	vecsizes.push_back(size);
}

void block_assembler :: set_no_of_rounds(int NO_OF_CLUSTERS){

	no_of_rounds = log2(NO_OF_CLUSTERS);
}

void block_assembler :: setmap(std::string seq, int i){
	std::string prefix = get_prefix(seq, klen);
	pmap.insert(std::pair <std::string, int> (prefix, i)) ;
}


void block_assembler :: getmap( int curr_klen){
	klen = curr_klen;
	for (int i = 0; i< v.size() ; i++){
			setmap(v[i].seq , i);
	}
}

void block_assembler::infomap(){
	std::cout << pmap.size() << " out of " << pmap.max_size () << "\n" ; 
}

void block_assembler::freemap(){
	pmap.clear();
}	

void block_assembler::find(int i,std::vector<int> &markU,std::vector<int> &markV){
	std::multimap <std::string, int> :: iterator itr ;
	std::pair <std::multimap<std::string,int>::iterator, std::multimap<std::string,int>::iterator > ret;	
	std::string seq,suffix;
	std::list <node> sList,pList,tempList;
	
	seq = u[i].seq;
	suffix = get_suffix(seq,klen); 
	sList = u[i].l;
	int count = pmap.count(suffix); 
	int loc;

	if(count == 1){
		itr = pmap.find(suffix);
		loc = itr->second ;	
		if(markV[loc] == 0){
			pList = v[loc].l;	
			sList.back().overlap = klen;
			sList.insert(sList.end(),pList.begin(),pList.end() );
			u[i].seq = join(seq, v[loc].seq, klen); u[i].l = sList;
			markU[i] = 1;
			markV[loc] = 1;			
		}

	}

	else if (count > 1){
		ret = pmap.equal_range(suffix);	
		int max = 0;
		int loc_max;
		int loc_hit = 0;
		for (itr=ret.first; itr!=ret.second; ++itr){  //find the extension which gives the longest contig
			loc = itr ->second;
			if(markV[loc] == 0){
				if(v[loc].seq.length() > max){
					max = v[loc].seq.length();
					loc_max = loc;  
					loc_hit = 1;           
				}				
			}
		}
		if(loc_hit ==1){
			markU[i] = 1;
			markV[loc_max] = 1;
			pList = v[loc_max].l;
			sList.back().overlap = klen;
			sList.insert(sList.end(),pList.begin(),pList.end() );
			u[i].seq = join(seq, v[loc].seq, klen); u[i].l = sList;					
		}


	}	

}

void block_assembler :: run_match(std::vector<int> &markU,std::vector<int> &markV){
	
	for (int i = 0; i< u.size(); i++){	
		if(markU[i] == 0)
			find(i,markU,markV);
		else
			continue;  //has already been extended in previous rounds with a bigger overlap length
	}

}

void block_assembler :: match_kernel(){

	std::vector<int> markU(u.size(),0);
	std::vector<int> markV(v.size(),0);
	std::vector<contig> tempU;


	for(int len = 150; len >=13; len--){
		getmap(len);
		run_match(markU,markV);
		std::cout << "olen: " << len << "\n";
		// std::cout << "olen and vnew size:" << len << "\t" <<  vnew.size() << "\n";
		freemap();
	}


	// for(int i = 0; i< u.size(); i++){

	// 	if (markU[i] == 1 ){
	// 		tempU.push_back(u[i]);
	// 	}

	// 	else{

	// 		if(markU[i] == 0 && u[i].seq.length() > (maxlenU * th)){
	// 		// if(markU[i] == 0 && u[i].seq.length() > (maxlenU * th)){
	// 			tempU.push_back(u[i]);
	// 		}
	// 	}

	// }

	// u.clear();
	// u = tempU;
	// tempU.clear();


	for(int i = 0; i < v.size(); i++){
		if(markV[i] == 0){
		 	//if( v[i].seq.length() > (maxlenV * th) ){
				u.push_back(v[i]);

			//}
		}
	}


}


void block_assembler :: getU(int iter){

	long long int start,size;
	u.clear();


	if(iter == 0){
		start = 0;
		size = vecsizes[0];
	}
	else{
		start = vecsizes[iter - 1];
		size = vecsizes[iter] - vecsizes[iter - 1];			
	}
	std::cout << "u size " << size<< "\n";
	u.insert(std::end(u), std::begin(s) + start, std::begin(s) + start + size);
}

long long int block_assembler :: get_max_length_u(){
	long long int max = 0;

	for(int i = 0; i< u.size();i++){

		if(u[i].seq.length() > max)
			max = u[i].seq.length();
	}

	return max;

}

long long int block_assembler :: get_max_length_v(){
	long long int max = 0;

	for(int i = 0; i< v.size();i++){

		if(v[i].seq.length() > max)
			max = v[i].seq.length();
	}

	return max;

}

void block_assembler :: getV(int iter){

	long long int start,size;
	v.clear();

	start = vecsizes[iter];
	size = vecsizes[iter + 1] - vecsizes[iter];			

	std::cout << "v size " << size<< "\n";
	v.insert(std::end(v), std::begin(s) + start, std::begin(s) + start + size);
}

void block_assembler :: block_assembly(){
	std::vector<contig> temps;
	std::vector<long long int>tempvecsizes;
	long long int size = 0;		
	std::cout << "No of blocks " << vecsizes.size() << "\n";
	for(int iter = 0 ; iter < vecsizes.size(); iter += 2){
		
		getU(iter);
		getV(iter);		
		maxlenU = get_max_length_u();
		maxlenV = get_max_length_v();
		std::cout << "maxlenU " << maxlenU << "\t"<< "maxlenV " << maxlenV << "\n";

		match_kernel();  // try for extension


		size += u.size(); 
		std::cout << "Size " << size << "\n";
		tempvecsizes.push_back(size);
		temps.insert(std::end(temps), std::begin(u), std::end(u)); //push assembled blocks in a temp vector
		u.clear();
		v.clear();

	}

	s = temps ;                          //reassignment
	vecsizes = tempvecsizes;
	temps.clear();
	tempvecsizes.clear();		

}



