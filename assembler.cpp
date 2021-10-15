#include "assembler.h"
void assembler::getmap(int curr_klen){
	klen = curr_klen;
	for (int i = 0; i< v.size() ; i++){
		if(v[i].seq.length() > klen)
			setmap(v[i].seq , i);
	}  
}

void assembler::setmap(std::string seq, int i){
	std::string prefix = get_prefix(seq, klen);
	pmap.insert(std::pair <std::string, int> (prefix, i)) ;
}

void assembler::infomap(){
	std::cout << pmap.size() << " out of " << pmap.max_size () << "\n" ; 
}

void assembler::freemap(){
	pmap.clear();
}	

void assembler::displaymap(){
	std::multimap <std::string, int> :: iterator itr ;
	std::cout << "Map:\n";
	for (itr = pmap.begin(); itr != pmap.end(); ++itr)                            
        std::cout  <<  itr->first << '\t' << itr->second << '\n';	
}

void assembler::run_match(){
	std::vector<int> mark(v.size(),0);
	std::vector<contig> vnew;
	// if(klen < 110)
	// 	std::cout << v.size() << "\n";

	for (int i = 0; i< v.size(); i++){	

		find(v[i].seq,v[i].l,i,mark,vnew);                   // look for a match in the map
	}

	for (int i = 0; i<v.size();i++){
		if(mark[i] == 0){
			vnew.push_back(v[i]);
		}		
	}
	
	v.clear();
	v = vnew;
	mark.clear();
	vnew.clear();

}

void assembler::find (std::string seq, std::list <node> & sList,int i, std::vector <int> &mark, std::vector<contig> &vnew){	
	std::multimap <std::string, int> :: iterator itr ;
	std::pair <std::multimap<std::string,int>::iterator, std::multimap<std::string,int>::iterator > ret;
	std::list<node> pList,tempList;
	int count;
	std::string suffix;

	suffix = get_suffix(seq,klen);               // take suffix of the read
	count = pmap.count(suffix);                     

	if(count == 1){
		
		itr = pmap.find(suffix);
		int loc = itr->second ;

		// if(loc != i){
			pList = v[loc].l;
			// if(sList.back() == pList.front())
			// 	pList.pop_front();
			sList.back().overlap = klen;
			sList.insert(sList.end(),pList.begin(),pList.end() );
			vnew.push_back({join(seq, v[loc].seq, klen),sList});
			mark[loc] = 1;
			mark [i] = 1;
			pList.clear();	
		// }

	}

	else if (count > 1){

		ret = pmap.equal_range(suffix);
		mark [i] = 1;
		int loc;
		int jnPt = 0;
		std::string firstSNo,currSNo;			

// check if i is a junction node or not


		for (itr=ret.first; itr!=ret.second; ++itr){
			loc = itr->second;
			if (loc != i){
				currSNo = v[loc].l.begin()->SNo;
				if(itr == ret.first)
					firstSNo = currSNo;
				if(firstSNo != currSNo){
					jnPt = 1;
					junctionPts.push_back(sList.back().SNo);

				}
			}
		}


	    for (itr=ret.first; itr!=ret.second; ++itr){
			tempList = sList;
			loc = itr->second ;
			// if (loc != i){
				pList = v[loc].l;
				// if(tempList.back() == pList.front())
				// 	pList.pop_front();	
				tempList.back().overlap = klen;
				tempList.back().isJnpt = jnPt;
				tempList.insert(tempList.end(),pList.begin(),pList.end() );
				vnew.push_back({join(seq, v[loc].seq, klen),tempList});	
				mark[loc] = 1;					
			// }
			pList.clear();
			tempList.clear();

		}
	}

	else
		mark[i] = 0;
}

// void assembler::run_cycle(){

// 	sort( junctionPts.begin(), junctionPts.end() );
// 	junctionPts.erase( unique( junctionPts.begin(), junctionPts.end() ), junctionPts.end() );

// 	for (int i =0; i< junctionPts.size(); i++){
// 		std::string curr_jnpt = junctionPts[i];
// 		store_vectors(curr_jnpt);
// 		flushLen(i);
// 		flushNov(i);
// 		cjn.clear();
// 	}

// }

// void assembler::store_vectors(std::string curr_jnpt){
// 	std::list<node> l1;	
// 	std::list <node> :: iterator it;
//     maxLen = 0;
// 	maxNOv = 0;

// 	for (int i = 0; i < v.size(); i++){
// 		l1 = v[i].l;
// 		for(it = l1.begin(); it != l1.end(); ++it){
			
// 			if (it-> SNo == curr_jnpt){
				
// 				if(maxLen < it->nfnodes)
// 					maxLen = it->nfnodes;

// 				if(maxNOv < it->nOverlap)
// 					maxNOv = it->nOverlap;		

// 				cjn.push_back({i,it->nfnodes,it->nOverlap});

// 			} 
// 		}

// 	}

// }

// void assembler::flushLen(int iter){
// 	std::vector<int> mark(v.size(),0);
// 	std::vector<contig> vnew;
	
// 	for (int i =0; i< cjn.size(); i++){

// 		if( cjn[i].len < (thLen * maxLen)  ){
// 				mark[i] = 1;
// 		}
// 	}

// 	for (int i = 0; i<v.size();i++){

// 		if(mark[i] == 0){
// 			vnew.push_back(v[i]);
// 		}		
// 	}
	
// 	v.clear();
// 	v = vnew;
// 	mark.clear();
// 	vnew.clear();

// }	
// void assembler::flushNov(int iter){
// 	std::vector<int> mark(v.size(),0);
// 	std::vector<contig> vnew;
	
// 	for (int i =0; i< cjn.size(); i++){

// 		if( cjn[i].nOv < (thNov * maxNOv) ){
// 				mark[i] = 1;
// 		}
// 	}

// 	for (int i = 0; i<v.size();i++){

// 		if(mark[i] == 0){
// 			vnew.push_back(v[i]);
// 		}		
// 	}
	
// 	v.clear();
// 	v = vnew;
// 	mark.clear();
// 	vnew.clear();

// }