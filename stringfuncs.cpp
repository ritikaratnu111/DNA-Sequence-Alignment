#include <iostream>
#include<string>

char base_complement(char base)
{
	char comp;
	switch (base)
   {
       case 'A': comp = 'T';
               break;
       case 'C': comp = 'G';
                break;
       case 'G': comp = 'C';
               break;
       default: comp = 'A';
                break;  
   }
   return comp;

} 
std::string rev_complement(std::string s){
	std::string rcomp = s;
	for (int i=0; i<s.length(); i++)
		rcomp[i] = base_complement(s[i]);
	return rcomp;
}

std::string join (std::string s1, std::string s2, int len){
	return(s1 + s2.substr( len , s2.length() ) );
}


std::string get_prefix(std::string read ,int len){	
	if (len <= read.length()){
		return read.substr(0,len) ;
	}

	else 
		std::cout << "\n" << "Input length longer than string" ;
}

std::string get_suffix(std::string read, int len){
	if (len <= read.length()){
		return read.substr(read.length()-len,read.length()) ;
	}

	else 
		std::cout << "\n" << "Input length longer than string" ;
}	
