#ifndef STRUCTDEFS_H
#define STRUCTDEFS_H

struct node{ 
   std::string SNo;
   int overlap;
   int isJnpt; //not used in assembly, only for pruning
   int nfnodes; // no of nodes that follow this node in the contig //not used in assembly
   float nOverlap; // normalised overlap of the following nodes //not used in assembly
}; 

struct contig{ 
   std::string seq;
   std::list <node> l;
}; 



#endif // STRUCTDEFS_H