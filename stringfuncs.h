#ifndef FUNCS_H
#define FUNCS_H

char base_complement(char base);

std::string rev_complement(std::string s);

std::string join (std::string s1, std::string s2, int len);

std::string get_prefix(std::string read ,int len);

std::string get_suffix(std::string read, int len);

#endif // FUNCS_H