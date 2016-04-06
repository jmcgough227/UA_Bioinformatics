#ifndef PARSER_H
#define PARSER_H

#include <map>

using namespace std;

struct Sequence {
	string file_name;
	map<string, string> genes;
	
	Sequence(ifstream &file, string name);
};

#endif