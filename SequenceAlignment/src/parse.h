#ifndef PARSE
#define PARSE

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

namespace Parser {

struct Sequence {
	string file_name;
	map<string, string> genes;
	
	Sequence(ifstream &file, string name);
};

void load_sequences(string dir, vector<Sequence> &s);
void compare_sequences(vector<Sequence> &s);






}
#endif
