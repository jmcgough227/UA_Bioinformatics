#ifndef ALIGNMENT_UTILITIES
#define ALIGNMENT_UTILITIES

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <locale>
#include <utility>
#include <map>

struct alignmentStats {
  int insertion;
  int deletion;
  
  int synMutation;
  int nonSynMutation;
  
  int match;
  int mismatch;
  
  alignmentStats() {
		insertion = 0;
		deletion = 0;
		synMutation = 0;
		nonSynMutation = 0;
		match = 0;
		mismatch = 0;
  }
};

alignmentStats gatherStats( const std::vector<std::string>& alignment );

void operator+=(alignmentStats &a, const alignmentStats &b);

std::string reverse(const std::string&);
int partition(const std::vector<int>&, const std::vector<int>&);

bool getSequence(const std::string&, std::string&);

std::vector<int> NWScore(const std::string&, const std::string&);
std::vector<std::string> NWAlignment(const std::string&, const std::string&);
std::vector<std::string> hirschbergAlignment(const std::string&, const std::string&);
const int ScoreMatrix[6][6] = {
      /* X   A   C   G   T   U
/*X*/{  -1, -1, -1, -1, -1, -1},
/*A*/{  -1,  1, -1, -1, -1, -1},
/*C*/{  -1, -1,  1, -1, -1, -1},
/*G*/{  -1, -1, -1,  1, -1, -1},
/*T*/{  -1, -1, -1, -1,  1, -1},
/*U*/{  -1, -1, -1, -1, -1,  1}};


const char BaseIndices[128] = 
{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //ascii 0-31
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,//ascii 32-64
1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,4,5,0,0,0,0,0,          //A-Z ascii 65-90
0,0,0,0,0,0, 													  //ascii 91-96
1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,4,5,0,0,0,0,0,         //a-z ascii 97-122
0,0,0,0,0}; 												    //ascii 123-127


bool aminoCheck ( const std::string& str1, const std::string& str2 );

inline int Score(char c1, char c2) { 
  return ScoreMatrix[BaseIndices[c1]][BaseIndices[c2]];
}

#endif
