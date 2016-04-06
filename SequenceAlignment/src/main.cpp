#include "alignment_utilities.h"
#include "parse.h"
#include "time.h"
#include <chrono>
#include <vector>

using namespace std;
using namespace Parser;

int main(int argc, char *argv[]){
	vector<Sequence> s;
	
	chrono::time_point<chrono::system_clock> start, end;
  start = chrono::system_clock::now();
  
	load_sequences("./sequences/CDS/", s);
	compare_sequences(s);
	
	end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end - start;
  cout << "Output written to: output.txt" << endl;
  cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
	
	return 0;

}
