#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <unistd.h>
#include <ctime>
#include <sstream>

using namespace std;

// Overload the stream operator for easily printing out a vector of char
ostream& operator<<(std::ostream& os, const vector<char>& data) {
  for (auto& c : data)
  	os << c;
}

// Dumps a file into a vector of chars
bool load_file(string name, vector<char> &data) {
	ifstream file;
	
	file.open(name.c_str());
	
	if(!file.is_open())
		return false;
	
	char c;
	while(file >> c) {
		data.push_back(c);
	}

	file.close();
	
	return true;
}

// Holds a count of how many normal strands and how many prime strands
// are of a given length (put into a map)
struct Strand {
	int count_normal;
	int count_prime;
};

// Determines whether a strand will copy successfully
// This takes into account the length of the strand to copy, and whether it is
// a normal or a prime strand
bool copy_strand(int length, bool normal_strand, int start_pos, int end_pos, int &pos) {
	// To determine whether the strand is copied successfully, pick a random
	// number between [0, length / fail_rate]. If the number is less than
	// the randomly picked number, assume that the copy failed. Use this
	// random number as the position that the failure happened at.
	
	int real_length;
	
	if(normal_strand)
		real_length = length - start_pos + 1;
	else
		real_length = end_pos + 1;
	
	// The fail rate of the PCR algorithm
	double fail_rate = .05;
	
	// Maximum number to pick
	int max_num = ((double) end_pos - start_pos + 1) / fail_rate;
	int num = rand() % max_num;
	
	// Find the position that it stops copying at
	if(num < real_length) {
		if(normal_strand)
			pos = start_pos + num;
		else
			pos = end_pos - num;
			
		return false;
	}
	else {
			if(normal_strand)
				pos = length - 1;
			else
				pos = 0;
				
		return true;
	}
}

// Performs the actual PCR algorithm
void run_pcr(int total_length, int start_pos, int end_pos, map<int, Strand> &strand_old) {
	// The minimum length for a strand to be considered "good"
	int good_strand_length = (end_pos - start_pos + 1);
	
	// Make a copy of the map
	map<int, Strand> strands = strand_old;
	
	// Iterate over all the strands
	for(auto &str : strand_old) {
		Strand &s = str.second;
		int len = str.first;
		
		if(len >= good_strand_length) {;
			// Copy the normal strands, which will become prime strands
			for(int i = 0;i < s.count_normal;i++) {
				int pos;
				
				copy_strand(total_length, true, start_pos, end_pos, pos);
				
				int slength = pos - start_pos + 1;
				
				if(slength > good_strand_length && len != total_length)
					slength = good_strand_length;
				
				strands[slength].count_prime++;
			}
			
			// Copy the prime strands, which will become the normal strands
			for(int i = 0;i < s.count_prime;i++) {
				int pos;
			
				copy_strand(total_length, false, start_pos, end_pos, pos);
				
				int slength = end_pos - pos + 1;
				
				if(slength > good_strand_length && len != total_length)
					slength = good_strand_length;
				
				strands[slength].count_normal++;
			}
		}
	}
	
	// Move the copied map back to the original
	strand_old = strands;
}


// Calculate the primers based on the input file
bool get_primers(string& left_primer, string& right_primer) {
	ifstream file;
	file.open("primer3_output.txt");

	if (!file.is_open())
		return false;
	
	string str;
	while(getline(file, str)) {
		if (str.substr(0,str.find('=')) == "PRIMER_LEFT_0_SEQUENCE")
			left_primer = str.substr(str.find('=') + 1);
		else if (str.substr(0,str.find('=')) == "PRIMER_RIGHT_0_SEQUENCE")
			right_primer = str.substr(str.find('=') + 1);
	}

	file.close();
	return true;
}

// Create a config file to calculate the primers
bool create_config_file(const int& start_pos, const int& end_pos, const vector<char>& data) {
	ofstream file;
	file.open("primer_config.txt", ofstream::out | ofstream::trunc);

	if (!file.is_open())
		return false;

	file << "SEQUENCE_ID=ebola\n";
	file << "SEQUENCE_TEMPLATE=" << data;
	file << "\nSEQUENCE_INCLUDED_REGION=" << start_pos << "," << end_pos << "\n";
	file << "PRIMER_NUM_RETURN=1\n";
	file << "PRIMER_TASK=generic\n";
	file << "PRIMER_PICK_LEFT_PRIMER=1\n";
	file << "PRIMER_PICK_RIGHT_PRIMER=1\n";
	file << "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" << getcwd(NULL,0) 
			<< "/primer3-2.3.6/src/primer3_config/\n=";

	file.close();

	return true;
}

int main(int argc, char *argv[]) {
	map<int, Strand> strands;
	vector<char> sequence;

	// Get the command line arguments
	string file_name = argv[1];
	if (!load_file(file_name, sequence)) {
		cout << "Error: sequence file could not be loaded\n";
		return 1;
	}
	
	int start_pos = atoi(argv[2]);
	int end_pos = atoi(argv[3]);
	int count = atoi(argv[4]);
	
	if (!create_config_file(start_pos, end_pos, sequence)) {
		cout << "Error: primer config file could not be created\n";
		return 1;
	}

	// Calculate the primers
	string cwd = getcwd(NULL,0);
	string primer3_command = "primer3-2.3.6/src/primer3_core -output=" 
		  + cwd + "/primer3_output.txt primer_config.txt &";

	system(primer3_command.c_str());

	string left_primer = "";
	string right_primer = "";

	if (!get_primers(left_primer, right_primer)) {
		cout << "Error: primer3 output file could not be loaded\n";
		return 1;
	}

	// Seed the random number generator
	srand(time(NULL));
	
	// Initialize the map with counts of 0
	for(int i = 0;i < sequence.size()+1;i++) {
		Strand s;
		s.count_normal = 0;
		s.count_prime = 0;
		
		strands[i] = s;
	}
	
	// Add in the first strand
	strands[sequence.size()].count_prime = 1;
	strands[sequence.size()].count_normal = 1;
	
	// Run PCR
	for(int i = 0;i < count;i++) {
		run_pcr(sequence.size(), start_pos, end_pos, strands);
	}
	
	double total = 0;
	double n = 0;
	stringstream s;

	ofstream out_file("output");

	s << endl << endl << "Length\t\tNormal Count\tPrime Count" << endl;
	s << "-------------------------------------------" << endl;

	// Calculate the percentage of strands which were copied successfully
	// Also, write the output file
	for(auto &i : strands) {
		if(i.first >= (end_pos - start_pos + 1)) {
			s << i.second.count_normal << i.second.count_prime;
			total += i.second.count_normal + i.second.count_prime;
		}
		
		s << i.first << "\t\t" << i.second.count_normal << "\t\t" << i.second.count_prime << endl;
		n += i.second.count_normal + i.second.count_prime;
	}
	
	cout << "Percent passed: " << 100.0 * (total / n) << endl;
	cout << "Please see file 'output' for the distrubution of strands." << endl;
	
	
	out_file << "Percent passed: " << 100.0 * (total / n) << endl;
	
	out_file << "left primer: " << left_primer << endl;
	out_file << "right primer: " << right_primer << endl;
	
	
	out_file << s.str();
}

