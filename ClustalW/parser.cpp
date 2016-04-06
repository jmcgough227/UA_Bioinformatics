#include <map>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#include "parser.h"

// Parses a single sequence file
Sequence::Sequence(ifstream &file, string name) {
	string comment_line;
	string new_comment_line;
	
	// Read in the first comment line
	getline(file, comment_line);
	
	bool done = false;
	
	int gp_count = 0;
	
	do {
		string line;
		string gene_data;
		
		do {
			if(getline(file, line)) {
		
				// Another comment line
				if(line[0] == '>') {
					new_comment_line = line;
					break;
				}
				else {
					// Add the current line to the gene data
					gene_data += line;
				}
			}
			else {
				done = true;
				break;
			}
		} while(true);
		
		// Find the name of the gene in the comment string
		string gene_name;
		size_t pos = comment_line.find("gene=");
		
		// Make sure we actually found it!
		if(pos == string::npos) {
			cout << "Error! Did not find the gene name!";
		}
		else {
			// Consume the "gene="
			pos += 5;
			
			// Grab the name
			while(pos < comment_line.length() && comment_line[pos] != ']') {
				gene_name += comment_line[pos++];
			}
			
			if(gene_name == "GP") {
				if(gp_count < 2) {
					genes[gene_name + to_string(++gp_count)] = gene_data;
				}
			}
			else {
				genes[gene_name] = gene_data;
			}
			
		}
		
		comment_line = new_comment_line;
	} while(!done);
	
	file_name = name;
}