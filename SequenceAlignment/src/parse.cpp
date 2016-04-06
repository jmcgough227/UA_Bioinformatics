#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "parse.h"
#include "alignment_utilities.h"

using namespace std;

namespace Parser {

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

// Loads all of the sequences
// The first sequence is the Sudan 1976, and the second is the Zaire 1976
void load_sequences(string dir, vector<Sequence> &s) {
	// Load in the Sudan and Zaire sequences
	ifstream sudan((dir + "Sudan1976.txt").c_str());
	ifstream zaire((dir + "Zaire1976.txt").c_str());
	
	if(!sudan.is_open() || !zaire.is_open()) {
		cout << "Error: failed to load base files (Zaire and Sudan)" << endl;
		throw -1;
	}
	
	s.push_back(Sequence(sudan, "Sudan 1976"));
	s.push_back(Sequence(zaire, "Zaire 1976"));
	
	// Load in the other sequences
	for(int i = 40;i <= 118;i++) {
		char buf[10];
		sprintf(buf, "%.3d", i);
		
		string name = "KM233" + (string)buf + ".txt";
		string file_name = dir + name;;
		ifstream file(file_name.c_str());
		
		if(!file.is_open()) {
			cout << "Error: failed to open file: " << file_name << endl;
			throw -1;
		}
		
		s.push_back(Sequence(file, name));
		
	}
}

double average(vector<double> &data) {
	double sum = 0;
	
	for(int i = 0;i < data.size();i++) {
		sum += data[i];
	}
	
	return sum / data.size();
}

double std_deviation(vector<double> &data) {
	double a = average(data);
	double sum = 0;
	
	for(int i = 0;i < data.size();i++) {
		sum += pow(data[i] - a, 2);
	}
	
	return sqrt(sum / (data.size() - 1));
}

struct StatInfo {
	double average;
	double deviation;
	
	StatInfo(double a, double s) : average(a), deviation(s) {}
};

// Truncates a string to at most 6 digits
double truncate6(double n) {
	string s = to_string(n).substr(0, 6);
	
	return atof(s.c_str());
}

struct StatReport {
	string name;
	vector<StatInfo> data;
	
	void print(ofstream &file) {
		
		file << name << " Average" << "\t";
		file << truncate6(data[0].average) << "\t\t";
		file << truncate6(data[1].average)  << "\t\t";
		file << truncate6(data[2].average) << "\t\t\t";
		file << truncate6(data[3].average) << "\t\t\t\t";
		file << truncate6(data[4].average) << "\t\t";
		file << truncate6(data[5].average) << "\t\t";
		file << endl;
		
		file << name << " Std Dev" << "\t";
		file << truncate6(data[0].deviation) << "\t\t";
		file << truncate6(data[1].deviation) << "\t\t";
		file << truncate6(data[2].deviation) << "\t\t\t";
		file << truncate6(data[3].deviation) << "\t\t\t\t";
		file << truncate6(data[4].deviation) << "\t\t";
		file << truncate6(data[5].deviation) << "\t\t";
		file << endl << endl;
		
	}
	
	StatReport(string n, vector<StatInfo> d) : name(n), data(d) {}
};

// Calulates average and standard deviation of alignmentStats
vector<StatInfo> get_sequence_stats(vector<alignmentStats> &stats) {
	vector<StatInfo> data;
	vector<double> v;
	
	for(int i = 0;i < stats.size();i++)
		v.push_back(stats[i].insertion);
	
	data.push_back(StatInfo(average(v), std_deviation(v)));
	v.clear();
	
	for(int i = 0;i < stats.size();i++)
		v.push_back(stats[i].deletion);
	
	data.push_back(StatInfo(average(v), std_deviation(v)));
	v.clear();
	
	for(int i = 0;i < stats.size();i++)
		v.push_back(stats[i].synMutation);
	
	data.push_back(StatInfo(average(v), std_deviation(v)));
	v.clear();
	
	for(int i = 0;i < stats.size();i++)
		v.push_back(stats[i].nonSynMutation);
	
	data.push_back(StatInfo(average(v), std_deviation(v)));
	v.clear();
	
	for(int i = 0;i < stats.size();i++)
		v.push_back(stats[i].match);
	
	data.push_back(StatInfo(average(v), std_deviation(v)));
	v.clear();
	
	for(int i = 0;i < stats.size();i++)
		v.push_back(stats[i].mismatch);
	
	data.push_back(StatInfo(average(v), std_deviation(v)));
	
	return std::move(data);
} 

// Compares all of the sequences against Sudan 1976 and Zaire 1976
void compare_sequences(vector<Sequence> &s) {
	map<string, alignmentStats> sudan_total;
	map<string, alignmentStats> zaire_total;
	
	vector<alignmentStats> sudan_data;
	vector<alignmentStats> zaire_data;
	
	ofstream file("output.txt");
	
	cout << "Aligning to Sudan 1976: 0%";
	cout.flush();
	
	int p = 0;
	
	// Compare all the strands against Sudan 1976
	for(int i = 2;i < s.size();i++) {
		alignmentStats stat;
		
		for(auto gene : s[i].genes) {
			// If the gene exists in the sequence but not Sudan 1976, ignore it
			if(s[0].genes.count(gene.first) > 0) {
				string &sudan_gene_data = s[0].genes[gene.first];
				string &other_gene_data = gene.second;
				
				vector<string> alignment = hirschbergAlignment(sudan_gene_data, other_gene_data);
				
				alignmentStats new_stat = gatherStats(alignment);
				sudan_total[gene.first] += new_stat;
				stat += new_stat;
				//cout << "i = " << i << endl;
			}
		}
		
		sudan_data.push_back(stat);
		
		string p_string = to_string(p);
		
		for(int i = 0;i < p_string.length() + 1;i++)
			cout << "\b";
		
		p = i*100 / s.size();
		
		cout << p << "%";
		
		cout.flush();
		
	}
	
	cout << "\b\b\b" << "100%" << endl;
	
	cout << "Aligning to Zaire 1976: 0%";
	
	// Zaire
	p = 0;
	
	// Compare all the strands against Zaire 1976
	for(int i = 2;i < s.size();i++) {
		alignmentStats stat;
		
		for(auto gene : s[i].genes) {
			// If the gene exists in the sequence but not Zaire 1976, ignore it
			if(s[1].genes.count(gene.first) > 0) {
				string &zaire_gene_data = s[1].genes[gene.first];
				string &other_gene_data = gene.second;
				
				vector<string> alignment = hirschbergAlignment(zaire_gene_data, other_gene_data);
				
				alignmentStats new_stat = gatherStats(alignment);
				zaire_total[gene.first] += new_stat;
				stat += new_stat;
				
				//cout << "i = " << i << endl;
			}
		}
		
		zaire_data.push_back(stat);
		
		string p_string = to_string(p);
		
		for(int i = 0;i < p_string.length() + 1;i++)
			cout << "\b";
		
		p = i*100 / s.size();
		
		cout << p << "%";
		
		cout.flush();
		
	}
	
	cout << "\b\b\b" << "100%" << endl;

	// Print Sudan data
	file << "/////////////////////////////////////////////////////////////////////////////////////////"
		"////////////" << endl;
		
	file << "//                            Table 1 (Sudan 1976 against 2014 Strands)                            //"
		<<    "\n//////////////////////////////////////////////////////////////////////////////////////////"
		"///////////\n\n";

	file << "Sudan Ebola Strain\tInsertions\tDeletions\tSynonymous Mutations\tNon-synonymous Mutations\tMatch\t\tMismatch" << endl;
	file << "--------------------------------------------------------------------------------------"
		"--------------------------------------------------------" << endl;
	
	for(int i = 2;i < sudan_data.size() + 2;i++) {
		file << s[i].file_name << "\t\t";
		file << sudan_data[i - 2].insertion << "\t\t";
		file << sudan_data[i - 2].deletion << "\t\t";
		file << sudan_data[i - 2].synMutation << "\t\t\t";
		file << sudan_data[i - 2].nonSynMutation << "\t\t\t\t";
		file << sudan_data[i - 2].match << "\t\t";
		file << sudan_data[i - 2].mismatch << "\t\t";
		
		file << endl;
	}
	
	
	// Print Zaire data
	file << "\n\n\n";
	file << "/////////////////////////////////////////////////////////////////////////////////////////"
		"////////////" << endl;
		
	file << "//                            Table 2 (Zaire 1976 against 2014 Strands)                            //"
		<<    "\n//////////////////////////////////////////////////////////////////////////////////////////"
		"///////////\n\n";

	file << "Zaire Ebola Strain\tInsertions\tDeletions\tSynonymous Mutations\tNon-synonymous Mutations\tMatch\t\tMismatch" << endl;
	file << "--------------------------------------------------------------------------------------"
		"--------------------------------------------------------" << endl;
	
	for(int i = 2;i < zaire_data.size() + 2;i++) {
		file << s[i].file_name << "\t\t";
		file << zaire_data[i - 2].insertion << "\t\t";;
		file << zaire_data[i - 2].deletion << "\t\t";;
		file << zaire_data[i - 2].synMutation << "\t\t\t";;
		file << zaire_data[i - 2].nonSynMutation << "\t\t\t\t";;
		file << zaire_data[i - 2].match << "\t\t";;
		file << zaire_data[i - 2].mismatch << "\t\t";;
		
		file << endl;
	}
	
	// Print Sudan gene stats
	file << "\n\n\n";
	file << "/////////////////////////////////////////////////////////////////////////////////////////"
		"////////////" << endl;
		
	file << "//                               Table 3 (Sudan Individual Genes)                                  //"
		<<    "\n//////////////////////////////////////////////////////////////////////////////////////////"
		"///////////\n\n";

	file << "Gene\tInsertions\tDeletions\tSynonymous Mutations\tNon-synonymous Mutations\tMatch\t\tMismatch" << endl;
	file << "--------------------------------------------------------------------------------------"
		"--------------------------------------------------------" << endl;
	
	vector<alignmentStats> sudan_total_vector;
	vector<alignmentStats> zaire_total_vector;
	
	for(auto i : sudan_total) {
		file << i.first << "\t";
		file << sudan_total[i.first].insertion << "\t\t";
		file << sudan_total[i.first].deletion << "\t\t";
		file << sudan_total[i.first].synMutation << "\t\t\t";
		file << sudan_total[i.first].nonSynMutation << "\t\t\t\t";
		file << sudan_total[i.first].match << "\t\t";
		file << sudan_total[i.first].mismatch << "\t\t";
		
		sudan_total_vector.push_back(sudan_total[i.first]);
		
		file << endl;
	}
		
	// Print Zaire gene stats
	file << "\n\n\n";
	file << "/////////////////////////////////////////////////////////////////////////////////////////"
		"////////////" << endl;
		
	file << "//                               Table 4 (Zaire Individual Genes)                                  //"
		<<    "\n//////////////////////////////////////////////////////////////////////////////////////////"
		"///////////\n\n";

	file << "Gene\tInsertions\tDeletions\tSynonymous Mutations\tNon-synonymous Mutations\tMatch\t\tMismatch" << endl;
	file << "--------------------------------------------------------------------------------------"
		"--------------------------------------------------------" << endl;
	
	for(auto i : zaire_total) {
		file << i.first << "\t";
		file << zaire_total[i.first].insertion << "\t\t";
		file << zaire_total[i.first].deletion << "\t\t";
		file << zaire_total[i.first].synMutation << "\t\t\t";
		file << zaire_total[i.first].nonSynMutation << "\t\t\t\t";
		file << zaire_total[i.first].match << "\t\t";
		file << zaire_total[i.first].mismatch << "\t\t";
		
		zaire_total_vector.push_back(zaire_total[i.first]);
		
		file << endl;
	}
	
	file << "\n\n\n";
	file << "/////////////////////////////////////////////////////////////////////////////////////////"
		"////////////" << endl;
		
	file << "//                               Table 5 (Statistics for Tab 1-4)                                  //"
		<<    "\n//////////////////////////////////////////////////////////////////////////////////////////"
		"///////////\n\n";

	file << "Stat\t\t\tInsertions\tDeletions\tSynonymous Mutations\tNon-synonymous Mutations\tMatch\t\tMismatch" << endl;
	file << "--------------------------------------------------------------------------------------"
		"--------------------------------------------------------" << endl;
	
	// Calculate averages and standard deviation for the tables
	StatReport report("Sudan Seq", get_sequence_stats(sudan_data));
	report.print(file);
	
	report = StatReport("Zaire Seq", get_sequence_stats(zaire_data));
	report.print(file);
		
	report = StatReport("Sudan Gene", get_sequence_stats(sudan_total_vector));
	report.print(file);
	
	report = StatReport("Zaire Gene", get_sequence_stats(zaire_total_vector));
	report.print(file);
}
	
	
	
	
	
	
	
	
	
	
	
	
	








}
