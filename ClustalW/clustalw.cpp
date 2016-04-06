#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include <algorithm>
#include <quadmath.h>
#include <fstream>
#include <set>

#include "alignment_utilities.h"
#include "parser.h"

using namespace std;

// Set this to change the desired precision for the Float class
using FLOAT = __float128;

// Custom float class for extra precision
struct Float {
	FLOAT data;
	
	Float(FLOAT f) : data(f) {}
		
	Float operator+(const Float &f) {
		Float a = Float(f.data + data);
		
		return a;
	}
	
	Float operator/(const Float &f) {
		Float a = Float(data / f.data);
		
		return a;
	}
	
	Float operator-(const Float &f) {
		Float a = Float(data - f.data);
		
		return a;
	}
	
	Float operator+=(const Float &f) {
		data += f.data;
		
		return *this;
	}
	
	Float operator++(int a) {
		data++;
		
		return *this;
	}
	
	double getDouble() {
		return (double)data;
	}
	
	bool operator<(Float &f) {
		return data < f.data;
	}
	
	Float() {
		data = 0;
	}
};

// Calculates the p-distance between two sequences
Float distance(string &s1, string &s2) {
	Float count = 0;
	
	for(int i = 0;i < s1.length();i++) {
		if(s1[i] != s2[i]) {
			count++;
		}
	}
	
	return count / s1.length();
}


// Sentinal value for indicating the end of the table items list
const int TABLE_END = 1000000;

// Represents a 2-dimensional table where elements (aka nodes) can be deleted
// The table is designed to be symmetrical i.e. table[a][b] == table[b][a], so
// adding table[a][b] also adds table[b][a] and deleting table[a][b] also deletes
// table[b][a]
// Since nodes can be deleted, the table may look something like this:
//
//		1	3	7	9
//		--------------
//	1	x	x	x	x
//	3	x	x	x	x
//	7	x	x	x	x
//	9	x	x	x	x
//
// i.e. the elements do not have to be contiguous
struct Table {
	map<pair<int, int>, Float> data;
	vector<int> items;
	
	set<int> items_set;
	
	
	// Sets table[a][b] and table[b][a] to dist
	void setEntry(int a, int b, Float dist) {
		data[pair<int, int>(a, b)] = dist;
		data[pair<int, int>(b, a)] = dist;
		
		items_set.insert(a);
		items_set.insert(b);
		
#if 0
		bool found = false;
		
		for(int &i : items) {
			if(i == a) {
				found = true;
				break;
			}
		}
		
		if(!found) {
			items.push_back(a);
		}
		
		
		found = false;
		for(int &i : items) {
			if(i == b) {
				found = true;
				break;
			}
		}
		
		if(!found) {
			items.push_back(b);
		}
		
		std::sort(items.begin(), items.end());
#endif
	}
	
	// Gets the value at table[a][b]
	Float get(int a, int b) {
		return data[pair<int, int>(a, b)];
	}
	
	// Removes the node a
	void remove(int a) {
		/*for(int &i : items) {
			if(i == a) {
				i = TABLE_END;
				break;
			}
		}*/
	
#if 0
		vector<int> new_items;
		
		for(int &i : items) {
			if(i != a) {
				new_items.push_back(i);
			}
		}
#endif

		items_set.erase(a);
		
		
		//std::sort(items.begin(), items.end());
	}
	
	// Calulculates the r values for the table
	
#if 0
	vector<Float> calcR() {
		vector<Float> r(items.size());
		
		/*
		for(int i = 0;i < items.size() && items[i] < TABLE_END;i++) {
			r[i] = 0;
			
			for(int d = 0;d < items.size() && items[d] < TABLE_END;d++) {
				r[i] += Float(get(items[i], items[d]) / (count() - 2));
			}
		}
		*/
		
		for(int i : items) {
			for(int d : items) {
				r[i] += Float(get(i, d) / (count() - 2));
			}
		}
		
		return r;
	}
#endif

	map<int, Float> calcR() {
		map<int, Float> r;
		
		for(int i : items_set) {
			Float &rvalue = r[i];
			rvalue = 0;
			
			for(int d : items_set) {
				rvalue += Float(get(i, d) / (count() - 2));
			}
		}
		
		return r;
	}
		
	
	// Calculates the adjusted distance matrix for the table
	Table getDist() {
		Table t = *this;
		map<int, Float> r = calcR();
		
#if 0
		for(int i = 0;i < items.size() && items[i] < TABLE_END;i++) {
			for(int d = 0;d < items.size() && items[d] < TABLE_END;d++) {
				t.setEntry(items[i], items[d], get(items[i], items[d]) - (r[i] + r[d]));
			}
		}
#endif
		for(int i : items_set) {
			for(int d : items_set) {
				t.setEntry(i, d, get(i, d) - (r[i] + r[d]));
			}
			
		}
		
		return t;
	}
	
	// Prints out the table and possibly the r values
	void print(bool p_dist = false) {
		return;
		printf("\t\t");
		
		for(int i : items_set) {
			printf("%7d\t\t\t", i + 1);
		}
		
		printf("\n");
	
		for(int i : items_set) {
			printf("%7d\t\t", i + 1);
			
			for(int d : items_set) {
				if(i != d || p_dist)
					printf("%7f\t\t", get(i, d).getDouble());
				else
					printf("---------\t\t");
			}
			
			cout << endl;
		}
		
#if 0
		if(p_dist && count() > 2) {
			cout << "===R===" << endl;
			vector<Float> r = calcR();
			
			for(int i = 0;i < r.size();i++) {
				if(items[i] < TABLE_END) {
					cout << r[i].getDouble() << endl;
				}
			}
		}
#endif
	}
	
	// Adds a new node and calculates the distance from it to all the other nodes
	void add(int node, int i, int j) {
#if 0
		vector<int> items_copy = items;
		
		for(int m : items_copy) {
			if(m < TABLE_END) {
				setEntry(node, m, (get(i, m) + get(j, m) - get(i, j))/2);
			}
		}
#endif
		set<int> items_copy = items_set;
		
		for(int m : items_copy) {
			setEntry(node, m, (get(i, m) + get(j, m) - get(i, j))/2);
		}
	

	}
	
	// Given a node ID, this returns the position in the items list of the node
#if 0
	int getItemIndex(int a) {
		for(int i = 0;i < items.size();i++) {
			if(items[i] == a) {
				return i;
			}
		}
		
		return -1;
	}
#endif
	
	// Returns the location of the minimum value in the table
	pair<int, int> pickMin() {
#if 0
		int min_i = -1;
		int min_j = -1;
		Float min_d = 1000000;
		
		
		for(auto i : items) {
			for(auto d : items) {
				if(i != d && i < TABLE_END && d < TABLE_END && get(i, d) < min_d) {
					min_d = get(i, d);
					min_i = i;
					min_j = d;
				}
			}
		}
		
		return pair<int, int>(min_i, min_j);
#endif
		Float min_d = 10000000;
		int min_i = -1;
		int min_j = -1;
		
		for(int i : items_set) {
			for(int j : items_set) {
				if(i != j) {
					Float val = get(i, j);
					
					if(val < min_d) {
						min_d = val;
						min_i = i;
						min_j = j;
					}
				}
			}
		}
		
		return pair<int, int>(min_i, min_j);
	}
	
	// Returns the number of items in the table
	int count() {
		
#if 0
		int c = 0;
		
		for(int i : items) {
			if(i < TABLE_END)
				c++;
			else
				break;
		}
		
		return c;
#endif
		cout << "Size: " << items_set.size() << endl;
		return items_set.size();
	}
	
	pair<int, int> getLastNodes() {
		vector<int> nodes;
		
		for(int i : items_set) {
			nodes.push_back(i);
		}
		
		return pair<int, int>(nodes[0], nodes[1]);
	}
};

// Calculates the final distance matrix between each pair of nodes
double final_dist(int node, int start, int end, map<pair<int, int>, Float> edges) {
	double minv = 1000000;
	
	if(node == end) {
		return 0;
	}
	
	for(auto i : edges) {
		int me = node;
		int other;
		
		if(i.first.first == node) {
			other = i.first.second;
		}
		else {
			other = i.first.first;
			me = i.first.second;
		}
		
		if(me == node) {
			map<pair<int, int>, Float> new_edges = edges;
			new_edges.erase(i.first);
			
			double d = final_dist(other, start, end, new_edges) + i.second.getDouble();
			
			if(d < minv) {
				minv = d;
			}
		}
	}
	
	return minv;
}

void clustalw(vector<string> &seq, map<pair<int, int>, Float> &edges, bool print) {
	size_t n = seq.size();		// Total number of sequences
	string align[n][n][2];		// Aligned sequence strings (note: this uses gcc variable length array)
	
	cout << "ClustalW on " << seq.size() << " sequences" << endl;
	cout << "Aligning..." << endl;
	
	// Get the scores for aligning every pair of sequences
	for(int i = 0; i < n; i++) {
		for(int d = 0; d <= i; d++) {
			vector<string> output = hirschbergAlignment(seq[i], seq[d]);
		
			align[i][d][0] = output[0];
			align[i][d][1] = output[2];
			align[d][i][0] = output[0];
			align[d][i][1] = output[2];
		}
	}
	
	Float p_dist[n][n];
	Table p_table;
	
	cout << "Calculating distances..." << endl;
	
	// Calculate the p-distance matrix and add the distance to the p-distance table
	for(int i = 0; i < n; i++) {
		for(int d = 0; d <= i; d++) {
			p_dist[i][d] = distance(align[i][d][0], align[i][d][1]);
			p_dist[d][i] = p_dist[i][d];
			
			p_table.setEntry(i, d, p_dist[i][d]);
		}
	}
	
	cout << "INITIAL DISTANCES" << endl;
	p_table.print(true);
	
	int node_count = n;
	
	Table d_table;
	
	int iter = 0;
	
	cout << "Running Neighbor-Joining" << endl;
	
	int x;
	
	// Keep running Neighbor-Joining until there are only two nodes left in the table
	while(p_table.count() != 2) {
		cout << "Iteration: " << ++iter << endl;
		
		if(print) {
			cout << "ITERATION: " << ++iter << endl;
			cout << "==============P-dist==============" << endl;
			p_table.print(true);
		}
		
		// Get the new modified distance matrix
		d_table = p_table.getDist();
		
		if(print) {
			cout << endl << "==============D-dist==============" << endl;
			d_table.print();
		}
		
		Table L = d_table;

		// Select the position in L which has the minimum value. mini.first and mini.second are the
		// x and y position in the table, but they're also the two nodes which will be removed from the
		// the table and joined to the new node we're going to create
		pair<int, int> mini = L.pickMin();
		
		if(print) {
			cout << endl << "Min i: " << mini.first + 1 << endl;
			cout << "Min j: " << mini.second + 1 << endl;
		}
		
		// Calculate the r value for each row
		map<int, Float> r = p_table.calcR();
		
		// Calculate the edge length of the new node we're going to create
		Float edge_length = (p_table.get(mini.first, mini.second) + r[mini.first] - r[mini.second]) / 2;
		
		// Create the new node and remove the other two
		p_table.add(node_count, mini.first, mini.second);
		//cout << "Size before: " << p_table.count() << endl;
		p_table.remove(mini.first);
		p_table.remove(mini.second);
		//cout << "Size after: " << p_table.count() << endl;
		
		// Add the edges between the new node and the other two
		edges[pair<int, int>(mini.first, node_count)] = edge_length;
		edges[pair<int, int>(mini.second, node_count)] = p_table.get(mini.first, mini.second) - edge_length;
		
		node_count++;
		
		if(print)
			cout << endl << endl << endl << endl;
		
		//cin >> x;
		
	}
	
	if(print) {
#if 0
		cout << "ITERATION: " << ++iter << endl;
		cout << "==============P-dist==============" << endl << endl;
		p_table.print(true);
	
		cout << endl << "Remaining node A: " << p_table.items[0] + 1 << endl;
		cout << "Remaining node B: " << p_table.items[1] + 1 << endl << endl << endl;
#endif
	}
	
	// After running Neighbor-Joining, there will be two nodes left in the table
	// Between these nodes, create an edge between them
	pair<int, int> last_nodes = p_table.getLastNodes();
	Float edge_length = p_table.get(last_nodes.first, last_nodes.second);
	
	edges[pair<int, int>(last_nodes.first, last_nodes.second)] = edge_length;
	
	char output[1024];
	
	if(print) {
		cout << "==============Edges==============" << endl;
		for(auto i : edges) {
			FLOAT a = 1;
			FLOAT b = 3;
			
			quadmath_snprintf(output, 1000, "%.20Qf", i.second.data);
			
			cout << i.first.first + 1 << " <--> " << i.first.second + 1  << " = " << output << endl;
		}
		
		cout << "\n==============Final Distance Matrix==============\n";
		
		printf("\t\t");
		for(int i = 0;i < node_count;i++) {
			printf("%7d\t\t", i + 1);
		}
		
		printf("\n");
		
		for(int i = 0;i < node_count;i++) {
			printf("%7d\t\t", i + 1);
			
			for(int d = 0;d < node_count;d++) {
				printf("%.7f\t", final_dist(i, i, d, edges));
			}
			
			cout << endl;
		}
	}
}

string get_newick_recursive(map<pair<int, int>, Float> edges, int node, Float length) {
	int count = 0;
	
	string res = "(";
	
	for(auto i : edges) {
		int me = node;
		int other;
		
		if(i.first.first == node) {
			other = i.first.second;
		}
		else {
			other = i.first.first;
			me = i.first.second;
		}
		
		if(me == node) {
			map<pair<int, int>, Float> new_edges = edges;
			new_edges.erase(i.first);
			
			if(count != 0)
				res += ",";
			
			res += get_newick_recursive(new_edges, other, i.second);
			
			++count;
			
		}
	}
	
	//char str[2] = {node + 'A', 0};
	string str = to_string(node);
	
	if(count == 0) {
		return (string)str + ":" + to_string(length.getDouble());
	}
	else {
		if(length.getDouble() != 0)
			return res + ")" + (string)str + ":" + to_string(length.getDouble());
		else
			return res + ")" + (string)str;
	}
}
	
	
	

string get_newick(map<pair<int, int>, Float> edges, int start_node) {
	return get_newick_recursive(edges, start_node, 0) + ";";
}

int main(int argc, char *argv[]) {
	string input_file;
	bool print = false;
	bool get_switch = true;
	char cmd_switch;
	
	vector<string> s;
	
	// Parse the command line arguments
	for(int i = 1; i < argc; i++) {
		if(get_switch) {
			if(argv[i][0] != '-') {
				
#if 0
				if(file.is_open()) {
					Sequence seq(file, argv[i]);
					
					s.push_back(seq.genes["VP40"]);
					//cout << "Sequence: " << seq.genes["NP"] << endl;
				}
				else {
					cout << "Error: failed to load file '" << argv[i] << "'" << endl;
					throw -1;
				}
#endif
				string seq;
				if(!getSequence(argv[i], seq)) {
					cout << "ERROR LOADING SEQUENCE" << endl;
					throw -1;
				}
				
				s.push_back(seq);
				
				
			}
			else {
				if(argv[i][1] == 'i') {
					cmd_switch = 'i';
					get_switch = false;
				}
				else if(argv[i][1] == 'p') {
					print = true;
				}
				else {
					cout << "Error: unrecognized switch: '" << argv[i] << "'" << endl;
				}
			}
		}
		else {
			if(cmd_switch == 'i') {
				input_file = argv[i];
			}
			
			get_switch = true;
		}
	}
	
 	vector<string> s2 = {"ACGGAGA", "AGTTGACA", "ACTGACA", "CCGTTCAC", "AGGAGA", "TTTTTTTT"};
	map<pair<int, int>, Float> edges;
	
	clustalw(s, edges, print);
	
	cout << get_newick(edges, s.size()) << endl;
}