
# include <chrono>
#include <iostream>
#include <vector>
#include <queue>
#include <list>
#include <fstream>
#include <string>
#include <algorithm>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <random>
#include <functional>
#include <tuple>
#include <typeinfo>
#include<sstream>
#include <iomanip>

using namespace std;

int main() {

	string base = "C:\\Users\\Lucas\\Estudos\\PucRio\\Metaheuristicas\\Project 2\\Instances\\";
	cout << "Write jobsize machinesize file name" << endl;
	int js, ms;
	cin >> js >> ms;

	string filename = "Tai_" + to_string(js) + "_" + to_string(ms) + ".txt";	
	filename = base + filename;

	cout << "Write initial counter" << endl;
	int counter; cin >> counter;

	ifstream inst_file;
	ofstream outfile;

	inst_file.open(filename, ios::in);
	string line;

	while (getline(inst_file, line)) {
		stringstream out_file_name;
		out_file_name << "Ta" << setw(2) << setfill('0') << counter << ".txt";
		outfile.open(out_file_name.str());
		outfile << line << "\n";

		getline(inst_file, line);
		outfile << line << "\n";

		std::stringstream ss(line);
		int n;
		ss >> n;

		//Times
		getline(inst_file, line);
		outfile << line << "\n";

		for (int i = 0; i < n; i++) {
			getline(inst_file, line);
			outfile << line << "\n";
		}

		//Machines
		getline(inst_file, line);
		outfile << line << "\n";

		for (int i = 0; i < n; i++) {
			getline(inst_file, line);
			outfile << line << "\n";
		}

		outfile.close();
		counter++;
	}

	inst_file.close();
	cout << "Done." << endl;
	cin >> line;
}