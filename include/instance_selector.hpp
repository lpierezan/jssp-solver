#pragma once
#include "global.hpp"

class instance_selector
{
    using string = std::string;
    /* Instance format:
	<header line>
	<n jobs> <n machines> ....
	Times
	<Time matrix: T_i_j = process time of operation j of job i.
	Machiens
	<Machine matrix: M_i_j = machine that operation j of job i uses.
	*/

    bool is_unique_inst = false, generated = false;
    string unique_inst_name;

public:
    string insts_folder_path = "./instances/";

    int ib, ie, i_cur;

    string inst_name, path, inst_prefix;

    instance_selector(int _ib = 1, int _ie = 10, string prefix = "Ta");

    instance_selector(string _unique_inst_name);

private:
    string generate_path();

public:
    string getNext();

    string signature();
};