#include "global.hpp"
#include "instance_selector.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

instance_selector::instance_selector(int _ib, int _ie, string prefix)
{
    ib = _ib;
    ie = _ie;
    i_cur = ib - 1;
    inst_prefix = prefix;
}

instance_selector::instance_selector(string _unique_inst_name)
{
    is_unique_inst = true;
    this->unique_inst_name = _unique_inst_name;
}


string instance_selector::generate_path()
{
    stringstream s;
    s << inst_prefix << setw(2) << setfill('0') << i_cur << ".txt";
    return insts_folder_path + s.str();
}

string instance_selector::getNext()
{
    if (is_unique_inst)
    {
        if (generated)
            return "";
        generated = true;

        return insts_folder_path + unique_inst_name;
    }
    i_cur++;

    if (i_cur > ie || i_cur < ib)
        return "";
    auto ret = generate_path();

    return ret;
}

string instance_selector::signature()
{
    stringstream s;
    if (is_unique_inst)
    {
        s << unique_inst_name.substr(0, unique_inst_name.find('.'));
    }
    else
    {
        s << inst_prefix << setw(2) << setfill('0') << i_cur;
    }

    return s.str();
}