#include "experiments.hpp"
#include <fstream>
#include <iostream>

using namespace std;

int main()
{
    
    // =========== IO config =================
    string output_path = "output/jssp_out.txt";
    std::ofstream out(output_path);

    if (out)
    {
        cout << "Output file: " << output_path << endl;
    }
    else
    {
        cout << "Unable to create output file: " << output_path << endl;
        exit(1);
    }

    auto coutbuf = std::cout.rdbuf(out.rdbuf());
    // ========================================

    // ============ Run Experiment ============
    
    experiments::do_inst_mh_solve();

    // ========================================

    std::cout.rdbuf(coutbuf);
}
