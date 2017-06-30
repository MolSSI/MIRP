#include "mirp_bin/boys_test_file.hpp"
#include <iostream>
#include <sstream>

namespace mirp {

int boys_max_m(const boys_data & data)
{  
    int maxm = 0;
    for(auto & it : data.values)
        maxm = std::max(maxm, it.m);
    return maxm;
}
    

boys_data boys_read_input_file(const std::string & filename)
{
    using std::ifstream;

    ifstream infile(filename, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filename + "\" for reading");

    infile.exceptions(ifstream::failbit);

    std::string line;
    boys_data ret;

    while(std::getline(infile, line).good())
    {
        if(line.length() == 0)
            continue;
        else if(line[0] == '#')
            ret.header += line + "\n";
        else
        {
            std::stringstream ss(line);
            ss.exceptions(std::stringstream::failbit);

            boys_data_entry ent;
            ss >> ent.m >> ent.t;
            ret.values.push_back(ent);
        }
    }

    return ret;
}


boys_data boys_read_file(const std::string & filename)
{
    using std::ifstream;

    ifstream infile(filename, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filename + "\" for reading");

    infile.exceptions(ifstream::failbit);

    std::string line;
    boys_data ret;
    bool have_ndigits = false;

    while(std::getline(infile, line).good())
    {
        if(line.length() == 0)
            continue;
        else if(line[0] == '#')
            ret.header += line + "\n";
        else
        {
            std::stringstream ss(line);
            ss.exceptions(std::stringstream::failbit);

            if(!have_ndigits)
            {
                ss >> ret.ndigits;
                have_ndigits = true;    
            }
            else
            {
                boys_data_entry ent;
                ss >> ent.m >> ent.t >> ent.value;
                ret.values.push_back(ent);
            }
        }
    }

    return ret;
}
    

void boys_write_file(const std::string & filename, const boys_data & data)
{
    using std::ofstream;

    ofstream outfile;
    outfile.open(filename, ofstream::out | ofstream::trunc);

    if(!outfile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filename + "\" for writing");

    outfile.exceptions(ofstream::badbit | ofstream::failbit);

    outfile << data.header;
    outfile << data.ndigits << "\n";
    for(const auto & it : data.values)
        outfile << it.m << " " << it.t << " " << it.value << "\n";
}
    
    
} // closing namespace mirp

