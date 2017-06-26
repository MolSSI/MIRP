
#include "test_files.hpp"
#include <iostream>
#include <sstream>

namespace mirp {

boys_data read_boys_file(const std::string & filename)
{
    using std::ifstream;


    ifstream infile;
    infile.open(filename, ifstream::in);

    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filename + "\"");

    //infile.exceptions(ifstream::badbit | ifstream::failbit);

    std::string line;
    boys_data ret;
    bool have_ndigits = false;

    while(std::getline(infile, line).good())
    {
        if(line.length() == 0)
            continue;
        if(line[0] == '#')
            continue;

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

    return ret;
}
    
} // closing namespace mirp

