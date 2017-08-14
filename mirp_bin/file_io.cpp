/*! \file
 *
 * \brief Helper functions for reading/writing test files
 */

#include "mirp_bin/file_io.hpp"
#include "mirp_bin/data_entry.hpp"
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace mirp {

integral_single_data integral_single_read_file(const std::string & filepath, int n, bool is_input)
{
    using std::ifstream;

    ifstream infile(filepath, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for reading");

    std::string line;
    integral_single_data data;
    int ngaussians = 0;

    bool have_ndigits = false;

    // We will fill in the shells here, then add to the data
    integral_single_data_entry ent;

    while(std::getline(infile, line).good())
    {
        if(line.length() == 0)
            continue;
        else if(line[0] == '#')
            data.header += line + "\n";
        else
        {
            std::stringstream ss(line);
            ss.exceptions(std::stringstream::failbit);

            // Read in ndigits if this isn't an input file
            // (and we haven't read it yet)
            if(!is_input && !have_ndigits)
            {
                ss >> data.ndigits;
                have_ndigits = true;
            }
            else if(!is_input && ngaussians == n)
            {
                // If this isn't an input file, and we have already read
                // 4 gaussians, read the integral and add the entry to the data
                ss >> ent.integral;
                data.values.push_back(ent);
                ent.g.clear();
                ngaussians = 0;
            }
            else
            {
                gaussian_single g;
                ss >> g.lmn[0] >> g.lmn[1] >> g.lmn[2]
                   >> g.xyz[0] >> g.xyz[1] >> g.xyz[2]
                   >> g.alpha;
                ent.g.push_back(std::move(g));

                ngaussians++;

                // If this is an input file and we just read 4 gaussians,
                // add the entry to the data
                if(is_input && ngaussians == 4)
                {
                    data.values.push_back(ent);
                    ent.g.clear();
                    ngaussians = 0;
                }
            }
        }
    }

    // sanity check
    for(const auto & ent : data.values)
    {
        if(ent.g.size() != static_cast<size_t>(n))
            throw std::runtime_error("Inconsistent number of gaussians");
    }

    return data;
}


void integral_single_write_file(const std::string & filepath, const integral_single_data & data)
{
    using std::ofstream;

    ofstream outfile;
    outfile.open(filepath, ofstream::out | ofstream::trunc);

    if(!outfile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for writing");

    outfile.exceptions(ofstream::badbit | ofstream::failbit);

    outfile << data.header;
    outfile << data.ndigits << "\n";
    for(const auto & ent : data.values)
    {
        for(const auto & g : ent.g)
        {
            outfile << g.lmn[0] << " " << g.lmn[1] << " " << g.lmn[2] << " "
                    << g.xyz[0] << " " << g.xyz[1] << " " << g.xyz[2] << " "
                    << g.alpha << "\n";
        }

        outfile << ent.integral << "\n\n";
    }
}


} // closing namespace mirp

