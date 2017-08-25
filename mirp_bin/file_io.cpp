/*! \file
 *
 * \brief Helper functions for reading/writing test files
 */

#include "mirp/shell.h"
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
                // n gaussians, read the integral and add the entry to the data
                ss >> ent.integral;
                data.entries.push_back(ent);
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

                // If this is an input file and we just read n gaussians,
                // add the entry to the data
                // (if not an input, there is more stuff to read)
                if(is_input && ngaussians == 4)
                {
                    data.entries.push_back(ent);
                    ent.g.clear();
                    ngaussians = 0;
                }
            }
        }
    }

    // sanity check
    for(const auto & ent : data.entries)
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
    for(const auto & ent : data.entries)
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



integral_data integral_read_file(const std::string & filepath,
                                 int n,
                                 bool is_input)
{
    using std::ifstream;

    ifstream infile(filepath, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for reading");

    std::string line;
    integral_data data;
    int ngaussians = 0;

    bool have_ndigits = false;

    // We will fill in the shells here, then add to the data
    integral_data_entry ent;

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
                // n gaussians, read the integrals and add the entry to the data
                size_t nintegral = 1;
                int dummy; // integral index
                for(const auto & it : ent.g)
                    nintegral *= MIRP_NCART(it.am)*it.ngeneral;

                ent.integrals.resize(nintegral);
                ss >> dummy >> ent.integrals[0];
                for(size_t i = 1; i < nintegral; i++)
                {
                    std::getline(infile, line);
                    std::stringstream ss2(line);
                    ss2 >> dummy >> ent.integrals[i];
                }

                data.entries.push_back(ent);
                ent.g.clear();
                ngaussians = 0;
            }
            else
            {
                gaussian_shell g;
                ss >> g.Z >> g.am >> g.nprim >> g.ngeneral
                   >> g.xyz[0] >> g.xyz[1] >> g.xyz[2];

                g.alpha.resize(g.nprim);
                g.coeff.resize(g.nprim*g.ngeneral);

                for(int i = 0; i < g.nprim; i++)
                {
                    std::getline(infile, line);
                    std::stringstream ss2(line);
                    ss2 >> g.alpha[i];
                    for(int j = 0; j < g.ngeneral; j++)
                        ss2 >> g.coeff[j*g.nprim+i];
                }

                ent.g.push_back(std::move(g));

                ngaussians++;

                // If this is an input file and we just read n gaussians,
                // add the entry to the data
                // (if not an input, there is more stuff to read)
                if(is_input && ngaussians == 4)
                {
                    data.entries.push_back(ent);
                    ent.g.clear();
                    ngaussians = 0;
                }
            }
        }
    }

    // sanity check
    for(const auto & ent : data.entries)
    {
        if(ent.g.size() != static_cast<size_t>(n))
            throw std::runtime_error("Inconsistent number of gaussians");
    }

    return data;

}


void integral_write_file(const std::string & filepath, const integral_data & data)
{
    using std::ofstream;

    ofstream outfile;
    outfile.open(filepath, ofstream::out | ofstream::trunc);

    if(!outfile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for writing");

    outfile.exceptions(ofstream::badbit | ofstream::failbit);

    outfile << data.header;
    outfile << data.ndigits << "\n";
    for(const auto & ent : data.entries)
    {
        for(const auto & g : ent.g)
        {
            outfile << g.Z << " " << g.am << " " << g.nprim << " " << g.ngeneral << " "
                    << g.xyz[0] << " " << g.xyz[1] << " " << g.xyz[2] << "\n";

            for(int i = 0; i < g.nprim; i++)
            {
                outfile << g.alpha[i];
                for(int j = 0; j < g.ngeneral; j++)
                    outfile << " " << g.coeff[j*g.nprim + i];
                outfile << "\n";
            }
        }


        size_t ncart = 1;
        for(const auto & it : ent.g)
            ncart *= MIRP_NCART(it.am) * it.ngeneral;

        for(size_t i = 0; i < ncart; i++)
            outfile << i << "  " << ent.integrals[i] << "\n";

        outfile << "\n";
    }
}



} // closing namespace mirp

