/*! \file
 *
 * \brief Helper functions for reading/writing test files
 */

#include "mirp_bin/testfile_io.hpp"
#include "mirp_bin/data_entry.hpp"
#include "mirp_bin/test_common.hpp"
#include <mirp/shell.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace mirp {

integral_single_data testfile_read_integral_single(const std::string & filepath, int n, bool is_input)
{
    using std::ifstream;

    if(n <= 0)
        throw std::logic_error("Cannot read negative or zero gaussians from a file");


    // Used in errors
    std::stringstream sserr;
    sserr << "Error reading file " << filepath << ": ";


    ifstream infile(filepath);
    if(!infile.is_open())
    {
        sserr << "Cannot open file";
        throw std::runtime_error(sserr.str());
    }

    integral_single_data data;
    size_t nentry = 0;

    // read in the header comments
    while(infile.peek() == '#')
    {
        std::string line;
        std::getline(infile, line);
        data.header += line + "\n";
    }

    // Read the expected number of entries
    infile >> nentry;

    // Read in the number of digits and the working prec
    if(!is_input)
    {
        infile >> data.ndigits >> data.working_prec;
        if(!infile.good())
        {
            sserr << "Error reading metadata (nentry, ndigits, working_prec)";
            throw std::runtime_error(sserr.str());
        }
    }

    // read the actual data
    while(infile.good())
    {
        // check if there is more data
        if(!file_skip(infile, '#'))
            break;

        // We will fill in the shells here, then add to the data
        integral_single_data_entry ent;

        // read in n gaussians
        for(int i = 0; i < n; i++)
        {
            if(!file_skip(infile, '#'))
            {
                if(infile.eof())
                {
                    sserr << "Unexpected end of file while reading gaussian " << i << "/" << n
                          << " for entry " << data.entries.size() << "\n";
                    throw std::runtime_error(sserr.str());
                }
                else
                {
                    sserr << "Error while reading gaussian " << i << "/" << n
                          << " for entry " << data.entries.size() << "\n";
                    throw std::runtime_error(sserr.str());
                }
            } 


            gaussian_single_str g;

            infile >> g.lmn[0] >> g.lmn[1] >> g.lmn[2]
                   >> g.xyz[0] >> g.xyz[1] >> g.xyz[2]
                   >> g.alpha;
 
            if(infile.bad() || infile.fail())
            {
                sserr << "Error while reading gaussian " << i << "/" << n
                      << " for entry " << data.entries.size() << "\n";
                throw std::runtime_error(sserr.str());
            }
                
            ent.g.push_back(std::move(g));
        }

        // If this is not an input file, read the integral also
        if(!is_input)
        {
            // the first one reads in the rest of the line (which should
            // be empty). The second will contain the integral
            file_skip(infile, '#');
            std::getline(infile, ent.integral);

            if(infile.bad() || infile.fail())
            {
                sserr << "Error while integral for entry " << data.entries.size() << "\n";
                throw std::runtime_error(sserr.str());
            }
        }

        // add the entry
        data.entries.push_back(std::move(ent));
    }

    if(data.entries.size() != nentry)
    {
        sserr << "Number of entries not consistent: Expected " << nentry
              << " but got " << data.entries.size() << "\n";
        throw std::runtime_error(sserr.str());
    }

    std::cout << "Read " << data.entries.size() << " values from " << filepath << "\n";
    return data;
}


void testfile_write_integral_single(const std::string & filepath, const integral_single_data & data)
{
    using std::ofstream;

    ofstream outfile;
    outfile.open(filepath, ofstream::out | ofstream::trunc);

    if(!outfile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for writing");

    outfile.exceptions(ofstream::badbit | ofstream::failbit);

    outfile << data.header;
    outfile << data.entries.size() << " "
            << data.ndigits << " "
            << data.working_prec << "\n";

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



integral_data testfile_read_integral(const std::string & filepath,
                                     int n,
                                     bool is_input)
{
    using std::ifstream;

    if(n <= 0)
        throw std::logic_error("Cannot read negative or zero gaussians from a file");


    // Used in errors
    std::stringstream sserr;
    sserr << "Error reading file " << filepath << ": ";


    ifstream infile(filepath, ifstream::in);
    if(!infile.is_open())
    {
        sserr << "Cannot open file";
        throw std::runtime_error(sserr.str());
    }

    integral_data data;
    size_t nentry = 0;

    // read in the header comments
    while(infile.peek() == '#')
    {
        std::string line;
        std::getline(infile, line);
        data.header += line + "\n";
    }


    // Read the expected number of entries
    infile >> nentry;

    // Read in the number of digits and the working prec
    if(!is_input)
    {
        infile >> data.ndigits >> data.working_prec;
        if(!infile.good())
        {
            sserr << "Error reading metadata (nentry, ndigits, working_prec)";
            throw std::runtime_error(sserr.str());
        }
    }


    // read the actual data
    while(infile.good())
    {
        // check if there is more data
        if(!file_skip(infile, '#'))
            break;

        // We will fill in the shells here, then add to the data
        integral_data_entry ent;

        // read in n gaussians
        for(int i = 0; i < n; i++)
        {
            if(!file_skip(infile, '#'))
            {
                if(infile.eof())
                {
                    sserr << "Unexpected end of file while reading gaussian " << i << "/" << n
                          << " for entry " << data.entries.size() << "\n";
                    throw std::runtime_error(sserr.str());
                }
                else
                {
                    sserr << "Error while reading gaussian " << i << "/" << n
                          << " for entry " << data.entries.size() << "\n";
                    throw std::runtime_error(sserr.str());
                }
            } 

            gaussian_shell_str g;
            infile >> g.Z >> g.am >> g.nprim >> g.ngeneral
                   >> g.xyz[0] >> g.xyz[1] >> g.xyz[2];

            if(infile.bad() || infile.fail())
            {
                sserr << "Error while reading gaussian " << i << "/" << n
                      << " for entry " << data.entries.size() << "\n";
                throw std::runtime_error(sserr.str());
            }

            if(!file_skip(infile, '#'))
            {
                sserr << "Unexpected EOF or error after reading info for gaussian " << i << "/" << n
                      << " for entry " << data.entries.size() << "\n";
                throw std::runtime_error(sserr.str());
            }

            g.alpha.resize(g.nprim);
            g.coeff.resize(g.nprim*g.ngeneral);

            for(int j = 0; j < g.nprim; j++)
            {
                infile >> g.alpha[j];
                for(int k = 0; k < g.ngeneral; k++)
                    infile >> g.coeff[k*g.nprim+j];
            }

            if(infile.bad() || infile.fail())
            {
                sserr << "Error while reading exponents and coefficients for gaussian " << i << "/" << n
                      << " for entry " << data.entries.size() << "\n";
                throw std::runtime_error(sserr.str());
            }

            ent.g.push_back(std::move(g));
        }


        // read in the integral values 
        if(!is_input)
        {
            // Read and discard the rest of the line
            infile.ignore(max_length, '\n');

            // number of integrals we should be reading
            size_t nintegral = 1;
            for(const auto & it : ent.g)
                nintegral *= MIRP_NCART(it.am)*it.ngeneral;

            ent.integrals.resize(nintegral);

            // If this isn't an input file, read the integrals
            int dummy; // integral index (not needed)

            for(size_t i = 0; i < nintegral; i++)
            {
                infile >> dummy;
                std::getline(infile, ent.integrals[i]);
            }

            if(infile.bad() || infile.fail())
            {
                sserr << "Error while reading integrals for entry " << data.entries.size() << "\n";
                throw std::runtime_error(sserr.str());
            }
        }

        // add to the entries
        data.entries.push_back(std::move(ent));
    }

    if(data.entries.size() != nentry)
    {
        sserr << "Number of entries not consistent: Expected " << nentry
              << " but got " << data.entries.size() << "\n";
        throw std::runtime_error(sserr.str());
    }

    std::cout << "Read " << data.entries.size() << " values from " << filepath << "\n";
    return data;
}


void testfile_write_integral(const std::string & filepath, const integral_data & data)
{
    using std::ofstream;

    ofstream outfile;
    outfile.open(filepath, ofstream::out | ofstream::trunc);

    if(!outfile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for writing");

    outfile.exceptions(ofstream::badbit | ofstream::failbit);

    outfile << data.header;
    outfile << data.entries.size() << " "
            << data.ndigits << " "
            << data.working_prec << "\n";

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

