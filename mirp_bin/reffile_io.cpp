/*! \file
 *
 * \brief Helper functions for reading/writing reference data files
 */

#include "mirp_bin/reffile_io.hpp"
#include "test_common.hpp"
#include <cstdio>
#include <stdexcept>

namespace mirp {


double read_hexdouble(std::istream & fs)
{
    std::string tmp;
    fs >> tmp;
    return std::strtod(tmp.c_str(), nullptr);
}


void write_hexdouble(double d, std::ostream & fs)
{
    char buf[64];
    int n = snprintf(buf, 64, "%a", d);

    if(n < 0 || n >= 64)
        throw std::logic_error("Buffer not big enough in write_hexdouble");

    fs << buf;
}


void reffile_write_basis(const std::vector<gaussian_shell> & shells, std::ostream & fs)
{
    fs << shells.size() << "\n";

    for(const auto & s : shells)
    {
        fs << s.am << " " << s.nprim << " " << s.ngeneral << " ";

        write_hexdouble(s.xyz[0], fs);
        fs << " ";
        write_hexdouble(s.xyz[1], fs);
        fs << " ";
        write_hexdouble(s.xyz[2], fs);
        fs << "\n"; 

        write_hexdouble(s.alpha[0], fs);

        for(int a = 1; a < s.nprim; a++)
        {
            fs << " ";
            write_hexdouble(s.alpha[a], fs);
        }

        fs << "\n";
        for(int n = 0; n < s.ngeneral; n++)
        {
            write_hexdouble(s.coeff[n*s.nprim+0], fs);

            for(int a = 1; a < s.nprim; a++)
            {
                fs << " ";
                write_hexdouble(s.coeff[n*s.nprim+a], fs);
            }

            fs << "\n";
        }
    }
}


std::vector<gaussian_shell> reffile_read_basis(std::istream & fs)
{
    size_t nshell;

    file_skip(fs, '#');
    fs >> nshell;

    std::vector<gaussian_shell> shells(nshell);

    for(size_t i = 0; i < nshell; i++)
    {
        auto & s = shells[i];

        file_skip(fs, '#');
        fs >> s.am >> s.nprim >> s.ngeneral;

        file_skip(fs, '#');
        s.xyz[0] = read_hexdouble(fs);
        s.xyz[1] = read_hexdouble(fs);
        s.xyz[2] = read_hexdouble(fs);

        s.alpha.resize(s.nprim);
        s.coeff.resize(s.nprim*s.ngeneral);

        file_skip(fs, '#');
        for(int j = 0; j < s.nprim; j++)
            s.alpha[j] = read_hexdouble(fs);

        for(int n = 0; n < s.ngeneral; n++)
        {
            file_skip(fs, '#');
            for(int j = 0; j < s.nprim; j++)
                s.coeff[n*s.nprim+j] = read_hexdouble(fs);
        }
    }

    return shells;
}


} // closing namespace mirp

