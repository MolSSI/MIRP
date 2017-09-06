/*! \file
 *
 * \brief Helper functions for reading/writing reference data files
 */

#include "mirp_bin/reffile_io.hpp"

namespace mirp {


double read_double(std::istream & fs)
{
    std::string tmp;
    fs >> tmp;
    return std::strtod(tmp.c_str(), nullptr);
}


void reffile_write_basis(const std::vector<gaussian_shell> & shells, std::ostream & fs)
{
    fs << shells.size() << "\n";

    for(const auto & s : shells)
    {
        fs << s.am << " " << s.nprim << " " << s.ngeneral << " "
            << s.xyz[0] << " " << s.xyz[1] << " " << s.xyz[2] << "\n";
        fs << s.alpha[0];
        for(int a = 1; a < s.nprim; a++)
            fs << " " << s.alpha[a];
        fs << "\n";
        for(int n = 0; n < s.ngeneral; n++)
        {
            fs << s.coeff[n*s.nprim+0];
            for(int a = 1; a < s.nprim; a++)
                fs << " " << s.coeff[n*s.nprim+a];
            fs << "\n";
        }
    }
}


std::vector<gaussian_shell> reffile_read_basis(std::istream & fs)
{
    size_t nshell;
    fs >> nshell;

    std::vector<gaussian_shell> shells(nshell);

    for(size_t i = 0; i < nshell; i++)
    {
        auto & s = shells[i];
        fs >> s.am >> s.nprim >> s.ngeneral;
        
        s.xyz[0] = read_double(fs);
        s.xyz[1] = read_double(fs);
        s.xyz[2] = read_double(fs);

        s.alpha.resize(s.nprim);
        s.coeff.resize(s.nprim*s.ngeneral);

        for(int j = 0; j < s.nprim; j++)
            s.alpha[j] = read_double(fs);

        for(int n = 0; n < s.ngeneral; n++)
        {
            for(int j = 0; j < s.nprim; j++)
                s.coeff[n*s.nprim+j] = read_double(fs);
        }
    }

    return shells;
}


} // closing namespace mirp

