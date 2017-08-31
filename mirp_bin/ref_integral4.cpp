/*! \file
 *
 * \brief Functions related to four-center integral reference files
 */

#include "mirp_bin/ref_integral4.hpp"
#include "mirp_bin/reffile_io.hpp"
#include "mirp_bin/test_common.hpp"

#include <mirp/shell.h>

#include <fstream>

namespace mirp {

void integral4_create_reference(const std::string & xyz_filepath,
                                const std::string & basis_filepath,
                                const std::string & output_filepath,
                                const std::string & header,
                                cb_integral4_exact cb)
{
    std::vector<gaussian_shell> shells = read_construct_basis(xyz_filepath, basis_filepath);

    const size_t nshell = shells.size();

    std::ofstream fs(output_filepath);
    if(!fs.is_open())
        throw std::runtime_error("Error opening output file for writing");

    fs << header << "\n";
    fs << std::hexfloat;
    reffile_write_basis(shells, fs);

    std::vector<double> integrals;

    for(size_t p = 0; p < nshell; p++)
    for(size_t r = 0; r < nshell; r++)
    for(size_t q = 0; q <= p; q++)
    for(size_t s = 0; s <= r; s++)
    {
        const size_t pq = (p*(p+1))/2 + q;
        const size_t rs = (r*(r+1))/2 + r;

        if(pq < rs)
            continue;

        const auto & s1 = shells[p];
        const auto & s2 = shells[q];
        const auto & s3 = shells[r];
        const auto & s4 = shells[s];


        const size_t ncart = MIRP_NCART4(s1.am, s2.am, s3.am, s4.am);
        const size_t ngen = s1.ngeneral * s2.ngeneral * s3.ngeneral * s4.ngeneral;
        const size_t nintegrals = ncart * ngen;

        integrals.resize(nintegrals);

        cb(integrals.data(),
           s1.am, s1.xyz.data(), s1.nprim, s1.ngeneral, s1.alpha.data(), s1.coeff.data(),
           s2.am, s2.xyz.data(), s2.nprim, s2.ngeneral, s2.alpha.data(), s2.coeff.data(),
           s3.am, s3.xyz.data(), s3.nprim, s3.ngeneral, s3.alpha.data(), s3.coeff.data(),
           s4.am, s4.xyz.data(), s4.nprim, s4.ngeneral, s4.alpha.data(), s4.coeff.data());

        fs << p << " " << q << " " << r << " " << s;
        for(size_t i = 0; i < nintegrals; i++)
            fs << " " <<  integrals[i];
        fs << "\n";
    }
}


long integral4_test_reference(const std::string & ref_filepath,
                              cb_integral4_exact cb)
{
    std::ifstream fs(ref_filepath);
    if(!fs.is_open())
        throw std::runtime_error("Error opening input file");


    file_skip_comments(fs, '#');
    std::vector<gaussian_shell> shells = reffile_read_basis(fs);

    std::string line;
    long nfailed = 0;
    long ncomputed = 0;

    //#pragma omp parallel for schedule(dynamic) collapse(2)
    while(fs.good())
    {
        size_t p, q, r, s;
        fs >> p >> q >> r >> s;

        if(!fs.good())
            break;

        const gaussian_shell & s1 = shells.at(p);
        const gaussian_shell & s2 = shells.at(q);
        const gaussian_shell & s3 = shells.at(r);
        const gaussian_shell & s4 = shells.at(s);

        const size_t ncart = MIRP_NCART4(s1.am, s2.am, s3.am, s4.am);
        const size_t ngen = s1.ngeneral * s2.ngeneral * s3.ngeneral * s4.ngeneral;
        const size_t nintegrals = ncart * ngen;

        std::vector<double> integrals(nintegrals);
        std::vector<double> integrals_file(nintegrals);
        std::string strtmp;

        for(size_t i = 0; i < nintegrals; i++)
        {
            fs >> strtmp;
            integrals_file[i] = std::strtod(strtmp.c_str(), nullptr);
        }

        cb(integrals.data(),
           s1.am, s1.xyz.data(), s1.nprim, s1.ngeneral, s1.alpha.data(), s1.coeff.data(),
           s2.am, s2.xyz.data(), s2.nprim, s2.ngeneral, s2.alpha.data(), s2.coeff.data(),
           s3.am, s3.xyz.data(), s3.nprim, s3.ngeneral, s3.alpha.data(), s3.coeff.data(),
           s4.am, s4.xyz.data(), s4.nprim, s4.ngeneral, s4.alpha.data(), s4.coeff.data());

        for(size_t i = 0; i < nintegrals; i++)
        {
            if(integrals[i] != integrals_file[i])
            {
                printf("Failed entry: ( %2d %2d | %2d %2d ) %4lu %4lu %4lu %4lu %7lu  -> %26.18e %26.18e\n",
                     s1.am, s2.am, s3.am, s4.am, 
                     p, q, r, s, i,
                     integrals[i], integrals_file[i]);
                nfailed++;
            }
        }
        ncomputed += nintegrals;
    }

    print_results(nfailed, ncomputed);

    return nfailed;
}

} // close namespace mirp

