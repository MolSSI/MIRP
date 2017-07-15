/*! \file
 *
 * \brief Functions related to testing single ERI
 */

#include "mirp_bin/eri_single_test.hpp"
#include <mirp/kernels/eri.h>
#include <mirp/math.h>
#include <mirp/arb_help.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>

/* Anonymous namespace for some helper functions */
namespace {

#if 0
/* Runs a Boys function test using interval arithmetic
 *
 * This just wraps the calculation of Boys function, increasing the
 * working precision until the desired target precision is reached. It
 * then handles the comparison with the reference data
 *
 * The number of failing tests is returned
 */
long boys_run_test_interval(const mirp::boys_data & data, long extra_m, long target_prec)
{
    using namespace mirp;

    const int max_m = boys_max_m(data) + extra_m;

    arb_t t_mp, vref_mp;
    arb_init(t_mp);
    arb_init(vref_mp);

    auto F_mp = std::unique_ptr<arb_t[]>(new arb_t[max_m+1]);
    for(int i = 0; i <= max_m; i++)
        arb_init(F_mp[i]);

    long nfailed = 0;

    for(const auto & it : data.values)
    {
        slong working_prec = target_prec;

        do {
            /* Increase the precision until we have achieved the target accuracy
               (with a bit of a safety factor */
            working_prec += 16;
            arb_set_str(t_mp, it.t.c_str(), working_prec);
            mirp_boys_interval(F_mp.get(), it.m+extra_m, t_mp, working_prec);


        } while(arb_rel_accuracy_bits(F_mp[it.m]) < (target_prec + 4));

        /* Round the calculated value and the reference value to the
           target precision */
        arb_set_round(F_mp[it.m], F_mp[it.m], target_prec);
        arb_set_str(vref_mp, it.value.c_str(), target_prec);

        if(!arb_overlaps(F_mp[it.m], vref_mp))
        {
            std::cout << "Entry failed test: m = " << it.m << " t = " << it.t << "\n";
            char * s1 = arb_get_str(F_mp[it.m], 1000, 0);
            char * s2 = arb_get_str(vref_mp, 1000, 0);
            std::cout << "   Calculated: " << s1 << "\n";
            std::cout << "    Reference: " << s2 << "\n\n";
            free(s1);
            free(s2);
            nfailed++;
        }
    }

    arb_clear(t_mp);
    arb_clear(vref_mp);
    for(int i = 0; i <= max_m; i++)
        arb_clear(F_mp[i]);

    return nfailed;
}


/* Runs a Boys function test using double precision
 *
 * This just wraps the calculation of Boys function and the
 * comparison of the computed data with the reference values.
 *
 * The number of failing tests is returned
 */
long boys_run_test_double(const mirp::boys_data & data, long extra_m)
{
    using namespace mirp;

    const int max_m = boys_max_m(data) + extra_m;

    std::vector<double> F_dbl(max_m+1);

    long nfailed = 0;

    for(const auto & it : data.values)
    {
        double t_dbl = std::strtod(it.t.c_str(), nullptr);
        double vref_dbl = std::strtod(it.value.c_str(), nullptr);

        mirp_boys_double(F_dbl.data(), it.m+extra_m, t_dbl);

        if(vref_dbl != F_dbl[it.m])
        {
            std::cout << "Entry failed test: m = " << it.m << " t = " << it.t << "\n";
            double reldiff = std::fabs(vref_dbl - F_dbl[it.m]);
            reldiff /= std::max(std::fabs(vref_dbl), std::fabs(F_dbl[it.m]));

            auto old_cout_prec = std::cout.precision(17);
            std::cout << "   Calculated: " << F_dbl[it.m] << "\n";
            std::cout << "    Reference: " << vref_dbl << "\n";
            std::cout << "Relative Diff: " << reldiff << "\n\n";
            std::cout.precision(old_cout_prec);
            nfailed++;
        }
    }

    return nfailed;
}
#endif

} // close anonymous namespace


namespace mirp {

eri_single_data eri_single_read_input_file(const std::string & filepath)
{
    using std::ifstream;

    ifstream infile(filepath, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for reading");

    infile.exceptions(ifstream::failbit);

    std::string line;
    eri_single_data ret;
    int ngaussians = 0;

    // We will fill in the shells here, then add to ret
    eri_single_data_entry ent;

    while(std::getline(infile, line).good())
    {
        if(line.length() == 0)
            continue;
        else if(line[0] == '#')
            ret.header += line + "\n";
        else
        {
            // read in a gaussian_single
            std::stringstream ss(line);
            ss.exceptions(std::stringstream::failbit);

            gaussian_single & g = ent.g[ngaussians];
            ss >> g.lmn[0] >> g.lmn[1] >> g.lmn[2]
               >> g.xyz[0] >> g.xyz[1] >> g.xyz[2]
               >> g.alpha;

            ngaussians++;
            if(ngaussians > 3)
            {
                ret.values.push_back(ent);
                ngaussians = 0;
            }
        }
    }

    return ret;
}


eri_single_data eri_single_read_file(const std::string & filepath)
{
    throw std::runtime_error("TODO");
/*
    using std::ifstream;

    ifstream infile(filepath, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for reading");

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
*/
}
    

void eri_single_write_file(const std::string & filepath, const eri_single_data & data)
{
    throw std::runtime_error("TODO");
/*
    using std::ofstream;

    ofstream outfile;
    outfile.open(filepath, ofstream::out | ofstream::trunc);

    if(!outfile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for writing");

    outfile.exceptions(ofstream::badbit | ofstream::failbit);

    outfile << data.header;
    outfile << data.ndigits << "\n";
    for(const auto & it : data.values)
        outfile << it.m << " " << it.t << " " << it.value << "\n";
*/
}


long eri_single_run_test(const std::string & filepath, const std::string & floattype, long extra_m, long target_prec)
{
    throw std::runtime_error("TODO");
/*
    boys_data data;
    try {
        data = boys_read_file(filepath);
    }
    catch(std::exception & ex)
    {
        std::cout << "Unable to read data from file \"" << filepath << "\": ";
        std::cout << ex.what() << "\n";
        return 2;
    }

    std::cout << "Read " << data.values.size() << " values from " << filepath << "\n";

    long nfailed = 0;

    if(floattype == "interval")
        nfailed = boys_run_test_interval(data, extra_m, target_prec);
    else if(floattype == "double")
        nfailed = boys_run_test_double(data, extra_m);
    else
    {
        std::string errmsg;
        errmsg = "Unknown floating-point type \"" + floattype + "\"";
        throw std::runtime_error(errmsg);
    }


    double percent_passed = 100.0 - (100.0 * (double)nfailed / (double)data.values.size());
    std::cout << nfailed << " / " << data.values.size() << " failed ("
              << percent_passed << "% passed)\n";

    return nfailed;
*/
}


void eri_single_create_test(const std::string & input_filepath,
                            const std::string & output_filepath,
                            long ndigits, const std::string & header)
{
    eri_single_data data = eri_single_read_input_file(input_filepath);
    data.ndigits = ndigits;
    data.header += header;

    /* Temporaries */
    arb_t xyz1[3], xyz2[3], xyz3[3], xyz4[3];
    arb_t alpha1, alpha2, alpha3, alpha4;
    mirp_init_arb_arr(xyz1, 3);
    mirp_init_arb_arr(xyz2, 3);
    mirp_init_arb_arr(xyz3, 3);
    mirp_init_arb_arr(xyz4, 3);
    arb_init(alpha1);
    arb_init(alpha2);
    arb_init(alpha3);
    arb_init(alpha4);

    /* Actual computed integral */
    arb_t integral;
    arb_init(integral);

    /* Target precision/accuracy, in bits, with a safety factor
       of 4 extra decimal digits */
    const slong target_prec = (ndigits+4) / MIRP_LOG_10_2;

    for(auto & it : data.values)
    {
        slong working_prec = target_prec;
        bool sufficient_accuracy = false;

        do {
            working_prec += 16;

            for(int i = 0; i < 3; i++)
            {
                arb_set_str(xyz1[i], it.g[0].xyz[i].c_str(), working_prec); 
                arb_set_str(xyz2[i], it.g[1].xyz[i].c_str(), working_prec); 
                arb_set_str(xyz3[i], it.g[2].xyz[i].c_str(), working_prec); 
                arb_set_str(xyz4[i], it.g[3].xyz[i].c_str(), working_prec); 
            }

            arb_set_str(alpha1, it.g[0].alpha.c_str(), working_prec);
            arb_set_str(alpha2, it.g[1].alpha.c_str(), working_prec);
            arb_set_str(alpha3, it.g[2].alpha.c_str(), working_prec);
            arb_set_str(alpha4, it.g[3].alpha.c_str(), working_prec);


            mirp_single_eri_interval(integral,
                                     it.g[0].lmn.data(), xyz1, alpha1, 
                                     it.g[1].lmn.data(), xyz2, alpha2, 
                                     it.g[2].lmn.data(), xyz3, alpha3, 
                                     it.g[3].lmn.data(), xyz4, alpha4,
                                     working_prec);

            if(arb_rel_accuracy_bits(integral) >= target_prec)
                sufficient_accuracy = true;
                
        } while(!sufficient_accuracy);

        char * s = arb_get_str(integral, ndigits, ARB_STR_NO_RADIUS);
        it.integral = s;
        free(s);
    }

    //boys_write_file(output_filepath, data);
    mirp_clear_arb_arr(xyz1, 3);
    mirp_clear_arb_arr(xyz2, 3);
    mirp_clear_arb_arr(xyz3, 3);
    mirp_clear_arb_arr(xyz4, 3);
    arb_clear(alpha1);
    arb_clear(alpha2);
    arb_clear(alpha3);
    arb_clear(alpha4);
}

    
} // closing namespace mirp

