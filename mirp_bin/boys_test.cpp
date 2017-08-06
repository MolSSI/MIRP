/*! \file
 *
 * \brief Functions related to testing the Boys function
 */

#include "mirp_bin/boys_test.hpp"
#include <mirp/kernels/boys.h>
#include <mirp/math.h>
#include <mirp/arb_help.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <memory>

/* Anonymous namespace for some helper functions */
namespace {

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

    long nfailed = 0;

    const int max_m = boys_max_m(data) + extra_m;

    arb_t t_mp, vref_mp;
    arb_init(t_mp);
    arb_init(vref_mp);

    auto F_mp = std::unique_ptr<arb_t[]>(new arb_t[max_m+1]);
    mirp_init_arb_arr(F_mp.get(), max_m+1);

    for(const auto & ent : data.values)
    {
        slong working_prec = target_prec;

        do {
            /* Increase the precision until we have achieved the target accuracy
               (with a bit of a safety factor */
            working_prec += 16;
            arb_set_str(t_mp, ent.t.c_str(), working_prec);
            mirp_boys_interval(F_mp.get(), ent.m+extra_m, t_mp, working_prec);


        } while(arb_rel_accuracy_bits(F_mp[ent.m]) < (target_prec + 4));

        /* Round the calculated value and the reference value to the
           target precision */
        arb_set_round(F_mp[ent.m], F_mp[ent.m], target_prec);
        arb_set_str(vref_mp, ent.value.c_str(), target_prec);

        if(!arb_overlaps(F_mp[ent.m], vref_mp))
        {
            std::cout << "Entry failed test: m = " << ent.m << " t = " << ent.t << "\n";
            char * s1 = arb_get_str(F_mp[ent.m], 1000, 0);
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
    mirp_clear_arb_arr(F_mp.get(), max_m+1);

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

    long nfailed = 0;

    const int max_m = boys_max_m(data) + extra_m;
    std::vector<double> F_dbl(max_m+1);

    for(const auto & ent : data.values)
    {
        double t_dbl = std::strtod(ent.t.c_str(), nullptr);
        double vref_dbl = std::strtod(ent.value.c_str(), nullptr);

        mirp_boys_double(F_dbl.data(), ent.m+extra_m, t_dbl);

        if(vref_dbl != F_dbl[ent.m])
        {
            std::cout << "Entry failed test: m = " << ent.m << " t = " << ent.t << "\n";
            double reldiff = std::fabs(vref_dbl - F_dbl[ent.m]);
            reldiff /= std::max(std::fabs(vref_dbl), std::fabs(F_dbl[ent.m]));

            auto old_cout_prec = std::cout.precision(17);
            std::cout << "   Calculated: " << F_dbl[ent.m] << "\n";
            std::cout << "    Reference: " << vref_dbl << "\n";
            std::cout << "Relative Diff: " << reldiff << "\n\n";
            std::cout.precision(old_cout_prec);
            nfailed++;
        }
    }

    return nfailed;
}

} // close anonymous namespace


namespace mirp {

int boys_max_m(const boys_data & data)
{  
    int max_m = 0;
    for(auto & ent : data.values)
        max_m = std::max(max_m, ent.m);
    return max_m;
}
    

boys_data boys_read_file(const std::string & filepath, bool is_input)
{
    using std::ifstream;

    ifstream infile(filepath, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for reading");

    std::string line;
    boys_data data;
    bool have_ndigits = false;

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

            if(!is_input && !have_ndigits)
            {
                ss >> data.ndigits;
                have_ndigits = true;    
            }
            else
            {
                boys_data_entry ent;
                ss >> ent.m >> ent.t;

                if(!is_input)
                    ss  >> ent.value;

                data.values.push_back(ent);
            }
        }
    }

    return data;
}
    

void boys_write_file(const std::string & filepath, const boys_data & data)
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
        outfile << ent.m << " " << ent.t << " " << ent.value << "\n";
}


long boys_run_test(const std::string & filepath, const std::string & floattype, long extra_m, long target_prec)
{
    boys_data data;
    try {
        data = boys_read_file(filepath, false);
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
        std::string err;
        err = "Unknown floating-point type \"" + floattype + "\"";
        throw std::runtime_error(err);
    }


    double percent_passed = 100.0 - (100.0 * (double)nfailed / (double)data.values.size());
    std::cout << nfailed << " / " << data.values.size() << " failed ("
              << percent_passed << "% passed)\n";

    return nfailed;
}


void boys_create_test(const std::string & input_filepath,
                      const std::string & output_filepath,
                      long ndigits, const std::string & header)
{
    boys_data data = boys_read_file(input_filepath, true);
    data.ndigits = ndigits;
    data.header += header;

    const int max_m = boys_max_m(data);

    arb_t t_mp;
    auto F_mp = std::unique_ptr<arb_t[]>(new arb_t[max_m+1]);
    arb_init(t_mp);
    mirp_init_arb_arr(F_mp.get(), max_m+1);

    /* Target precision/accuracy, in bits, with a safety factor
       of 4 extra decimal digits */
    const slong target_prec = (ndigits+4) / MIRP_LOG_10_2;

    for(auto & ent : data.values)
    {
        slong working_prec = target_prec;
        bool sufficient_accuracy = false;
        
        do {
            working_prec += 16;

            arb_set_str(t_mp, ent.t.c_str(), working_prec);
            mirp_boys_interval(F_mp.get(), ent.m, t_mp, working_prec); 

            if(arb_rel_accuracy_bits(F_mp[ent.m]) >= target_prec)
                sufficient_accuracy = true;
                
        } while(!sufficient_accuracy);

        char * s = arb_get_str(F_mp[ent.m], ndigits, ARB_STR_NO_RADIUS);
        ent.value = s;
        free(s);
    }

    boys_write_file(output_filepath, data);
    arb_clear(t_mp);
    mirp_clear_arb_arr(F_mp.get(), max_m+1);
}

    
} // closing namespace mirp

