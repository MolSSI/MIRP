/*! \file
 *
 * \brief Functions related to testing the Boys function
 */

#include "mirp_bin/boys_test.hpp"
#include "mirp_bin/test_common.hpp"
#include <mirp/kernels/boys.h>
#include <mirp/math.h>
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
 *
 * \todo This function is not exception safe
 */
long boys_run_test(const mirp::boys_data & data, long extra_m, long target_prec)
{
    using namespace mirp;

    /* Number of binary digits contained in the reference value strings
       (and the number of binary digits of accuracy, taking into account
       that the number printed is +/- 1 decimal ulp) */
    long integral_bits = data.ndigits / MIRP_LOG_10_2;
    long round_bits = (data.ndigits-1) / MIRP_LOG_10_2;

    long nfailed = 0;

    const int max_m = boys_max_m(data) + extra_m;

    arb_t vref_mp;
    arb_init(vref_mp);

    arb_ptr F_mp = _arb_vec_init(max_m+1);

    for(const auto & ent : data.values)
    {
        /* 16 extra bits (~4-5 decimal digits) for safety */
        mirp_boys_target_str(F_mp, ent.m + extra_m, ent.t.c_str(), target_prec+16);

        /* Convert the reference to arb_t. The reference is printed to +/- 1 ulp (decimal). */
        arb_set_str(vref_mp, ent.value.c_str(), integral_bits + 16);
        arf_mag_add_ulp(arb_radref(vref_mp), arb_radref(vref_mp), arb_midref(vref_mp), round_bits);
        arb_set_round(vref_mp, vref_mp, target_prec);

        /* Rounding the reference value to the target precision results in
         * an interval. Does that interval contain our (more precise) result? */
        if(!arb_contains(vref_mp, F_mp + ent.m))
        {
            std::cout << "Entry failed test: m = " << ent.m << " t = " << ent.t << "\n";
            char * s1 = arb_get_str(F_mp + ent.m, 2*data.ndigits, 0);
            char * s2 = arb_get_str(vref_mp, 2*data.ndigits, 0);
            std::cout << "   Calculated: " << s1 << "\n";
            std::cout << "    Reference: " << s2 << "\n\n";
            free(s1);
            free(s2);
            nfailed++;
        }
    }

    arb_clear(vref_mp);
    _arb_vec_clear(F_mp, max_m+1);

    return nfailed;
}


/* Runs a Boys function test using double precision
 *
 * This just wraps the calculation of Boys function and the
 * comparison of the computed data with the reference values.
 *
 * The number of failing tests is returned
 */
long boys_run_test_d(const mirp::boys_data & data, long extra_m)
{
    using namespace mirp;
    using detail::almost_equal;

    long nfailed = 0;

    const int max_m = boys_max_m(data) + extra_m;
    std::vector<double> F_dbl(max_m+1);

    for(const auto & ent : data.values)
    {
        double t_dbl = std::strtod(ent.t.c_str(), nullptr);
        double vref_dbl = std::strtod(ent.value.c_str(), nullptr);

        mirp_boys_d(F_dbl.data(), ent.m+extra_m, t_dbl);

        if(!almost_equal(F_dbl[ent.m], vref_dbl, 1e-13))
        {
            double reldiff = std::fabs(vref_dbl - F_dbl[ent.m]);
            reldiff /= std::fmax(std::fabs(vref_dbl), std::fabs(F_dbl[ent.m]));

            std::cout << "Entry failed test: m = " << ent.m << " t = " << ent.t << "\n";
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


/* Runs a Boys function test using 'exact' double precision
 *
 * This is just a simple test of the wrapper. The function compares
 * the 'exact' version with the result of an interval version
 * with a high working precision.
 *
 * This, therefore, just ensures that the wrappers are written correctly.
 */
long boys_run_test_exact(const mirp::boys_data & data, long extra_m)
{
    using namespace mirp;

    long nfailed = 0;

    const int max_m = boys_max_m(data) + extra_m;
    std::vector<double> F_dbl(max_m+1);

    /* For comparison */
    arb_t t_mp;
    arb_init(t_mp);

    arb_ptr F_mp = _arb_vec_init(max_m+1);

    for(const auto & ent : data.values)
    {
        double t_dbl = std::strtod(ent.t.c_str(), nullptr);

        /* Compute using the "exact" code */
        mirp_boys_exact(F_dbl.data(), ent.m+extra_m, t_dbl);

        /* Compute using the interval arithmetic code */
        /* 256 bits should be enough for testing... */
        arb_set_d(t_mp, t_dbl);
        mirp_boys(F_mp, ent.m+extra_m, t_mp, 256);

        /* Make sure we really didn't lose a whole bunch of precision */
        if(arb_rel_accuracy_bits(F_mp + ent.m) < 64)
            throw std::logic_error("Not enough bits in testing boys exact function. Contact the developer");

        double vref_dbl = std::strtod(ent.value.c_str(), nullptr);
        double vref2_dbl = arf_get_d(arb_midref(F_mp + ent.m), ARF_RND_NEAR);

        if(F_dbl[ent.m] != vref_dbl && F_dbl[ent.m] != vref2_dbl)
        {
            std::cout << "Entry failed test: m = " << ent.m << " t = " << ent.t << "\n";
            auto old_cout_prec = std::cout.precision(17);
            std::cout << "     Calculated: " << F_dbl[ent.m] << "\n";
            std::cout << "      Reference: " << vref2_dbl << "\n";
            std::cout << " File Reference: " << vref_dbl << "\n\n";
            std::cout.precision(old_cout_prec);
            nfailed++;
        }
    }

    arb_clear(t_mp);
    _arb_vec_clear(F_mp, max_m+1);

    return nfailed;
}

} // close anonymous namespace


namespace mirp {

int boys_max_m(const boys_data & data)
{
    int max_m = 0;
    for(auto & ent : data.values)
        max_m = std::fmax(max_m, ent.m);
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


long boys_run_test_main(const std::string & filepath,
                        const std::string & floattype,
                        long extra_m,
                        long target_prec)
{
    boys_data data = boys_read_file(filepath, false);

    std::cout << "Read " << data.values.size() << " values from " << filepath << "\n";

    long nfailed = 0;

    if(floattype == "interval")
        nfailed = boys_run_test(data, extra_m, target_prec);
    else if(floattype == "exact")
        nfailed = boys_run_test_exact(data, extra_m);
    else if(floattype == "double")
        nfailed = boys_run_test_d(data, extra_m);
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
    arb_ptr F_mp = _arb_vec_init(max_m+1);
    arb_init(t_mp);

    /* Target precision/accuracy, in bits, with a safety factor
       of 5 extra decimal digits */
    const slong target_prec = (ndigits+5) / MIRP_LOG_10_2;

    for(auto & ent : data.values)
    {
        mirp_boys_target_str(F_mp, ent.m, ent.t.c_str(), target_prec);

        char * s = arb_get_str(F_mp + ent.m, ndigits, ARB_STR_NO_RADIUS);
        ent.value = s;
        free(s);
    }

    boys_write_file(output_filepath, data);
    arb_clear(t_mp);
    _arb_vec_clear(F_mp, max_m+1);
}


} // closing namespace mirp

