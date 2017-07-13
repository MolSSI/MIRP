/*! \file
 *
 * \brief Functions related to testing the Boys function
 */

#include "mirp_bin/boys_test.hpp"
#include <mirp/kernels/boys.h>
#include <mirp/math.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>

/* Anonymous namespace for some helper functions */
namespace {

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

} // close anonymous namespace


namespace mirp {

int boys_max_m(const boys_data & data)
{  
    int maxm = 0;
    for(auto & it : data.values)
        maxm = std::max(maxm, it.m);
    return maxm;
}
    

boys_data boys_read_input_file(const std::string & filepath)
{
    using std::ifstream;

    ifstream infile(filepath, ifstream::in);
    if(!infile.is_open())
        throw std::runtime_error(std::string("Unable to open file \"") + filepath + "\" for reading");

    infile.exceptions(ifstream::failbit);

    std::string line;
    boys_data ret;

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

            boys_data_entry ent;
            ss >> ent.m >> ent.t;
            ret.values.push_back(ent);
        }
    }

    return ret;
}


boys_data boys_read_file(const std::string & filepath)
{
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
    for(const auto & it : data.values)
        outfile << it.m << " " << it.t << " " << it.value << "\n";
}


long boys_run_test(const std::string & filename, const std::string & floattype, long extra_m, long target_prec)
{
    boys_data data;
    try {
        data = boys_read_file(filename);
    }
    catch(std::exception & ex)
    {
        std::cout << "Unable to read data from file \"" << filename << "\": ";
        std::cout << ex.what() << "\n";
        return 2;
    }

    std::cout << "Read " << data.values.size() << " values from " << filename << "\n";

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
}


void boys_create_test(const std::string & infile, const std::string & outfile, long ndigits, const std::string & header)
{
    boys_data data = boys_read_input_file(infile);
    data.ndigits = ndigits;
    data.header += header;

    const int maxm = boys_max_m(data);

    arb_t t_mp;
    auto F_mp = std::unique_ptr<arb_t[]>(new arb_t[maxm+1]);
    arb_init(t_mp);
    for(int i = 0; i <= maxm; i++)
        arb_init(F_mp[i]);

    /* Target precision/accuracy, in bits, with a safety factor
       of 8 extra decimal digits */
    const slong target_prec = (ndigits+4) / MIRP_LOG_10_2;

    for(auto & it : data.values)
    {
        slong working_prec = target_prec;
        bool sufficient_accuracy = false;
        
        do {
            working_prec += 16;

            arb_set_str(t_mp, it.t.c_str(), working_prec);
            mirp_boys_interval(F_mp.get(), it.m, t_mp, working_prec); 

            if(arb_rel_accuracy_bits(F_mp[it.m]) >= target_prec)
                sufficient_accuracy = true;
                
        } while(!sufficient_accuracy);

        char * s = arb_get_str(F_mp[it.m], ndigits, ARB_STR_NO_RADIUS);
        it.value = s;
        free(s);
    }

    boys_write_file(outfile, data);
    arb_clear(t_mp);
    for(int i = 0; i <= maxm; i++)
        arb_clear(F_mp[i]);
}

    
} // closing namespace mirp

