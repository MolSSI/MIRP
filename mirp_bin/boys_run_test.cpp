
#include "mirp_bin/boys_test_file.hpp"
#include <mirp/mirp.h>
#include <memory>
#include <cmath>
#include <iostream>

namespace mirp {

static long boys_run_test_interval(const boys_data & data, long extra_m, long target_prec, long working_prec)
{
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
        arb_set_str(t_mp, it.t.c_str(), working_prec);
        arb_set_str(vref_mp, it.value.c_str(), target_prec);
        mirp_boys_interval(F_mp.get(), it.m+extra_m, t_mp, working_prec);

        arb_set_round(F_mp[it.m], F_mp[it.m], target_prec);

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

static long boys_run_test_double(const boys_data & data, long extra_m)
{
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


long boys_run_test(const std::string & filename, const std::string & floattype, long extra_m, long target_prec, long working_prec)
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
        nfailed = boys_run_test_interval(data, extra_m, target_prec, working_prec);
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

} // close namespace mirp
