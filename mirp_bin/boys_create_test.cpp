
#include "mirp_bin/boys_test_file.hpp"
#include <mirp/mirp.h>
#include <mirp/math.h>
#include <memory>
#include <cmath>
#include <iostream>

using namespace mirp;
    
namespace mirp {

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

} // close namespace mirp
