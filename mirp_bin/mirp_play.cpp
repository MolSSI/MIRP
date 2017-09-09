/*! \file
 *
 * \brief mirp_run_test main function
 */

#include "mirp_bin/cmdline.hpp"
#include "mirp_bin/test_boys.hpp"
#include "mirp_bin/test_integral4.hpp"

#include <mirp/kernels/all.h>

#include <sstream>
#include <iostream>

using namespace mirp;


void test_overlap(int m, const char * t, slong prec1, slong prec2)
{
    arb_ptr F_1 = _arb_vec_init(m+1);
    arb_ptr F_2 = _arb_vec_init(m+1);


    mirp_boys_str(F_1, m, t, prec1);
    mirp_boys_str(F_2, m, t, prec2);

    char * s_1 = arb_get_str(F_1 + m, 75, ARB_STR_MORE); 
    char * s_2 = arb_get_str(F_2 + m, 75, ARB_STR_MORE); 


    printf("[%4ld] %s\n", prec1, s_1);
    printf("[%4ld] %s\n", prec2, s_2);
    const char * ovl = arb_overlaps(F_1 + m, F_2 + m) ? "Overlaps" : "DOES NOT OVERLAP";
    printf("%s\n", ovl);

    free(s_1);
    free(s_2);

    _arb_vec_clear(F_1, m+1);
    _arb_vec_clear(F_2, m+1);
}


int main(int argc, char ** argv)
{
    const char * t = "9e1";
    const int m = 2;

    const int max_prec = 256;
    arb_ptr F[max_prec];


    for(int prec = 1; prec <= max_prec; prec++)
    {
        F[prec-1] = _arb_vec_init(m+1);
        mirp_boys_str(F[prec-1], m, t, prec);
    }

    for(int i = 1; i < max_prec; i++)
    {
        int overlap = arb_overlaps(F[i]+m, F[i-1]+m);
        char * s = arb_get_str(F[i]+m, 75, ARB_STR_MORE);

        printf("[%4d] (%d) %s\n", i+1, overlap, s);

        free(s);
    }

    for(int prec = 1; prec <= max_prec; prec++)
        _arb_vec_clear(F[prec-1], m+1);

    test_overlap(m, t, 128, 256);

    return 0;
}
