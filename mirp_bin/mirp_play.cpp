/*! \file
 *
 * \brief Sandbox for testing features
 */

#include "mirp_bin/cmdline.hpp"
#include <mirp/mirp.h>
#include <mirp/arb_help.h>
#include <mirp/math.h>
#include <iostream>

using namespace mirp;

int main(void)
{
    double result[1296];
    const double xyz_O[3] = { 0.0, 0.0, 0.0 };
    const double alpha_O_s[26] = { 54451891.10000000,
                                   13888721.40000000,
                                   4408873.420000000,
                                   1513672.400000000,
                                   553624.0620000000,
                                   211319.2440000000,
                                    83690.7483000000,
                                    34026.0317000000,
                                    14084.0491000000,
                                     5932.9989900000,
                                     2559.3576800000,
                                     1134.9182500000,
                                      515.3053460000,
                                      243.1742760000,
                                      118.6884990000,
                                       59.6159312000,
                                       30.2863196000,
                                       15.2365445000,
                                        7.6844570400,
                                        3.6410450600,
                                        1.7032524100,
                                        0.7325956700,
                                        0.2872888000,
                                        0.1149155200,
                                        0.0459662000,
                                        0.0183864800 };


    const double coeff_O_s[26] = { 0.00220712,
                                   0.00216389,
                                   0.00595756,
                                   0.00907607,
                                   0.01798328,
                                   0.03028123,
                                   0.05473004,
                                   0.09429869,
                                   0.16307628,
                                   0.25493454,
                                   0.31985104,
                                   0.23476028,
                                   0.05867659,
                                   0.00047976,
                                   0.00069094,
                                  -0.00043125,
                                   0.00022266,
                                  -0.00012497,
                                   0.00007405,
                                  -0.00003837,
                                   0.00001878,
                                  -0.00000860,
                                   0.00000419,
                                  -0.00000221,
                                   0.00000102,
                                  -0.00000029 };


    int ngeneral1 = 1;
    int ngeneral2 = 1;
    int ngeneral3 = 1;
    int ngeneral4 = 1;
    int nprim1 = 6;
    int nprim2 = 6;
    int nprim3 = 6;
    int nprim4 = 6;
    int am1 = 2;
    int am2 = 2;
    int am3 = 2;
    int am4 = 2;
    const size_t ncart1 = MIRP_NCART(am1);
    const size_t ncart2 = MIRP_NCART(am2);
    const size_t ncart3 = MIRP_NCART(am3);
    const size_t ncart4 = MIRP_NCART(am4);
    const size_t ncart1234 = ncart1*ncart2*ncart3*ncart4;
    size_t ngeneral1234 = ngeneral1*ngeneral2*ngeneral3*ngeneral4;

    mirp_eri_double(result,
                    am1, xyz_O, nprim1, ngeneral1, alpha_O_s, coeff_O_s, 
                    am2, xyz_O, nprim2, ngeneral2, alpha_O_s, coeff_O_s, 
                    am3, xyz_O, nprim3, ngeneral3, alpha_O_s, coeff_O_s, 
                    am4, xyz_O, nprim4, ngeneral4, alpha_O_s, coeff_O_s); 

    for(size_t n = 0; n < ngeneral1234; n++)
    for(size_t i = 0; i < ncart1234; i++)
        printf("DOUBLE: General %4lu Element %4lu : %30.16e\n", n, i, result[n*ncart1234+i]);


    const mpfr_prec_t final_prec = 128;
    const mpfr_prec_t working_prec = 128;

    arb_t result_interval[1296];
    arb_t xyz_O_interval[3], alpha_O_s_interval[nprim1], coeff_O_s_interval[nprim1*ngeneral1];
    mirp_init_arb_arr(result_interval, 1296);
    mirp_init_arb_arr(xyz_O_interval, 3);
    mirp_init_arb_arr(alpha_O_s_interval, nprim1);
    mirp_init_arb_arr(coeff_O_s_interval, nprim1*ngeneral1);

    for(int i = 0; i < 3; i++)
        arb_set_d(xyz_O_interval[i], xyz_O[i]);
    for(int i = 0; i < nprim1; i++)
        arb_set_d(alpha_O_s_interval[i], alpha_O_s[i]);
    for(int i = 0; i < nprim1*ngeneral1; i++)
        arb_set_d(coeff_O_s_interval[i], coeff_O_s[i]);

    mirp_eri_interval(result_interval,
                      am1, xyz_O_interval, nprim1, ngeneral1, alpha_O_s_interval, coeff_O_s_interval,
                      am2, xyz_O_interval, nprim2, ngeneral2, alpha_O_s_interval, coeff_O_s_interval,
                      am3, xyz_O_interval, nprim3, ngeneral3, alpha_O_s_interval, coeff_O_s_interval,
                      am4, xyz_O_interval, nprim4, ngeneral4, alpha_O_s_interval, coeff_O_s_interval,
                      working_prec);

    for(size_t n = 0; n < ngeneral1234; n++)
    for(size_t i = 0; i < ncart1234; i++)
    {
        char * s = arb_get_str(result_interval[n*ncart1234+i], 1000, 0);
        printf("INTERVAL: General %4lu Element %4lu : %s\n", n, i, s);
        free(s);
    }

    mirp_clear_arb_arr(result_interval, 1296);
    mirp_clear_arb_arr(xyz_O_interval, 3);
    mirp_clear_arb_arr(alpha_O_s_interval, nprim1);
    mirp_clear_arb_arr(coeff_O_s_interval, nprim1*ngeneral1);

    return 0;
}

