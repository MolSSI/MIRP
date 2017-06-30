#include "mirp_bin/cmdline.hpp"
#include <mirp/mirp.h>
#include <mirp/mpfr_help.h>
#include <mirp/arb_help.h>
#include <mirp/math.h>
#include <iostream>

using namespace mirp;

int main(void)
{
    double result[81];
    const double xyz_O[3] = { .1009671920798554, 0.07139585036771477, 0.0 };
    const double alpha_O_s[3] = { 130.7093200, 23.8088610, 6.4436083 };
    const double coeff_O_s[6] = { 4.2519432829437198, 4.1122937184311841, 1.2816225325813408, 
                                  4.2519432829437198, 4.1122937184311841, 1.2816225325813408 };

    int ngeneral1 = 2;
    int ngeneral2 = 2;
    int ngeneral3 = 2;
    int ngeneral4 = 2;
    int nprim1 = 3;
    int nprim2 = 3;
    int nprim3 = 3;
    int nprim4 = 3;
    int am1 = 0;
    int am2 = 0;
    int am3 = 0;
    int am4 = 0;
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
    const mpfr_prec_t working_prec = 256;
    mpfr_t result_mp[81];
    mpfr_t xyz_O_mp[3], alpha_O_s_mp[3], coeff_O_s_mp[6];
    mirp_init_mpfr_arr(result_mp, 81, final_prec);
    mirp_init_mpfr_arr(xyz_O_mp, 3, working_prec);
    mirp_init_mpfr_arr(alpha_O_s_mp, 3, working_prec);
    mirp_init_mpfr_arr(coeff_O_s_mp, 6, working_prec);

    for(int i = 0; i < 3; i++)
        mpfr_set_d(xyz_O_mp[i], xyz_O[i], MPFR_RNDN);
    for(int i = 0; i < nprim1; i++)
        mpfr_set_d(alpha_O_s_mp[i], alpha_O_s[i], MPFR_RNDN);
    for(int i = 0; i < nprim1*ngeneral1; i++)
        mpfr_set_d(coeff_O_s_mp[i], coeff_O_s[i], MPFR_RNDN);

    mirp_eri_mp(result_mp,
                am1, xyz_O_mp, nprim1, ngeneral1, alpha_O_s_mp, coeff_O_s_mp,
                am2, xyz_O_mp, nprim2, ngeneral2, alpha_O_s_mp, coeff_O_s_mp,
                am3, xyz_O_mp, nprim3, ngeneral3, alpha_O_s_mp, coeff_O_s_mp,
                am4, xyz_O_mp, nprim4, ngeneral4, alpha_O_s_mp, coeff_O_s_mp,
                working_prec);

    for(size_t n = 0; n < ngeneral1234; n++)
    for(size_t i = 0; i < ncart1234; i++)
        mpfr_printf("MP: General %4lu Element %4lu : %Re\n", n, i, result_mp[n*ncart1234+i]);

    mirp_clear_mpfr_arr(result_mp, 81);
    mirp_clear_mpfr_arr(xyz_O_mp, 3);
    mirp_clear_mpfr_arr(alpha_O_s_mp, 3);
    mirp_clear_mpfr_arr(coeff_O_s_mp, 6);

    arb_t result_interval[81];
    arb_t xyz_O_interval[3], alpha_O_s_interval[3], coeff_O_s_interval[6];
    mirp_init_arb_arr(result_interval, 81);
    mirp_init_arb_arr(xyz_O_interval, 3);
    mirp_init_arb_arr(alpha_O_s_interval, 3);
    mirp_init_arb_arr(coeff_O_s_interval, 6);

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

    mirp_clear_arb_arr(result_interval, 81);
    mirp_clear_arb_arr(xyz_O_interval, 3);
    mirp_clear_arb_arr(alpha_O_s_interval, 3);
    mirp_clear_arb_arr(coeff_O_s_interval, 6);

    return 0;
}

