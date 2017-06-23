#include "mirp/mirp.h"
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    /* Arbitrary precision as string */
    const char * x1_str = "0.1009671921";
    const char * y1_str = "0.0713958504";
    const char * z1_str = "0.0";
    const char * x2_str = "1.009671921";
    const char * y2_str = "0.713958504";
    const char * z2_str = "0.0";
    const char * x3_str = "2.009671921";
    const char * y3_str = "1.713958504";
    const char * z3_str = "0.0";
    const char * x4_str = "2.209671921";
    const char * y4_str = "1.813958504";
    const char * z4_str = "0.0";
    const char * xyz1_str[3] = { x1_str, y1_str, z1_str };
    const char * xyz2_str[3] = { x2_str, y2_str, z2_str };
    const char * xyz3_str[3] = { x3_str, y3_str, z3_str };
    const char * xyz4_str[3] = { x4_str, y4_str, z4_str };
    const char * alpha1_str = "13.0709311";
    const char * alpha2_str = "10.0709311";
    const char * alpha3_str = "1.30709311";
    const char * alpha4_str = "0.130709311";
    char * result_mp_str = NULL;

    /* Double precision */
    double xyz1[3] = { atof(x1_str), atof(y1_str), atof(z1_str) };
    double xyz2[3] = { atof(x2_str), atof(y2_str), atof(z2_str) };
    double xyz3[3] = { atof(x3_str), atof(y3_str), atof(z3_str) };
    double xyz4[3] = { atof(x4_str), atof(y4_str), atof(z4_str) };
    double alpha1 = atof(alpha1_str);
    double alpha2 = atof(alpha2_str);
    double alpha3 = atof(alpha3_str);
    double alpha4 = atof(alpha4_str);
    double result = 0.0;

    /* Arbitrary precision as mpfr_t */
    mpfr_t alpha1_mp, alpha2_mp, alpha3_mp, alpha4_mp;
    mpfr_t xyz1_mp[3], xyz2_mp[3], xyz3_mp[3], xyz4_mp[3];

    mpfr_init_set_str(alpha1_mp, alpha1_str, 10, MPFR_RNDN);
    mpfr_init_set_str(alpha2_mp, alpha2_str, 10, MPFR_RNDN);
    mpfr_init_set_str(alpha3_mp, alpha3_str, 10, MPFR_RNDN);
    mpfr_init_set_str(alpha4_mp, alpha4_str, 10, MPFR_RNDN);

    mpfr_init_set_str(xyz1_mp[0], x1_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz1_mp[1], y1_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz1_mp[2], z1_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz2_mp[0], x2_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz2_mp[1], y2_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz2_mp[2], z2_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz3_mp[0], x3_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz3_mp[1], y3_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz3_mp[2], z3_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz4_mp[0], x4_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz4_mp[1], y4_str, 10, MPFR_RNDN);
    mpfr_init_set_str(xyz4_mp[2], z4_str, 10, MPFR_RNDN);

    mpfr_t result_mp;
    mpfr_init2(result_mp, 256);

    /* Actually calculate */
    mirp_single_eri(&result,
                      1, 1, 1, alpha1, xyz1,
                      1, 1, 1, alpha2, xyz2,
                      1, 1, 1, alpha3, xyz3,
                      1, 1, 1, alpha4, xyz4);

    mirp_single_eri_mp(result_mp, 1, 1, 1, alpha1_mp, xyz1_mp,
                         1, 1, 1, alpha2_mp, xyz2_mp,
                         1, 1, 1, alpha3_mp, xyz3_mp,
                         1, 1, 1, alpha4_mp, xyz4_mp,
                         256);

    mirp_single_eri_mp_str(&result_mp_str,
                             1, 1, 1, alpha1_str, xyz1_str,
                             1, 1, 1, alpha2_str, xyz2_str,
                             1, 1, 1, alpha3_str, xyz3_str,
                             1, 1, 1, alpha4_str, xyz4_str,
                             256);

    printf("    result val: %-26.16e\n", result); 
    mpfr_printf(" mp result val: %Re\n", result_mp);
    printf("str result val: %s\n", result_mp_str);

    mpfr_free_str(result_mp_str);

    mpfr_clear(result_mp);
    mpfr_clears(alpha1_mp, alpha2_mp, alpha3_mp, alpha4_mp, (mpfr_ptr)0);
    mpfr_clears(xyz1_mp[0], xyz1_mp[1], xyz1_mp[2], (mpfr_ptr)0);
    mpfr_clears(xyz2_mp[0], xyz2_mp[1], xyz2_mp[2], (mpfr_ptr)0);
    mpfr_clears(xyz3_mp[0], xyz3_mp[1], xyz3_mp[2], (mpfr_ptr)0);
    mpfr_clears(xyz4_mp[0], xyz4_mp[1], xyz4_mp[2], (mpfr_ptr)0);

    return 0;
}
