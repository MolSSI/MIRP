/*! \file
 *
 * \brief A simple example of using MIRP as a library
 */


/* This example computes the ERI of water with HF/STO-3G
 * using arbitrary precision/interval arithmetic as well
 * as in double precision
 */

#include <mirp/mirp.h>
#include <stdio.h>

static void strarr_to_d(const char ** strarr, double * darr, int length)
{
    for(int i = 0; i < length; i++)
        darr[i] = atof(strarr[i]);
}

static void strarr_to_interval(const char ** strarr, arb_ptr iarr, int length, slong prec)
{
    for(int i = 0; i < length; i++)
        arb_set_str(iarr + i, strarr[i], prec);
}

static void print_integrals(const char * title,
                            const double * integrals_d,
                            arb_srcptr integrals_interval,
                            int length)
{
    printf("\n-- %s --\n", title);

    for(int i = 0; i < length; i++)
    {
        /* We will max out at 100 digits. Should be enough for this example */
        char * arbstr = arb_get_str(integrals_interval + i, 100, 0);
        printf("[%5d] %30.17e   %s\n", i, integrals_d[i], arbstr);
        free(arbstr);
    }

    printf("\n");
}

int main(void)
{
    /***************************************************
     * The infamous water STO-3G example
     ***************************************************/

    /* Holds only one shell quartet at a time. Only needs to
     * be 81, since we only have up to p shells
     */
    double integrals_d[81];
    arb_ptr integrals_interval = _arb_vec_init(81);

    /* Geometry, in bohr
     * Taken from CCCBDB (water HF/STO-3G) and converted to
     * bohr (keeping to the 4 decimal places shown in the CCCBDB)
     */
    const char * xyz_O[3] =  { "0.0",     "0.0",  "0.2402" };
    const char * xyz_H1[3] = { "0.0",  "1.4324", "-0.9609" };
    const char * xyz_H2[3] = { "0.0", "-1.4324", "-0.9609" };

    /******************************
     * Exponents and Coefficients
     ******************************/
    /* Note: MIRP supports general contractions, but not general contractions
     * with different angular momentum (ie, sp and spd shells). So we have to
     * split up the sp shell in STO-3G */
    const char * alpha_O_s1[3] = { "130.7093200", "23.8088610", "6.4436083" };
    const char * coeff_O_s1[3] = { "0.15432897",  "0.53532814", "0.44463454" };
    const char * alpha_O_s2[3] = { "5.0331513",   "1.1695961",  "0.3803890" };
    const char * coeff_O_s2[3] = { "-0.09996723", "0.39951283", "0.70011547" };
    const char * alpha_O_p[3] =  { "5.0331513",   "1.1695961",  "0.3803890" };
    const char * coeff_O_p[3] =  { "0.15591627",  "0.60768372", "0.39195739" };

    const char * alpha_H1_s[3] = { "3.42525091", "0.62391373", "0.16885540" };
    const char * coeff_H1_s[3] = { "0.15432897", "0.53532814", "0.44463454" };
    const char * alpha_H2_s[3] = { "3.42525091", "0.62391373", "0.16885540" };
    const char * coeff_H2_s[3] = { "0.15432897", "0.53532814", "0.44463454" };


    /**********************************************
     * Conversion of the above to double precision
     **********************************************/
    double xyz_O_d[3];
    double xyz_H1_d[3];
    double xyz_H2_d[3];
    strarr_to_d(xyz_O,  xyz_O_d,  3);
    strarr_to_d(xyz_H1, xyz_H2_d, 3);
    strarr_to_d(xyz_H2, xyz_H1_d, 3);

    double alpha_O_s1_d[3], coeff_O_s1_d[3];
    double alpha_O_s2_d[3], coeff_O_s2_d[3];
    double alpha_O_p_d[3],  coeff_O_p_d[3];
    strarr_to_d(alpha_O_s1, alpha_O_s1_d, 3);
    strarr_to_d(alpha_O_s2, alpha_O_s2_d, 3);
    strarr_to_d(alpha_O_p,  alpha_O_p_d,  3);
    strarr_to_d(coeff_O_s1, coeff_O_s1_d, 3);
    strarr_to_d(coeff_O_s2, coeff_O_s2_d, 3);
    strarr_to_d(coeff_O_p,  coeff_O_p_d,  3);

    double alpha_H1_s_d[3], coeff_H1_s_d[3];
    double alpha_H2_s_d[3], coeff_H2_s_d[3];
    strarr_to_d(alpha_H1_s, alpha_H1_s_d, 3);
    strarr_to_d(alpha_H2_s, alpha_H2_s_d, 3);
    strarr_to_d(coeff_H1_s, coeff_H1_s_d, 3);
    strarr_to_d(coeff_H2_s, coeff_H2_s_d, 3);


    /**********************************************
     * Conversion of the above to arb_t
     **********************************************/
    arb_ptr xyz_O_interval  = _arb_vec_init(3);
    arb_ptr xyz_H1_interval = _arb_vec_init(3);
    arb_ptr xyz_H2_interval = _arb_vec_init(3);
    strarr_to_interval(xyz_O,  xyz_O_interval,  3, 256);
    strarr_to_interval(xyz_H1, xyz_H2_interval, 3, 256);
    strarr_to_interval(xyz_H2, xyz_H1_interval, 3, 256);

    arb_ptr alpha_O_s1_interval = _arb_vec_init(3);
    arb_ptr alpha_O_s2_interval = _arb_vec_init(3);
    arb_ptr alpha_O_p_interval  = _arb_vec_init(3);
    arb_ptr coeff_O_s1_interval = _arb_vec_init(3);
    arb_ptr coeff_O_s2_interval = _arb_vec_init(3);
    arb_ptr coeff_O_p_interval  = _arb_vec_init(3);
    strarr_to_interval(alpha_O_s1, alpha_O_s1_interval, 3, 256);
    strarr_to_interval(alpha_O_s2, alpha_O_s2_interval, 3, 256);
    strarr_to_interval(alpha_O_p,  alpha_O_p_interval,  3, 256);
    strarr_to_interval(coeff_O_s1, coeff_O_s1_interval, 3, 256);
    strarr_to_interval(coeff_O_s2, coeff_O_s2_interval, 3, 256);
    strarr_to_interval(coeff_O_p,  coeff_O_p_interval,  3, 256);

    arb_ptr alpha_H1_s_interval = _arb_vec_init(3);
    arb_ptr alpha_H2_s_interval = _arb_vec_init(3);
    arb_ptr coeff_H1_s_interval = _arb_vec_init(3);
    arb_ptr coeff_H2_s_interval = _arb_vec_init(3);
    strarr_to_interval(alpha_H1_s, alpha_H1_s_interval, 3, 256);
    strarr_to_interval(alpha_H2_s, alpha_H2_s_interval, 3, 256);
    strarr_to_interval(coeff_H1_s, coeff_H1_s_interval, 3, 256);
    strarr_to_interval(coeff_H2_s, coeff_H2_s_interval, 3, 256);


    /*********************************************************************************************
     * Computing integrals
     * (not going to do them all, just a few examples)
     * NOTE: In my notation, { i, j, k, l } refer to the shell indices, NOT the angular momentum 
     *
     *  Ordering: 0 = Oxygen s1
     *            1 = Oxygen s2
     *            2 = Oxygen p
     *            3 = Hydrogen1 s
     *            4 = Hydrogen2 s
     *********************************************************************************************/


    /* {0, 0, 0, 0} */ 
    mirp_eri_d(integrals_d,
               0, xyz_O_d, 3, 1, alpha_O_s1_d, coeff_O_s1_d,
               0, xyz_O_d, 3, 1, alpha_O_s1_d, coeff_O_s1_d,
               0, xyz_O_d, 3, 1, alpha_O_s1_d, coeff_O_s1_d,
               0, xyz_O_d, 3, 1, alpha_O_s1_d, coeff_O_s1_d);

    mirp_eri(integrals_interval,
             0, xyz_O_interval, 3, 1, alpha_O_s1_interval, coeff_O_s1_interval,
             0, xyz_O_interval, 3, 1, alpha_O_s1_interval, coeff_O_s1_interval,
             0, xyz_O_interval, 3, 1, alpha_O_s1_interval, coeff_O_s1_interval,
             0, xyz_O_interval, 3, 1, alpha_O_s1_interval, coeff_O_s1_interval, 512);

    print_integrals("{0, 0, 0, 0}", integrals_d, integrals_interval, 1);


    /* {0, 1, 0, 1} */
    mirp_eri_d(integrals_d,
               0, xyz_O_d, 3, 1, alpha_O_s1_d, coeff_O_s1_d,
               0, xyz_O_d, 3, 1, alpha_O_s2_d, coeff_O_s2_d,
               0, xyz_O_d, 3, 1, alpha_O_s1_d, coeff_O_s1_d,
               0, xyz_O_d, 3, 1, alpha_O_s2_d, coeff_O_s2_d);

    mirp_eri(integrals_interval,
             0, xyz_O_interval, 3, 1, alpha_O_s1_interval, coeff_O_s1_interval,
             0, xyz_O_interval, 3, 1, alpha_O_s2_interval, coeff_O_s2_interval,
             0, xyz_O_interval, 3, 1, alpha_O_s1_interval, coeff_O_s1_interval,
             0, xyz_O_interval, 3, 1, alpha_O_s2_interval, coeff_O_s2_interval, 512);

    print_integrals("{0, 1, 0, 1}", integrals_d, integrals_interval, 1);


    /* {3, 2, 4, 1} */
    mirp_eri_d(integrals_d,
               0, xyz_H1_d, 3, 1, alpha_H1_s_d, coeff_H1_s_d,
               1, xyz_O_d,  3, 1, alpha_O_p_d,  coeff_O_p_d,
               0, xyz_H2_d, 3, 1, alpha_H2_s_d, coeff_H2_s_d,
               0, xyz_O_d,  3, 1, alpha_O_s2_d, coeff_O_s2_d);

    mirp_eri(integrals_interval,
             0, xyz_H1_interval, 3, 1, alpha_H1_s_interval, coeff_H1_s_interval,
             1, xyz_O_interval,  3, 1, alpha_O_p_interval,  coeff_O_p_interval,
             0, xyz_H2_interval, 3, 1, alpha_H2_s_interval, coeff_H2_s_interval,
             0, xyz_O_interval,  3, 1, alpha_O_s2_interval, coeff_O_s2_interval, 512);

    print_integrals("{3, 2, 4, 1}", integrals_d, integrals_interval, 3);


    /* {2, 2, 2, 2} (this is the (pp|pp) quartet */ 
    mirp_eri_d(integrals_d,
               1, xyz_O_d, 3, 1, alpha_O_p_d, coeff_O_p_d,
               1, xyz_O_d, 3, 1, alpha_O_p_d, coeff_O_p_d,
               1, xyz_O_d, 3, 1, alpha_O_p_d, coeff_O_p_d,
               1, xyz_O_d, 3, 1, alpha_O_p_d, coeff_O_p_d);

    mirp_eri(integrals_interval,
             1, xyz_O_interval, 3, 1, alpha_O_p_interval, coeff_O_p_interval,
             1, xyz_O_interval, 3, 1, alpha_O_p_interval, coeff_O_p_interval,
             1, xyz_O_interval, 3, 1, alpha_O_p_interval, coeff_O_p_interval,
             1, xyz_O_interval, 3, 1, alpha_O_p_interval, coeff_O_p_interval, 512);

    print_integrals("{2, 2, 2, 2}", integrals_d, integrals_interval, 81);


    /* Cleanup */
    _arb_vec_clear(integrals_interval, 81);
    _arb_vec_clear(xyz_O_interval, 3);
    _arb_vec_clear(xyz_H1_interval, 3);
    _arb_vec_clear(xyz_H2_interval, 3);
    _arb_vec_clear(alpha_O_s1_interval, 3);
    _arb_vec_clear(alpha_O_s2_interval, 3);
    _arb_vec_clear(alpha_O_p_interval,  3);
    _arb_vec_clear(coeff_O_s1_interval, 3);
    _arb_vec_clear(coeff_O_s2_interval, 3);
    _arb_vec_clear(coeff_O_p_interval,  3);
    _arb_vec_clear(alpha_H1_s_interval, 3);
    _arb_vec_clear(alpha_H2_s_interval, 3);
    _arb_vec_clear(coeff_H1_s_interval, 3);
    _arb_vec_clear(coeff_H2_s_interval, 3);

    return 0;
}

