/*! \file
 *
 * \brief Some useful wrapping of functionality
 */

#include "mirp/wrappers.h"
#include "mirp/math.h"
#include "mirp/shell.h"

void mirp_integral4_single_exact(double * integral,
                                 const int * lmn1, const double * A, double alpha1,
                                 const int * lmn2, const double * B, double alpha2,
                                 const int * lmn3, const double * C, double alpha3,
                                 const int * lmn4, const double * D, double alpha4,
                                 cb_single4_interval cb)
{
    /* convert arguments to arb_t */
    arb_ptr A_mp = _arb_vec_init(3);
    arb_ptr B_mp = _arb_vec_init(3);
    arb_ptr C_mp = _arb_vec_init(3);
    arb_ptr D_mp = _arb_vec_init(3);

    arb_t alpha1_mp, alpha2_mp, alpha3_mp, alpha4_mp;
    arb_init(alpha1_mp);
    arb_init(alpha2_mp);
    arb_init(alpha3_mp);
    arb_init(alpha4_mp);

    arb_set_d(alpha1_mp, alpha1);
    arb_set_d(alpha2_mp, alpha2);
    arb_set_d(alpha3_mp, alpha3);
    arb_set_d(alpha4_mp, alpha4);

    for(int i = 0; i < 3; i++)
    {
        arb_set_d(A_mp + i, A[i]);
        arb_set_d(B_mp + i, B[i]);
        arb_set_d(C_mp + i, C[i]);
        arb_set_d(D_mp + i, D[i]);
    }

    /* Final integral output */
    arb_t integral_mp;
    arb_init(integral_mp);

    /* The target precision is the number of bits in double precision (53) + safety */
    const slong target_prec = 64;

    slong working_prec = target_prec;
    int sufficient_accuracy = 0;

    do {

        working_prec += 16;

        /* Call the callback */
        cb(integral_mp,
           lmn1, A_mp, alpha1_mp,
           lmn2, B_mp, alpha2_mp,
           lmn3, C_mp, alpha3_mp,
           lmn4, D_mp, alpha4_mp,
           working_prec);

        if(arb_rel_accuracy_bits(integral_mp) >= target_prec ||
           mirp_test_zero_prec(integral_mp, target_prec))
            sufficient_accuracy = 1;

    } while(!sufficient_accuracy);

    /* We get the value from the midpoint of the arb struct */
    *integral = arf_get_d(arb_midref(integral_mp), ARF_RND_NEAR);

    /* Cleanup */
    _arb_vec_clear(A_mp, 3);
    _arb_vec_clear(B_mp, 3);
    _arb_vec_clear(C_mp, 3);
    _arb_vec_clear(D_mp, 3);
    arb_clear(alpha1_mp);
    arb_clear(alpha2_mp);
    arb_clear(alpha3_mp);
    arb_clear(alpha4_mp);
    arb_clear(integral_mp);
}


void mirp_integral4_exact(double * integral,
                          int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1, 
                          int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2, 
                          int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3, 
                          int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4,
                          cb_integral4_interval cb)
{
    /* convert arguments to arb_t */
    arb_ptr A_mp = _arb_vec_init(3);
    arb_ptr B_mp = _arb_vec_init(3);
    arb_ptr C_mp = _arb_vec_init(3);
    arb_ptr D_mp = _arb_vec_init(3);

    arb_ptr alpha1_mp = _arb_vec_init(nprim1);
    arb_ptr alpha2_mp = _arb_vec_init(nprim2);
    arb_ptr alpha3_mp = _arb_vec_init(nprim3);
    arb_ptr alpha4_mp = _arb_vec_init(nprim4);
    arb_ptr coeff1_mp = _arb_vec_init(nprim1*ngeneral1);
    arb_ptr coeff2_mp = _arb_vec_init(nprim2*ngeneral2);
    arb_ptr coeff3_mp = _arb_vec_init(nprim3*ngeneral3);
    arb_ptr coeff4_mp = _arb_vec_init(nprim4*ngeneral4);

    for(int i = 0; i < 3; i++)
    {
        arb_set_d(A_mp + i, A[i]);
        arb_set_d(B_mp + i, B[i]);
        arb_set_d(C_mp + i, C[i]);
        arb_set_d(D_mp + i, D[i]);
    }

    for(int i = 0; i < nprim1; i++)
        arb_set_d(alpha1_mp + i, alpha1[i]);
    for(int i = 0; i < nprim2; i++)
        arb_set_d(alpha2_mp + i, alpha2[i]);
    for(int i = 0; i < nprim3; i++)
        arb_set_d(alpha3_mp + i, alpha3[i]);
    for(int i = 0; i < nprim4; i++)
        arb_set_d(alpha4_mp + i, alpha4[i]);

    for(int i = 0; i < nprim1 * ngeneral1; i++)
        arb_set_d(coeff1_mp + i, coeff1[i]);
    for(int i = 0; i < nprim2 * ngeneral2; i++)
        arb_set_d(coeff2_mp + i, coeff2[i]);
    for(int i = 0; i < nprim3 * ngeneral3; i++)
        arb_set_d(coeff3_mp + i, coeff3[i]);
    for(int i = 0; i < nprim4 * ngeneral4; i++)
        arb_set_d(coeff4_mp + i, coeff4[i]);


    /* Final integral output */
    const size_t ngen = ngeneral1 * ngeneral2 * ngeneral3 * ngeneral4;
    const size_t ncart = MIRP_NCART(am1) * MIRP_NCART(am2) * MIRP_NCART(am3) * MIRP_NCART(am4);
    const size_t nintegrals = ngen*ncart;
    arb_ptr integral_mp = _arb_vec_init(nintegrals);

    /* The target precision is the number of bits in double precision (53) + safety */
    const slong target_prec = 64;

    slong working_prec = target_prec;
    int sufficient_accuracy = 0;

    do {

        working_prec += 16;

        /* Call the callback */
        cb(integral_mp,
           am1, A_mp, nprim1, ngeneral1, alpha1_mp, coeff1_mp,
           am2, A_mp, nprim2, ngeneral2, alpha2_mp, coeff2_mp,
           am3, A_mp, nprim3, ngeneral3, alpha3_mp, coeff3_mp,
           am4, A_mp, nprim4, ngeneral4, alpha4_mp, coeff4_mp,
           working_prec);

        /* Check if all values have sufficient accuracy */
        sufficient_accuracy = 1;

        for(size_t j = 0; j < nintegrals; j++)
        {
            if(arb_rel_accuracy_bits(integral_mp + j) < target_prec &&
               !mirp_test_zero_prec(integral_mp, target_prec))
                sufficient_accuracy = 0;
        }

    } while(!sufficient_accuracy);

    /* We get the value from the midpoint of the arb struct */
    for(size_t i = 0; i < nintegrals; i++)
        integral[i] = arf_get_d(arb_midref(integral_mp + i), ARF_RND_NEAR);

    /* Cleanup */
    _arb_vec_clear(A_mp, 3);
    _arb_vec_clear(B_mp, 3);
    _arb_vec_clear(C_mp, 3);
    _arb_vec_clear(D_mp, 3);
    _arb_vec_clear(alpha1_mp, nprim1);
    _arb_vec_clear(alpha2_mp, nprim2);
    _arb_vec_clear(alpha3_mp, nprim3);
    _arb_vec_clear(alpha4_mp, nprim4);
    _arb_vec_clear(coeff1_mp, nprim1*ngeneral1);
    _arb_vec_clear(coeff2_mp, nprim2*ngeneral2);
    _arb_vec_clear(coeff3_mp, nprim3*ngeneral3);
    _arb_vec_clear(coeff4_mp, nprim4*ngeneral4);
    _arb_vec_clear(integral_mp, nintegrals);
}
