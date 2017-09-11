/*! \file
 *
 * \brief Some useful wrapping of functionality
 */

#include "mirp/pragma.h"
#include "mirp/math.h"
#include "mirp/shell.h"
#include "mirp/kernels/boys.h"
#include "mirp/kernels/integral4_wrappers.h"
#include <string.h> /* for memset */
#include <assert.h>



/*! \brief Compute all cartesian components of a single primitive integral
 *         (interval arithmetic)
 *
 * \copydetails mirp_cartloop4_d
 * \param [in] working_prec The working precision (binary digits/bits) to use in the calculation
 */
static void mirp_cartloop4(arb_ptr integral,
                           int am1, arb_srcptr A, const arb_t alpha1,
                           int am2, arb_srcptr B, const arb_t alpha2,
                           int am3, arb_srcptr C, const arb_t alpha3,
                           int am4, arb_srcptr D, const arb_t alpha4,
                           slong working_prec, cb_integral4_single cb)
{
    assert(am1 >= 0);
    assert(am2 >= 0);
    assert(am3 >= 0);
    assert(am4 >= 0);

    const long ncart1 = MIRP_NCART(am1);
    const long ncart2 = MIRP_NCART(am2);
    const long ncart3 = MIRP_NCART(am3);
    const long ncart4 = MIRP_NCART(am4);

    int lmn1[ncart1][3];
    int lmn2[ncart2][3];
    int lmn3[ncart3][3];
    int lmn4[ncart4][3];

    mirp_gaussian_fill_lmn(am1, (int*)lmn1);
    mirp_gaussian_fill_lmn(am2, (int*)lmn2);
    mirp_gaussian_fill_lmn(am3, (int*)lmn3);
    mirp_gaussian_fill_lmn(am4, (int*)lmn4);


    #ifdef _OPENMP
    #pragma omp parallel for collapse(4)
    #endif
    for(long i = 0; i < ncart1; i++)
    for(long j = 0; j < ncart2; j++)
    for(long k = 0; k < ncart3; k++)
    for(long l = 0; l < ncart4; l++)
    {
        const long idx = i*ncart4*ncart3*ncart2
                       + j*ncart4*ncart3
                       + k*ncart4
                       + l;

        cb(integral + idx,
           lmn1[i], A, alpha1,
           lmn2[j], B, alpha2,
           lmn3[k], C, alpha3,
           lmn4[l], D, alpha4,
           working_prec);
    }
}


/*! \brief Compute all cartesian components of a single primitive integral
 *         (four-center, double precision)
 *
 * The \p integral buffer is expected to be able to hold all primitive integrals
 * (ie, it can hold ncart(am1) * ncart(am2) * ncart(am3) * ncart(am4) elements).
 *
 * \param [out] integral
 *              Resulting integral integral
 * \param [in]  am1,am2,am3,am4
 *              Angular momentum of the four-centers
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four-centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four-centers
 * \param [in]  cb
 *              Callback that calculates a single cartesian component of a
 *              primitive integral
 */
static void mirp_cartloop4_d(double * integral,
                             int am1, const double * A, double alpha1,
                             int am2, const double * B, double alpha2,
                             int am3, const double * C, double alpha3,
                             int am4, const double * D, double alpha4,
                             cb_integral4_single_d cb)
{
    assert(am1 >= 0);
    assert(am2 >= 0);
    assert(am3 >= 0);
    assert(am4 >= 0);

    const long ncart1 = MIRP_NCART(am1);
    const long ncart2 = MIRP_NCART(am2);
    const long ncart3 = MIRP_NCART(am3);
    const long ncart4 = MIRP_NCART(am4);

    int lmn1[ncart1][3];
    int lmn2[ncart2][3];
    int lmn3[ncart3][3];
    int lmn4[ncart4][3];

    mirp_gaussian_fill_lmn(am1, (int*)lmn1);
    mirp_gaussian_fill_lmn(am2, (int*)lmn2);
    mirp_gaussian_fill_lmn(am3, (int*)lmn3);
    mirp_gaussian_fill_lmn(am4, (int*)lmn4);


    #ifdef _OPENMP
    #pragma omp parallel for collapse(4)
    #endif
    for(long i = 0; i < ncart1; i++)
    for(long j = 0; j < ncart2; j++)
    for(long k = 0; k < ncart3; k++)
    for(long l = 0; l < ncart4; l++)
    {
        const long idx = i*ncart4*ncart3*ncart2
                       + j*ncart4*ncart3
                       + k*ncart4
                       + l;

        cb(integral + idx,
           lmn1[i], A, alpha1,
           lmn2[j], B, alpha2,
           lmn3[k], C, alpha3,
           lmn4[l], D, alpha4);
    }
}


void mirp_loop_shell4(arb_ptr integral,
                      int am1, arb_srcptr A, int nprim1, int ngen1, arb_srcptr alpha1, arb_srcptr coeff1,
                      int am2, arb_srcptr B, int nprim2, int ngen2, arb_srcptr alpha2, arb_srcptr coeff2,
                      int am3, arb_srcptr C, int nprim3, int ngen3, arb_srcptr alpha3, arb_srcptr coeff3,
                      int am4, arb_srcptr D, int nprim4, int ngen4, arb_srcptr alpha4, arb_srcptr coeff4,
                      slong working_prec, cb_integral4_single cb)
{
    assert(am1 >= 0); assert(nprim1 > 0); assert(ngen1 > 0);
    assert(am2 >= 0); assert(nprim2 > 0); assert(ngen2 > 0);
    assert(am3 >= 0); assert(nprim3 > 0); assert(ngen3 > 0);
    assert(am4 >= 0); assert(nprim4 > 0); assert(ngen4 > 0);

    const long ncart1234 = MIRP_NCART4(am1, am2, am3, am4);
    const long ngen1234 = ngen1*ngen2*ngen3*ngen4;
    const long full_size = ncart1234*ngen1234;

    arb_ptr integral_buffer = _arb_vec_init(ncart1234);
    arb_ptr coeff1_norm = _arb_vec_init(nprim1 * ngen1);
    arb_ptr coeff2_norm = _arb_vec_init(nprim2 * ngen2);
    arb_ptr coeff3_norm = _arb_vec_init(nprim3 * ngen3);
    arb_ptr coeff4_norm = _arb_vec_init(nprim4 * ngen4);

    mirp_normalize_shell(am1, nprim1, ngen1, alpha1, coeff1, coeff1_norm, working_prec);
    mirp_normalize_shell(am2, nprim2, ngen2, alpha2, coeff2, coeff2_norm, working_prec);
    mirp_normalize_shell(am3, nprim3, ngen3, alpha3, coeff3, coeff3_norm, working_prec);
    mirp_normalize_shell(am4, nprim4, ngen4, alpha4, coeff4, coeff4_norm, working_prec);

    _arb_vec_zero(integral, full_size);


    for(int i = 0; i < nprim1; i++)
    for(int j = 0; j < nprim2; j++)
    for(int k = 0; k < nprim3; k++)
    for(int l = 0; l < nprim4; l++)
    {
        mirp_cartloop4(integral_buffer,
                       am1, A, alpha1 + i,
                       am2, B, alpha2 + j,
                       am3, C, alpha3 + k,
                       am4, D, alpha4 + l,
                       working_prec, cb);

        #ifdef _OPENMP
        #pragma omp parallel for collapse(4)
        #endif
        for(int m = 0; m < ngen1; m++)
        for(int n = 0; n < ngen2; n++)
        for(int o = 0; o < ngen3; o++)
        for(int p = 0; p < ngen4; p++)
        {
            /* A temporary variable (used to build up the coefficient) */
            arb_t coeff;
            arb_init(coeff);

            arb_mul(coeff, coeff1_norm+(m*nprim1+i), coeff2_norm+(n*nprim2+j), working_prec);
            arb_mul(coeff, coeff,                    coeff3_norm+(o*nprim3+k), working_prec);
            arb_mul(coeff, coeff,                    coeff4_norm+(p*nprim4+l), working_prec);

            const long start = ncart1234*(
                               m*ngen4*ngen3*ngen2
                             + n*ngen4*ngen3
                             + o*ngen4
                             + p);

            for(long q = 0; q < ncart1234; q++)
                arb_addmul(integral+start+q, integral_buffer+q, coeff, working_prec);

            arb_clear(coeff);
        }
    }

    _arb_vec_clear(integral_buffer, ncart1234);
    _arb_vec_clear(coeff1_norm, nprim1*ngen1);
    _arb_vec_clear(coeff2_norm, nprim2*ngen2);
    _arb_vec_clear(coeff3_norm, nprim3*ngen3);
    _arb_vec_clear(coeff4_norm, nprim4*ngen4);
}

void mirp_integral4_single_str(arb_t integral,
                               const int * lmn1, const char ** A, const char * alpha1,
                               const int * lmn2, const char ** B, const char * alpha2,
                               const int * lmn3, const char ** C, const char * alpha3,
                               const int * lmn4, const char ** D, const char * alpha4,
                               slong working_prec, cb_integral4_single cb)
{
    assert(lmn1[0] > 0); assert(lmn1[1] >= 0); assert(lmn1[2] >= 0);
    assert(lmn2[0] > 0); assert(lmn2[1] >= 0); assert(lmn2[2] >= 0);
    assert(lmn3[0] > 0); assert(lmn3[1] >= 0); assert(lmn3[2] >= 0);
    assert(lmn4[0] > 0); assert(lmn4[1] >= 0); assert(lmn4[2] >= 0);

    /* Similar to mirp_integral4_single_target, but should
     * always succeed */
    arb_ptr A_mp = _arb_vec_init(3);
    arb_ptr B_mp = _arb_vec_init(3);
    arb_ptr C_mp = _arb_vec_init(3);
    arb_ptr D_mp = _arb_vec_init(3);
    arb_t alpha1_mp, alpha2_mp, alpha3_mp, alpha4_mp;
    arb_init(alpha1_mp);
    arb_init(alpha2_mp);
    arb_init(alpha3_mp);
    arb_init(alpha4_mp);

    for(int i = 0; i < 3; i++)
    {
        arb_set_str(A_mp + i, A[i], working_prec);
        arb_set_str(B_mp + i, B[i], working_prec);
        arb_set_str(C_mp + i, C[i], working_prec);
        arb_set_str(D_mp + i, D[i], working_prec);
    }

    arb_set_str(alpha1_mp, alpha1, working_prec);
    arb_set_str(alpha2_mp, alpha2, working_prec);
    arb_set_str(alpha3_mp, alpha3, working_prec);
    arb_set_str(alpha4_mp, alpha4, working_prec);

    cb(integral,
       lmn1, A_mp, alpha1_mp,
       lmn2, B_mp, alpha2_mp,
       lmn3, C_mp, alpha3_mp,
       lmn4, D_mp, alpha4_mp,
       working_prec);


    _arb_vec_clear(A_mp, 3);
    _arb_vec_clear(B_mp, 3);
    _arb_vec_clear(C_mp, 3);
    _arb_vec_clear(D_mp, 3);
    arb_clear(alpha1_mp);
    arb_clear(alpha2_mp);
    arb_clear(alpha3_mp);
    arb_clear(alpha4_mp);
}


void mirp_integral4_str(arb_ptr integrals,
                        int am1, const char ** A, int nprim1, int ngen1, const char ** alpha1, const char ** coeff1,
                        int am2, const char ** B, int nprim2, int ngen2, const char ** alpha2, const char ** coeff2,
                        int am3, const char ** C, int nprim3, int ngen3, const char ** alpha3, const char ** coeff3,
                        int am4, const char ** D, int nprim4, int ngen4, const char ** alpha4, const char ** coeff4,
                        slong working_prec, cb_integral4 cb)
{
    assert(am1 >= 0); assert(nprim1 > 0); assert(ngen1 > 0);
    assert(am2 >= 0); assert(nprim2 > 0); assert(ngen2 > 0);
    assert(am3 >= 0); assert(nprim3 > 0); assert(ngen3 > 0);
    assert(am4 >= 0); assert(nprim4 > 0); assert(ngen4 > 0);

    arb_ptr A_mp = _arb_vec_init(3);
    arb_ptr B_mp = _arb_vec_init(3);
    arb_ptr C_mp = _arb_vec_init(3);
    arb_ptr D_mp = _arb_vec_init(3);
    arb_ptr alpha1_mp = _arb_vec_init(nprim1);
    arb_ptr alpha2_mp = _arb_vec_init(nprim2);
    arb_ptr alpha3_mp = _arb_vec_init(nprim3);
    arb_ptr alpha4_mp = _arb_vec_init(nprim4);
    arb_ptr coeff1_mp = _arb_vec_init(nprim1*ngen1);
    arb_ptr coeff2_mp = _arb_vec_init(nprim2*ngen2);
    arb_ptr coeff3_mp = _arb_vec_init(nprim3*ngen3);
    arb_ptr coeff4_mp = _arb_vec_init(nprim4*ngen4);

    for(int i = 0; i < 3; i++)
    {
        arb_set_str(A_mp + i, A[i], working_prec);
        arb_set_str(B_mp + i, B[i], working_prec);
        arb_set_str(C_mp + i, C[i], working_prec);
        arb_set_str(D_mp + i, D[i], working_prec);
    }

    for(int i = 0; i < nprim1; i++)
        arb_set_str(alpha1_mp + i, alpha1[i], working_prec);

    for(int i = 0; i < nprim2; i++)
        arb_set_str(alpha2_mp + i, alpha2[i], working_prec);

    for(int i = 0; i < nprim3; i++)
        arb_set_str(alpha3_mp + i, alpha3[i], working_prec);

    for(int i = 0; i < nprim4; i++)
        arb_set_str(alpha4_mp + i, alpha4[i], working_prec);

    for(int i = 0; i < nprim1*ngen1; i++)
        arb_set_str(coeff1_mp + i, coeff1[i], working_prec);

    for(int i = 0; i < nprim2*ngen2; i++)
        arb_set_str(coeff2_mp + i, coeff2[i], working_prec);

    for(int i = 0; i < nprim3*ngen3; i++)
        arb_set_str(coeff3_mp + i, coeff3[i], working_prec);

    for(int i = 0; i < nprim4*ngen4; i++)
        arb_set_str(coeff4_mp + i, coeff4[i], working_prec);


    cb(integrals,
       am1, A_mp, nprim1, ngen1, alpha1_mp, coeff1_mp,
       am2, B_mp, nprim2, ngen2, alpha2_mp, coeff2_mp,
       am3, C_mp, nprim3, ngen3, alpha3_mp, coeff3_mp,
       am4, D_mp, nprim4, ngen4, alpha4_mp, coeff4_mp,
       working_prec);


    _arb_vec_clear(A_mp, 3);
    _arb_vec_clear(B_mp, 3);
    _arb_vec_clear(C_mp, 3);
    _arb_vec_clear(D_mp, 3);
    _arb_vec_clear(alpha1_mp, nprim1);
    _arb_vec_clear(alpha2_mp, nprim2);
    _arb_vec_clear(alpha3_mp, nprim3);
    _arb_vec_clear(alpha4_mp, nprim4);
    _arb_vec_clear(coeff1_mp, nprim1*ngen1);
    _arb_vec_clear(coeff2_mp, nprim2*ngen2);
    _arb_vec_clear(coeff3_mp, nprim3*ngen3);
    _arb_vec_clear(coeff4_mp, nprim4*ngen4);
}



void mirp_integral4_single_exact(double * integral,
                                 const int * lmn1, const double * A, double alpha1,
                                 const int * lmn2, const double * B, double alpha2,
                                 const int * lmn3, const double * C, double alpha3,
                                 const int * lmn4, const double * D, double alpha4,
                                 cb_integral4_single cb)
{
    assert(lmn1[0] > 0); assert(lmn1[1] >= 0); assert(lmn1[2] >= 0);
    assert(lmn2[0] > 0); assert(lmn2[1] >= 0); assert(lmn2[2] >= 0);
    assert(lmn3[0] > 0); assert(lmn3[1] >= 0); assert(lmn3[2] >= 0);
    assert(lmn4[0] > 0); assert(lmn4[1] >= 0); assert(lmn4[2] >= 0);

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

    /* The target precision is the number of bits in
     * double precision (53) + lots of safety */
    const slong target_prec = 64;

    slong working_prec = target_prec;
    int suff_acc = 0;

    /* for comparisons */
    arf_t ubound, lbound;
    arf_init(ubound);
    arf_init(lbound);

    while(!suff_acc)
    {
        working_prec += target_prec;

        /* Call the callback */
        cb(integral_mp,
           lmn1, A_mp, alpha1_mp,
           lmn2, B_mp, alpha2_mp,
           lmn3, C_mp, alpha3_mp,
           lmn4, D_mp, alpha4_mp,
           working_prec);

        /* Do we have sufficient accuracy? We need at least
         * 53 bits + 11 bits safety OR the value is zero (has zero precision)
         * and the error bounds is exactly zero when converted to double precision */
        slong bits = arb_rel_accuracy_bits(integral_mp);

        suff_acc = 1;
        if(bits > 0 && bits < target_prec)
            suff_acc = 0;
        else if(bits <= 0)
        {
            arb_get_ubound_arf(ubound, integral_mp, working_prec);
            arb_get_lbound_arf(lbound, integral_mp, working_prec);


            PRAGMA_WARNING_PUSH
            PRAGMA_WARNING_IGNORE_FP_UNDERFLOW

            if(arf_cmpabs_d(lbound, MIRP_DBL_TRUE_MIN) > 0 || 
                arf_cmpabs_d(ubound, MIRP_DBL_TRUE_MIN) > 0)
                suff_acc = 0; 

            PRAGMA_WARNING_POP
        }
    }

    /* We get the value from the midpoint of the arb struct */
    *integral = arf_get_d(arb_midref(integral_mp), ARF_RND_NEAR);

    /* Cleanup */
    arf_clear(lbound);
    arf_clear(ubound);
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
                          int am1, const double * A, int nprim1, int ngen1, const double * alpha1, const double * coeff1,
                          int am2, const double * B, int nprim2, int ngen2, const double * alpha2, const double * coeff2,
                          int am3, const double * C, int nprim3, int ngen3, const double * alpha3, const double * coeff3,
                          int am4, const double * D, int nprim4, int ngen4, const double * alpha4, const double * coeff4,
                          cb_integral4 cb)
{
    assert(am1 >= 0); assert(nprim1 > 0); assert(ngen1 > 0);
    assert(am2 >= 0); assert(nprim2 > 0); assert(ngen2 > 0);
    assert(am3 >= 0); assert(nprim3 > 0); assert(ngen3 > 0);
    assert(am4 >= 0); assert(nprim4 > 0); assert(ngen4 > 0);

    /* convert arguments to arb_t */
    arb_ptr A_mp = _arb_vec_init(3);
    arb_ptr B_mp = _arb_vec_init(3);
    arb_ptr C_mp = _arb_vec_init(3);
    arb_ptr D_mp = _arb_vec_init(3);

    arb_ptr alpha1_mp = _arb_vec_init(nprim1);
    arb_ptr alpha2_mp = _arb_vec_init(nprim2);
    arb_ptr alpha3_mp = _arb_vec_init(nprim3);
    arb_ptr alpha4_mp = _arb_vec_init(nprim4);
    arb_ptr coeff1_mp = _arb_vec_init(nprim1*ngen1);
    arb_ptr coeff2_mp = _arb_vec_init(nprim2*ngen2);
    arb_ptr coeff3_mp = _arb_vec_init(nprim3*ngen3);
    arb_ptr coeff4_mp = _arb_vec_init(nprim4*ngen4);

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

    for(int i = 0; i < nprim1*ngen1; i++)
        arb_set_d(coeff1_mp + i, coeff1[i]);
    for(int i = 0; i < nprim2*ngen2; i++)
        arb_set_d(coeff2_mp + i, coeff2[i]);
    for(int i = 0; i < nprim3*ngen3; i++)
        arb_set_d(coeff3_mp + i, coeff3[i]);
    for(int i = 0; i < nprim4*ngen4; i++)
        arb_set_d(coeff4_mp + i, coeff4[i]);


    /* Final integral output */
    const long ngen = ngen1 * ngen2 * ngen3 * ngen4;
    const long ncart = MIRP_NCART4(am1, am2, am3, am4);
    const long nintegrals = ngen*ncart;
    arb_ptr integral_mp = _arb_vec_init(nintegrals);

    /* The target precision is the number of bits in double precision (53) + safety */
    const slong target_prec = 64;

    slong working_prec = target_prec;
    int suff_acc = 0;

    /* for comparisons */
    arf_t ubound, lbound;
    arf_init(ubound);
    arf_init(lbound);

    while(!suff_acc)
    {
        working_prec += target_prec;

        /* Call the callback */
        cb(integral_mp,
           am1, A_mp, nprim1, ngen1, alpha1_mp, coeff1_mp,
           am2, B_mp, nprim2, ngen2, alpha2_mp, coeff2_mp,
           am3, C_mp, nprim3, ngen3, alpha3_mp, coeff3_mp,
           am4, D_mp, nprim4, ngen4, alpha4_mp, coeff4_mp,
           working_prec);

        suff_acc = 1;
        for(long i = 0; i < nintegrals; i++)
        {
            /* Do we have sufficient accuracy? We need at least
             * 53 bits + 11 bits safety OR the value is zero (has zero precision)
             * and the error bounds is exactly zero when converted to double precision */
            slong bits = arb_rel_accuracy_bits(integral_mp + i);

            if(bits > 0 && bits < target_prec)
                suff_acc = 0;
            else if(bits <= 0)
            {
                arb_get_ubound_arf(ubound, integral_mp + i, working_prec);
                arb_get_lbound_arf(lbound, integral_mp + i, working_prec);

                PRAGMA_WARNING_PUSH
                PRAGMA_WARNING_IGNORE_FP_UNDERFLOW

                if(arf_cmpabs_d(lbound, MIRP_DBL_TRUE_MIN) > 0 || 
                   arf_cmpabs_d(ubound, MIRP_DBL_TRUE_MIN) > 0)
                    suff_acc = 0; 

                PRAGMA_WARNING_POP
            }
        }
    }


    /* We get the value from the midpoint of the arb struct */
    for(long i = 0; i < nintegrals; i++)
    {
        if(arb_rel_accuracy_bits(integral_mp + i) <= 0)
            integral[i] = 0.0;
        else
            integral[i] = arf_get_d(arb_midref(integral_mp + i), ARF_RND_NEAR);
    }

    /* Cleanup */
    arf_clear(lbound);
    arf_clear(ubound);
    _arb_vec_clear(A_mp, 3);
    _arb_vec_clear(B_mp, 3);
    _arb_vec_clear(C_mp, 3);
    _arb_vec_clear(D_mp, 3);
    _arb_vec_clear(alpha1_mp, nprim1);
    _arb_vec_clear(alpha2_mp, nprim2);
    _arb_vec_clear(alpha3_mp, nprim3);
    _arb_vec_clear(alpha4_mp, nprim4);
    _arb_vec_clear(coeff1_mp, nprim1*ngen1);
    _arb_vec_clear(coeff2_mp, nprim2*ngen2);
    _arb_vec_clear(coeff3_mp, nprim3*ngen3);
    _arb_vec_clear(coeff4_mp, nprim4*ngen4);
    _arb_vec_clear(integral_mp, nintegrals);
}


void mirp_loop_shell4_d(double * integral,
                        int am1, const double * A, int nprim1, int ngen1, const double * alpha1, const double * coeff1,
                        int am2, const double * B, int nprim2, int ngen2, const double * alpha2, const double * coeff2,
                        int am3, const double * C, int nprim3, int ngen3, const double * alpha3, const double * coeff3,
                        int am4, const double * D, int nprim4, int ngen4, const double * alpha4, const double * coeff4,
                        cb_integral4_single_d cb)
{
    assert(am1 >= 0); assert(nprim1 > 0); assert(ngen1 > 0);
    assert(am2 >= 0); assert(nprim2 > 0); assert(ngen2 > 0);
    assert(am3 >= 0); assert(nprim3 > 0); assert(ngen3 > 0);
    assert(am4 >= 0); assert(nprim4 > 0); assert(ngen4 > 0);

    const long ncart1 = MIRP_NCART(am1);
    const long ncart2 = MIRP_NCART(am2);
    const long ncart3 = MIRP_NCART(am3);
    const long ncart4 = MIRP_NCART(am4);
    const long ncart1234 = ncart1*ncart2*ncart3*ncart4;
    const long ngen1234 = ngen1*ngen2*ngen3*ngen4;
    const size_t full_size = (size_t)(ncart1234*ngen1234);

    double * integral_buffer = malloc(full_size * sizeof(double));
    double * coeff1_norm = malloc((size_t)(nprim1*ngen1) * sizeof(double));
    double * coeff2_norm = malloc((size_t)(nprim2*ngen2) * sizeof(double));
    double * coeff3_norm = malloc((size_t)(nprim3*ngen3) * sizeof(double));
    double * coeff4_norm = malloc((size_t)(nprim4*ngen4) * sizeof(double));

    mirp_normalize_shell_d(am1, nprim1, ngen1, alpha1, coeff1, coeff1_norm);
    mirp_normalize_shell_d(am2, nprim2, ngen2, alpha2, coeff2, coeff2_norm);
    mirp_normalize_shell_d(am3, nprim3, ngen3, alpha3, coeff3, coeff3_norm);
    mirp_normalize_shell_d(am4, nprim4, ngen4, alpha4, coeff4, coeff4_norm);

    memset(integral, 0, full_size * sizeof(double));
    for(int i = 0; i < nprim1; i++)
    for(int j = 0; j < nprim2; j++)
    for(int k = 0; k < nprim3; k++)
    for(int l = 0; l < nprim4; l++)
    {
        mirp_cartloop4_d(integral_buffer,
                              am1, A, alpha1[i],
                              am2, B, alpha2[j],
                              am3, C, alpha3[k],
                              am4, D, alpha4[l],
                              cb);

        #ifdef _OPENMP
        #pragma omp parallel for collapse(4)
        #endif
        for(int m = 0; m < ngen1; m++)
        for(int n = 0; n < ngen2; n++)
        for(int o = 0; o < ngen3; o++)
        for(int p = 0; p < ngen4; p++)
        {
            const double coeff = coeff1_norm[m*nprim1+i]
                               * coeff2_norm[n*nprim2+j]
                               * coeff3_norm[o*nprim3+k]
                               * coeff4_norm[p*nprim4+l];

            const long start = ncart1234*(
                               m*ngen4*ngen3*ngen2
                             + n*ngen4*ngen3
                             + o*ngen4
                             + p);

            for(long q = 0; q < ncart1234; q++)
                integral[start*ncart1234+q] += integral_buffer[q] * coeff;
        }
    }

    free(integral_buffer);
    free(coeff1_norm);
    free(coeff2_norm);
    free(coeff3_norm);
    free(coeff4_norm);
}

