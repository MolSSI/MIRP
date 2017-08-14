/*! \file
 *
 * \brief Some useful wrapping of functionality
 */

#pragma once

#include "mirp/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Compute all cartesian components of a contracted shell quartet
 *         (interval arithmetic)
 *
 * \copydetails mirp_loop4_d
 * \param [in] working_prec The working precision (binary digits/bits) to use in the calculation
 */
void mirp_loop_shell4(arb_ptr output,
                      int am1, arb_srcptr A, int nprim1, int ngeneral1, arb_srcptr alpha1, arb_srcptr coeff1,
                      int am2, arb_srcptr B, int nprim2, int ngeneral2, arb_srcptr alpha2, arb_srcptr coeff2,
                      int am3, arb_srcptr C, int nprim3, int ngeneral3, arb_srcptr alpha3, arb_srcptr coeff3,
                      int am4, arb_srcptr D, int nprim4, int ngeneral4, arb_srcptr alpha4, arb_srcptr coeff4,
                      slong working_prec, cb_integral4_single cb);



void mirp_loop_shell4_d(double * integral,
                        int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1,
                        int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2,
                        int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3,
                        int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4,
                        cb_integral4_single_d cb);


/*! \brief Compute a single 4-center integral to a target precision
 *
 * This function seeks to calculate the integral to a given precision. This works
 * by increasing the working precision until the result is the desired
 * precision (\p target_prec).
 *
 * If the inputs are not of high enough precision to reach the target
 * precision, the function will return nonzero, indicating that reaching the target
 * precision is impossible. Otherwise, the function returns zero.
 *
 * Note that the result may be of higher precision that the
 * desired target precision.
 *
 * \param [out] integral
 *              Resulting integral output
 * \param [in]  lmn1,lmn2,lmn3,lmn4
 *              Exponents of x, y, and z that signify angular momentum. Required
 *              to be 3 elements.
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four centers
 * \param [in]  cb
 *              The callback function that computes integrals in interval arithmetic
 */
int mirp_integral4_single_target(arb_t integral,
                                      const int * lmn1, arb_srcptr A, const arb_t alpha1,
                                      const int * lmn2, arb_srcptr B, const arb_t alpha2,
                                      const int * lmn3, arb_srcptr C, const arb_t alpha3,
                                      const int * lmn4, arb_srcptr D, const arb_t alpha4,
                                      slong target_prec, cb_integral4_single cb);


/*! \brief Compute a single 4-center integral to a target precision (string input)
 *
 * This function seeks to calculate the integral to a given
 * precision. This works by increasing the working precision
 * until the result is the desired precision (\p target_prec).
 *
 * The input is interpreted to be to infinite precision (ie,
 * there are an infinite number of trailing zeroes after the
 * decimal point).
 *
 * Note that the result may be of higher precision that the
 * desired target precision.
 *
 *
 * \param [out] integral
 *              Resulting integral output
 * \param [in]  lmn1,lmn2,lmn3,lmn4
 *              Exponents of x, y, and z that signify angular momentum. Required
 *              to be 3 elements.
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four centers
 * \param [in]  cb
 *              The callback function that computes integrals in interval arithmetic
 */
void mirp_integral4_single_target_str(arb_t integral,
                                           const int * lmn1, const char ** A, const char * alpha1,
                                           const int * lmn2, const char ** B, const char * alpha2,
                                           const int * lmn3, const char ** C, const char * alpha3,
                                           const int * lmn4, const char ** D, const char * alpha4,
                                           slong target_prec, cb_integral4_single cb);


/*! \brief Computes a single 4-center integral to double precision using interval arithmetic
 *
 * This function takes double precision as input and returns double precision
 * as output. Internally, it uses interval arithmetic to ensure that no
 * precision is lost
 *
 * \param [out] integral
 *              Resulting integral output
 * \param [in]  lmn1,lmn2,lmn3,lmn4
 *              Exponents of x, y, and z that signify angular momentum. Required
 *              to be 3 elements.
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four centers
 * \param [in]  cb
 *              The callback function that computes integrals in interval arithmetic
 */
void mirp_integral4_single_exact(double * integral,
                                 const int * lmn1, const double * A, double alpha1,
                                 const int * lmn2, const double * B, double alpha2,
                                 const int * lmn3, const double * C, double alpha3,
                                 const int * lmn4, const double * D, double alpha4,
                                 cb_integral4_single cb);



#define MIRP_WRAP_SHELL4(name) \
    static inline \
    void mirp_##name(arb_t integral, \
                     int am1, arb_srcptr A, int nprim1, int ngeneral1, arb_srcptr alpha1, arb_srcptr coeff1, \
                     int am2, arb_srcptr B, int nprim2, int ngeneral2, arb_srcptr alpha2, arb_srcptr coeff2, \
                     int am3, arb_srcptr C, int nprim3, int ngeneral3, arb_srcptr alpha3, arb_srcptr coeff3, \
                     int am4, arb_srcptr D, int nprim4, int ngeneral4, arb_srcptr alpha4, arb_srcptr coeff4, \
                     slong working_prec) \
    { \
        mirp_loop_shell4(integral, \
                         am1, A, nprim1, ngeneral1, alpha1, coeff1, \
                         am2, B, nprim2, ngeneral2, alpha2, coeff2, \
                         am3, C, nprim3, ngeneral3, alpha3, coeff3, \
                         am4, D, nprim4, ngeneral4, alpha4, coeff4, \
                         working_prec, mirp_##name##_single); \
    }

#define MIRP_WRAP_SHELL4_D(name) \
    static inline \
    void mirp_##name##_d(double * integral, \
                         int am1, const double * A, int nprim1, int ngeneral1, const double * alpha1, const double * coeff1, \
                         int am2, const double * B, int nprim2, int ngeneral2, const double * alpha2, const double * coeff2, \
                         int am3, const double * C, int nprim3, int ngeneral3, const double * alpha3, const double * coeff3, \
                         int am4, const double * D, int nprim4, int ngeneral4, const double * alpha4, const double * coeff4) \
    { \
        mirp_loop_shell4_d(integral, \
                           am1, A, nprim1, ngeneral1, alpha1, coeff1, \
                           am2, B, nprim2, ngeneral2, alpha2, coeff2, \
                           am3, C, nprim3, ngeneral3, alpha3, coeff3, \
                           am4, D, nprim4, ngeneral4, alpha4, coeff4, \
                           mirp_##name##_single_d); \
    }


#define MIRP_WRAP_SINGLE4_TARGET_PREC(name) \
    static inline \
    int mirp_##name##_single_target(arb_t integral, \
                                         const int * lmn1, arb_srcptr A, const arb_t alpha1, \
                                         const int * lmn2, arb_srcptr B, const arb_t alpha2, \
                                         const int * lmn3, arb_srcptr C, const arb_t alpha3, \
                                         const int * lmn4, arb_srcptr D, const arb_t alpha4, \
                                         slong target_prec) \
    { \
        return mirp_integral4_single_target(integral, \
                                                 lmn1, A, alpha1, \
                                                 lmn2, B, alpha2, \
                                                 lmn3, C, alpha3, \
                                                 lmn4, D, alpha4, \
                                                 target_prec, \
                                                 mirp_##name##_single); \
    }


#define MIRP_WRAP_SINGLE4_TARGET_PREC_STR(name) \
    static inline \
    void mirp_##name##_single_target_str(arb_t integral, \
                                              const int * lmn1, const char ** A, const char * alpha1, \
                                              const int * lmn2, const char ** B, const char * alpha2, \
                                              const int * lmn3, const char ** C, const char * alpha3, \
                                              const int * lmn4, const char ** D, const char * alpha4, \
                                              slong target_prec) \
    { \
        mirp_integral4_single_target_str(integral, \
                                              lmn1, A, alpha1, \
                                              lmn2, B, alpha2, \
                                              lmn3, C, alpha3, \
                                              lmn4, D, alpha4, \
                                              target_prec, \
                                              mirp_##name##_single); \
    }

#define MIRP_WRAP_SINGLE4_EXACT(name) \
    static inline \
    void mirp_##name##_single_exact(double * integral, \
                                    const int * lmn1, const double * A, double alpha1, \
                                    const int * lmn2, const double * B, double alpha2, \
                                    const int * lmn3, const double * C, double alpha3, \
                                    const int * lmn4, const double * D, double alpha4) \
    { \
        mirp_integral4_single_exact(integral, \
                                    lmn1, A, alpha1, \
                                    lmn2, B, alpha2, \
                                    lmn3, C, alpha3, \
                                    lmn4, D, alpha4, \
                                    mirp_##name##_single); \
    }





#ifdef __cplusplus
}
#endif

