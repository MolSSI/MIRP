/*! \file
 *
 * \brief Some useful wrapping of functionality
 */

#pragma once

#include "mirp/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Compute all cartesian integrals of a contracted shell quartet
 *         for an integral (four-center, double precision)
 *
 * This function takes in a pointer to a function that computes single,
 * primitive cartesian integrals, and uses it to compute all the cartesian
 * components for an contracted shell quartet.
 *
 * \param [out] integral
 *              Output for the computed integral
 * \param [in]  am1,am2,am3,am4
 *              Angular momentum for the four-centers
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four-centers (each of length 3)
 * \param [in]  nprim1,nprim2,nprim3,nprim4
 *              Number of primitive gaussians for each shell
 * \param [in]  ngen1,ngen2,ngen3,ngen4
 *              Number of general contractions for each shell
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the primitive gaussians on the four-centers
 *              (of lengths \p nprim1, \p nprim2, \p nprim3, \p nprim4 respectively)
 * \param [in]  coeff1,coeff2,coeff3,coeff4
 *              Coefficients for all primitives and for all general contractions
 *              for each shell (of lengths \p nprim1 * \p ngen1, \p nprim2 * \p ngen2,
 *              \p nprim3 * \p ngen3, \p nprim4 * \p ngen4 respectively)
 * \param [in]  cb
 *              Function that computes a single cartesian four-center integral
 *              with interval arithmetic
 */
void mirp_loop_shell4_d(double * integral,
                        int am1, const double * A, int nprim1, int ngen1, const double * alpha1, const double * coeff1,
                        int am2, const double * B, int nprim2, int ngen2, const double * alpha2, const double * coeff2,
                        int am3, const double * C, int nprim3, int ngen3, const double * alpha3, const double * coeff3,
                        int am4, const double * D, int nprim4, int ngen4, const double * alpha4, const double * coeff4,
                        cb_integral4_single_d cb);


/*! \brief Compute all cartesian integrals of a contracted shell quartet
 *         for an integral (four-center, interval arithmetic)
 *
 * \copydetails mirp_loop_shell4_d
 *
 * \param [in]  working_prec
 *              The working precision (binary digits/bits) to use
 *              in the calculation
 */
void mirp_loop_shell4(arb_ptr integral,
                      int am1, arb_srcptr A, int nprim1, int ngen1, arb_srcptr alpha1, arb_srcptr coeff1,
                      int am2, arb_srcptr B, int nprim2, int ngen2, arb_srcptr alpha2, arb_srcptr coeff2,
                      int am3, arb_srcptr C, int nprim3, int ngen3, arb_srcptr alpha3, arb_srcptr coeff3,
                      int am4, arb_srcptr D, int nprim4, int ngen4, arb_srcptr alpha4, arb_srcptr coeff4,
                      slong working_prec, cb_integral4_single cb);


/*! \brief Compute a single integral to a target precision (four-center)
 *
 * This function seeks to calculate the integral to a target precision. This works
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
 *              Output for the computed integral
 * \param [in]  lmn1,lmn2,lmn3,lmn4
 *              Exponents of x, y, and z that signify angular momentum. Required
 *              to be 3 elements.
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four-centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four-centers
 * \param [in]  target_prec
 *              Desired precision (in bits)
 * \param [in]  cb
 *              Function that computes a single cartesian four-center integral
 *              with interval arithmetic
 * \return Nonzero if there is a problem with obtaining the target precision (usually
 *         due to a lack of precision in the input parameters)
 */
int mirp_integral4_single_target(arb_t integral,
                                 const int * lmn1, arb_srcptr A, const arb_t alpha1,
                                 const int * lmn2, arb_srcptr B, const arb_t alpha2,
                                 const int * lmn3, arb_srcptr C, const arb_t alpha3,
                                 const int * lmn4, arb_srcptr D, const arb_t alpha4,
                                 slong target_prec, cb_integral4_single cb);


/*! \brief Computes all cartesian integrals of a contracted shell quartet
 *         to a target precision (four-center)
 *
 * This function seeks to calculate the integrals to a target precision. This works
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
 * \param [out] integrals
 *              Output for the computed integrals
 * \param [in]  am1,am2,am3,am4
 *              Angular momentum for the four-centers
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four-centers (each of length 3)
 * \param [in]  nprim1,nprim2,nprim3,nprim4
 *              Number of primitive gaussians for each shell
 * \param [in]  ngen1,ngen2,ngen3,ngen4
 *              Number of general contractions for each shell
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the primitive gaussians on the four-centers
 *              (of lengths \p nprim1, \p nprim2, \p nprim3, \p nprim4 respectively)
 * \param [in]  coeff1,coeff2,coeff3,coeff4
 *              Coefficients for all primitives and for all general contractions
 *              for each shell (of lengths \p nprim1 * \p ngen1, \p nprim2 * \p ngen2,
 *              \p nprim3 * \p ngen3, \p nprim4 * \p ngen4 respectively)
 * \param [in]  target_prec
 *              Desired target precision (in bits)
 * \param [in]  cb
 *              Function that computes a single cartesian four-center integral
 *              with interval arithmetic
 * \return Nonzero if there is a problem with obtaining the target precision (usually
 *         due to a lack of precision in the input parameters)
 */
int mirp_integral4_target(arb_ptr integrals,
                          int am1, arb_srcptr A, int nprim1, int ngen1, arb_srcptr alpha1, arb_srcptr coeff1,
                          int am2, arb_srcptr B, int nprim2, int ngen2, arb_srcptr alpha2, arb_srcptr coeff2,
                          int am3, arb_srcptr C, int nprim3, int ngen3, arb_srcptr alpha3, arb_srcptr coeff3,
                          int am4, arb_srcptr D, int nprim4, int ngen4, arb_srcptr alpha4, arb_srcptr coeff4,
                          slong target_prec, cb_integral4 cb);


/*! \brief Compute a single 4-center integral to a target precision (string input)
 *
 * \copydetails mirp_integral4_single_target
 *
 * Similar to mirp_integral4_single_target, except that the input parameters are strings.
 * The strings are interpreted to have infinite precision. Because of this, this function
 * should not fail due to inability to reach the target precision.
 */
void mirp_integral4_single_target_str(arb_t integral,
                                      const int * lmn1, const char ** A, const char * alpha1,
                                      const int * lmn2, const char ** B, const char * alpha2,
                                      const int * lmn3, const char ** C, const char * alpha3,
                                      const int * lmn4, const char ** D, const char * alpha4,
                                      slong target_prec, cb_integral4_single cb);


/*! \brief Compute a single 4-center integral to a target precision (string input)
 *
 * \copydetails mirp_integral4_target
 *
 * Similar to mirp_integral4_target, except that the input parameters are strings.
 * The strings are interpreted to have infinite precision. Because of this, this function
 * should not fail due to inability to reach the target precision.
 */
void mirp_integral4_target_str(arb_ptr integrals,
                               int am1, const char ** A, int nprim1, int ngen1, const char ** alpha1, const char ** coeff1,
                               int am2, const char ** B, int nprim2, int ngen2, const char ** alpha2, const char ** coeff2,
                               int am3, const char ** C, int nprim3, int ngen3, const char ** alpha3, const char ** coeff3,
                               int am4, const char ** D, int nprim4, int ngen4, const char ** alpha4, const char ** coeff4,
                               slong target_prec, cb_integral4 cb);


/*! \brief Computes a single integral to double precision using interval arithmetic (four-center)
 *
 * This function takes double precision as input and returns double precision
 * as output. Internally, it uses interval arithmetic to ensure that no
 * precision is lost
 *
 * \param [out] integral
 *              Output for the computed integral
 * \param [in]  lmn1,lmn2,lmn3,lmn4
 *              Exponents of x, y, and z that signify angular momentum. Required
 *              to be 3 elements.
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four-centers (each of length 3)
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the gaussian on the four-centers
 * \param [in]  cb
 *              Function that computes a single cartesian four-center integral
 *              with interval arithmetic
 */
void mirp_integral4_single_exact(double * integral,
                                 const int * lmn1, const double * A, double alpha1,
                                 const int * lmn2, const double * B, double alpha2,
                                 const int * lmn3, const double * C, double alpha3,
                                 const int * lmn4, const double * D, double alpha4,
                                 cb_integral4_single_target cb);


/*! \brief Compute all cartesian integrals of a contracted shell quartet
 *         for an integral to exact double precision (four-center)
 *
 * This function takes double precision as input and returns double precision
 * as output. Internally, it uses interval arithmetic to ensure that no
 * precision is lost
 *
 * \param [out] integrals
 *              Output for the computed integrals
 * \param [in]  am1,am2,am3,am4
 *              Angular momentum for the four-centers
 * \param [in]  A,B,C,D
 *              XYZ coordinates of the four-centers (each of length 3)
 * \param [in]  nprim1,nprim2,nprim3,nprim4
 *              Number of primitive gaussians for each shell
 * \param [in]  ngen1,ngen2,ngen3,ngen4
 *              Number of general contractions for each shell
 * \param [in]  alpha1,alpha2,alpha3,alpha4
 *              Exponents of the primitive gaussians on the four-centers
 *              (of lengths \p nprim1, \p nprim2, \p nprim3, \p nprim4 respectively)
 * \param [in]  coeff1,coeff2,coeff3,coeff4
 *              Coefficients for all primitives and for all general contractions
 *              for each shell (of lengths \p nprim1 * \p ngen1, \p nprim2 * \p ngen2,
 *              \p nprim3 * \p ngen3, \p nprim4 * \p ngen4 respectively)
 * \param [in]  cb
 *              Function that computes a single cartesian four-center integral
 *              with interval arithmetic
 */
void mirp_integral4_exact(double * integrals,
                          int am1, const double * A, int nprim1, int ngen1, const double * alpha1, const double * coeff1,
                          int am2, const double * B, int nprim2, int ngen2, const double * alpha2, const double * coeff2,
                          int am3, const double * C, int nprim3, int ngen3, const double * alpha3, const double * coeff3,
                          int am4, const double * D, int nprim4, int ngen4, const double * alpha4, const double * coeff4,
                          cb_integral4_target cb);


/*! \brief Create a function that computes all cartesian integrals
 *         of a contracted shell quartet (four-center, interval arithmetic)
 *
 *  A function computing single cartesian integrals is expected to exist and
 *  be named `mirp_{name}_single`
 *
 *  The created function is named `mirp_{name}`.
 *
 *  \sa mirp_loop_shell4
 */
#define MIRP_WRAP_SHELL4(name) \
    static inline \
    void mirp_##name(arb_t integral, \
                     int am1, arb_srcptr A, int nprim1, int ngen1, arb_srcptr alpha1, arb_srcptr coeff1, \
                     int am2, arb_srcptr B, int nprim2, int ngen2, arb_srcptr alpha2, arb_srcptr coeff2, \
                     int am3, arb_srcptr C, int nprim3, int ngen3, arb_srcptr alpha3, arb_srcptr coeff3, \
                     int am4, arb_srcptr D, int nprim4, int ngen4, arb_srcptr alpha4, arb_srcptr coeff4, \
                     slong working_prec) \
    { \
        mirp_loop_shell4(integral, \
                         am1, A, nprim1, ngen1, alpha1, coeff1, \
                         am2, B, nprim2, ngen2, alpha2, coeff2, \
                         am3, C, nprim3, ngen3, alpha3, coeff3, \
                         am4, D, nprim4, ngen4, alpha4, coeff4, \
                         working_prec, mirp_##name##_single); \
    }


/*! \brief Create a function that computes all cartesian integrals
 *         of a contracted shell quartet (four-center, double precision)
 *
 *  A function computing single cartesian integrals in double precision
 *  is expected to exist and be named `mirp_{name}_single_d`
 *
 *  The created function is named `mirp_{name}_d`.
 *
 *  \sa mirp_loop_shell4_d
 */
#define MIRP_WRAP_SHELL4_D(name) \
    static inline \
    void mirp_##name##_d(double * integral, \
                         int am1, const double * A, int nprim1, int ngen1, const double * alpha1, const double * coeff1, \
                         int am2, const double * B, int nprim2, int ngen2, const double * alpha2, const double * coeff2, \
                         int am3, const double * C, int nprim3, int ngen3, const double * alpha3, const double * coeff3, \
                         int am4, const double * D, int nprim4, int ngen4, const double * alpha4, const double * coeff4) \
    { \
        mirp_loop_shell4_d(integral, \
                           am1, A, nprim1, ngen1, alpha1, coeff1, \
                           am2, B, nprim2, ngen2, alpha2, coeff2, \
                           am3, C, nprim3, ngen3, alpha3, coeff3, \
                           am4, D, nprim4, ngen4, alpha4, coeff4, \
                           mirp_##name##_single_d); \
    }


/*! \brief Create a function that computes single cartesian integrals
 *         to a target precision (four-center) 
 *
 *  A function computing single cartesian integrals is expected to exist and
 *  be named `mirp_{name}_single`
 *
 *  The created function is named `mirp_{name}_single_target`.
 *
 *  \sa mirp_integral4_single_target
 */
#define MIRP_WRAP_SINGLE4_TARGET(name) \
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


/*! \brief Create a function that computes all cartesian integrals
 *         of a contracted shell quartet to a target precision (four-center)
 *
 *  A function computing all cartesian integrals of a contracted shell quartet
 *  in interval arithmetic is expected to exist and be named `mirp_{name}`
 *
 *  The created function is named `mirp_{name}_target`.
 *
 *  \sa mirp_integral4_target
 */
#define MIRP_WRAP_SHELL4_TARGET(name) \
    static inline \
    int mirp_##name##_target(arb_ptr integrals, \
                             int am1, arb_srcptr A, int nprim1, int ngen1, arb_srcptr alpha1, arb_srcptr coeff1, \
                             int am2, arb_srcptr B, int nprim2, int ngen2, arb_srcptr alpha2, arb_srcptr coeff2, \
                             int am3, arb_srcptr C, int nprim3, int ngen3, arb_srcptr alpha3, arb_srcptr coeff3, \
                             int am4, arb_srcptr D, int nprim4, int ngen4, arb_srcptr alpha4, arb_srcptr coeff4, \
                             slong target_prec) \
    { \
        return mirp_integral4_target(integrals, \
                                     am1, A, nprim1, ngen1, alpha1, coeff1, \
                                     am2, B, nprim2, ngen2, alpha2, coeff2, \
                                     am3, C, nprim3, ngen3, alpha3, coeff3, \
                                     am4, D, nprim4, ngen4, alpha4, coeff4, \
                                     target_prec, \
                                     mirp_##name); \
    }


/*! \brief Create a function that computes single cartesian integrals
 *         to a target precision (four-center, string inputs) 
 *
 *  A function computing single cartesian integrals is expected to exist and
 *  be named `mirp_{name}_single`
 *
 *  The created function is named `mirp_{name}_single_target_str`.
 *
 *  \sa mirp_integral4_single_target_str
 */
#define MIRP_WRAP_SINGLE4_TARGET_STR(name) \
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


/*! \brief Create a function that computes all cartesian integrals
 *         of a contracted shell quartet to a target precision
 *         (four-center, string inputs)
 *
 *  A function computing all cartesian integrals of a contracted shell quartet
 *  in interval arithmetic is expected to exist and be named `mirp_{name}`
 *
 *  The created function is named `mirp_{name}_target_str`.
 *
 *  \sa mirp_integral4_target_str
 */
#define MIRP_WRAP_SHELL4_TARGET_STR(name) \
    static inline \
    void mirp_##name##_target_str(arb_ptr integrals, \
                                  int am1, const char ** A, int nprim1, int ngen1, const char ** alpha1, const char ** coeff1, \
                                  int am2, const char ** B, int nprim2, int ngen2, const char ** alpha2, const char ** coeff2, \
                                  int am3, const char ** C, int nprim3, int ngen3, const char ** alpha3, const char ** coeff3, \
                                  int am4, const char ** D, int nprim4, int ngen4, const char ** alpha4, const char ** coeff4, \
                                  slong target_prec) \
    { \
        mirp_integral4_target_str(integrals, \
                                  am1, A, nprim1, ngen1, alpha1, coeff1, \
                                  am2, B, nprim2, ngen2, alpha2, coeff2, \
                                  am3, C, nprim3, ngen3, alpha3, coeff3, \
                                  am4, D, nprim4, ngen4, alpha4, coeff4, \
                                  target_prec, \
                                  mirp_##name); \
    }


/*! \brief Create a function that computes single cartesian integrals
 *         to exact double precision (four-center)
 *
 *  A function computing single cartesian integrals to a target precision is
 *  expected to exist and be named `mirp_{name}_single_target`
 *
 *  The created function is named `mirp_{name}_single_exact`.
 *
 *  \sa mirp_integral4_single_exact
 */
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
                                    mirp_##name##_single_target); \
    }


/*! \brief Create a function that computes all cartesian integrals
 *         of a contracted shell quartet to exact double precision (four-center)
 *
 *  A function computing all cartesian integrals of a contracted shell quartet
 *  to a target precision is expected to exist and be named `mirp_{name}_target`
 *
 *  The created function is named `mirp_{name}_exact`.
 *
 *  \sa mirp_integral4_exact
 */
#define MIRP_WRAP_SHELL4_EXACT(name) \
    static inline \
    void mirp_##name##_exact(double * integrals, \
                             int am1, const double * A, int nprim1, int ngen1, const double * alpha1, const double * coeff1, \
                             int am2, const double * B, int nprim2, int ngen2, const double * alpha2, const double * coeff2, \
                             int am3, const double * C, int nprim3, int ngen3, const double * alpha3, const double * coeff3, \
                             int am4, const double * D, int nprim4, int ngen4, const double * alpha4, const double * coeff4) \
    { \
        mirp_integral4_exact(integrals, \
                             am1, A, nprim1, ngen1, alpha1, coeff1, \
                             am2, B, nprim2, ngen2, alpha2, coeff2, \
                             am3, C, nprim3, ngen3, alpha3, coeff3, \
                             am4, D, nprim4, ngen4, alpha4, coeff4, \
                             mirp_##name##_target); \
    }



#ifdef __cplusplus
}
#endif

