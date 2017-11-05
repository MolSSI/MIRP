/*! \file
 *
 * \brief Structures that help with calling callbacks
 */

#pragma once

#include <array>
#include <vector>
#include <mirp/typedefs.h>

namespace mirp {

template<>
struct callback_helper<4>
{
    typedef cb_integral4                cb_type;
    typedef cb_integral4_str            cb_str_type;
    typedef cb_integral4_exact          cb_exact_type;

    typedef cb_integral4_single         cb_single_type;
    typedef cb_integral4_single_str     cb_single_str_type;
    typedef cb_integral4_single_exact   cb_single_exact_type;

    static void 
    call_str(arb_ptr integrals,
             std::array<int, 4> & am,
             std::array<std::array<const char *, 3>, 4> & xyz,
             std::array<int, 4> & nprim,
             std::array<int, 4> & ngeneral,
             std::array<std::vector<const char *>, 4> & alpha,
             std::array<std::vector<const char *>, 4> & coeff,
             slong working_prec,
             cb_str_type cb)
    {
        cb(integrals,
           am[0], xyz[0].data(), nprim[0], ngeneral[0], alpha[0].data(), coeff[0].data(),
           am[1], xyz[1].data(), nprim[1], ngeneral[1], alpha[1].data(), coeff[1].data(),
           am[2], xyz[2].data(), nprim[2], ngeneral[2], alpha[2].data(), coeff[2].data(),
           am[3], xyz[3].data(), nprim[3], ngeneral[3], alpha[3].data(), coeff[3].data(),
           working_prec);
    }


    static void
    call(arb_ptr integrals,
         std::array<int, 4> & am,
         std::array<arb_ptr, 4> & xyz,
         std::array<int, 4> & nprim,
         std::array<int, 4> & ngeneral,
         std::array<arb_ptr, 4> & alpha,
         std::array<arb_ptr, 4> & coeff,
         slong working_prec,
         cb_type cb)
    {
        cb(integrals,
           am[0], xyz[0], nprim[0], ngeneral[0], alpha[0], coeff[0],
           am[1], xyz[1], nprim[1], ngeneral[1], alpha[1], coeff[1],
           am[2], xyz[2], nprim[2], ngeneral[2], alpha[2], coeff[2],
           am[3], xyz[3], nprim[3], ngeneral[3], alpha[3], coeff[3],
           working_prec);
    }


    static void
    call_exact(double * integrals,
               std::array<int, 4> & am,
               std::array<std::array<double, 3>, 4> & xyz,
               std::array<int, 4> & nprim,
               std::array<int, 4> & ngeneral,
               std::array<std::vector<double>, 4> & alpha,
               std::array<std::vector<double>, 4> & coeff,
               cb_exact_type cb)
    {
        cb(integrals,
           am[0], xyz[0].data(), nprim[0], ngeneral[0], alpha[0].data(), coeff[0].data(),
           am[1], xyz[1].data(), nprim[1], ngeneral[1], alpha[1].data(), coeff[1].data(),
           am[2], xyz[2].data(), nprim[2], ngeneral[2], alpha[2].data(), coeff[2].data(),
           am[3], xyz[3].data(), nprim[3], ngeneral[3], alpha[3].data(), coeff[3].data());
    }


    static void
    call_single_arb(arb_t integral,
                    std::array<std::array<int, 3>, 4> & lmn,
                    std::array<arb_ptr, 4> & xyz,
                    std::array<arb_t, 4> & alpha,
                    slong working_prec,
                    cb_single_type cb)
    {
        cb(integral,
           lmn[0].data(), xyz[0], alpha[0],
           lmn[1].data(), xyz[1], alpha[1],
           lmn[2].data(), xyz[2], alpha[2],
           lmn[3].data(), xyz[3], alpha[3],
           working_prec);
    }


    static void
    call_single_str(arb_t integral,
                    std::array<std::array<int, 3>, 4> & lmn,
                    std::array<std::array<const char *, 3>, 4> & xyz,
                    std::array<const char *, 4> & alpha,
                    slong working_prec,
                    cb_single_str_type cb)
    {
        cb(integral,
           lmn[0].data(), xyz[0].data(), alpha[0],
           lmn[1].data(), xyz[1].data(), alpha[1],
           lmn[2].data(), xyz[2].data(), alpha[2],
           lmn[3].data(), xyz[3].data(), alpha[3],
           working_prec);
    }


    static void
    call_single_exact(double * integral,
                      std::array<std::array<int, 3>, 4> & lmn,
                      std::array<std::array<double, 3>, 4> & xyz,
                      std::array<double, 4> & alpha,
                      cb_single_exact_type cb)
    {
        cb(integral,
           lmn[0].data(), xyz[0].data(), alpha[0],
           lmn[1].data(), xyz[1].data(), alpha[1],
           lmn[2].data(), xyz[2].data(), alpha[2],
           lmn[3].data(), xyz[3].data(), alpha[3]);
    }

};

} // close namespace mirp

