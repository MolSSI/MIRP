/*! \file
 *
 * \brief Shell structures used in testing
 */

#pragma once

#include <vector>
#include <string>
#include <array>

namespace mirp {

    struct gaussian_single
    {
        std::array<int, 3> lmn;         //!< AM exponents
        std::array<std::string, 3> xyz; //!< Coordinates (in bohr)
        std::string alpha;              //!< Exponent of the gaussian
    };

    struct gaussian_shell
    {
        int am;                          //!< Angular momentum
        int ngeneral;                    //!< Number of general contractions
        std::array<std::string, 3> xyz;  //!< Coordinates (in bohr)
        std::vector<std::string> alpha;  //!< Exponents of the gaussians
        std::vector<std::string> coef;   //!< Contraction coefficients (unnormalized)
    };

} // close namespace mirp

