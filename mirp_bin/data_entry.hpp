/*! \file
 *
 * \brief Common data structures for testing
 */

#pragma once

#include <array>
#include <string>
#include <vector>

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


struct integral_single_data_entry
{
    std::vector<gaussian_single> g;  //!< Shell array for which the data is calculated
    std::string integral;            //!< Actual data (integral) for this gaussian
};


struct integral_single_data
{
    long ndigits;                        //!< The number of digits of accuracy in the file 
    std::string header;                  //!< Header or comments about the test
    std::vector<integral_single_data_entry> values;  //!< Actual data for the integrals
};


} // close namespace mirp

