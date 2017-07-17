/*! \file
 *
 * \brief Common data structures for testing
 */

#pragma once

#include <array>
#include <string>
#include <vector>

namespace mirp {

/*! \brief A single cartesian gaussian function */
struct gaussian_single
{
    std::array<int, 3> lmn;         //!< AM exponents
    std::array<std::string, 3> xyz; //!< Coordinates (in bohr)
    std::string alpha;              //!< Exponent of the gaussian
};


/*! \brief A shell of cartesian gaussian functions */
struct gaussian_shell
{
    int am;                          //!< Angular momentum
    int ngeneral;                    //!< Number of general contractions
    std::array<std::string, 3> xyz;  //!< Coordinates (in bohr)
    std::vector<std::string> alpha;  //!< Exponents of the gaussians
    std::vector<std::string> coef;   //!< Contraction coefficients (unnormalized)
};


/*! \brief Data for a single cartesian integral
 *
 * The length of integral_single_data_entry::g represents the number of
 * centers for the integral. 
 */
struct integral_single_data_entry
{
    std::vector<gaussian_single> g;  //!< Shell array for which the data is calculated
    std::string integral;            //!< Actual data (integral) for this gaussian
};


/*! \brief Contents of a data file for a single cartesian integral
 *
 * The data file holds data for multiple integrals, plus a descriptive
 * header and some metadata.
 */
struct integral_single_data
{
    long ndigits;                        //!< The number of digits of accuracy in the file 
    std::string header;                  //!< Header or comments about the test
    std::vector<integral_single_data_entry> values;  //!< Actual data for the integrals
};


} // close namespace mirp

