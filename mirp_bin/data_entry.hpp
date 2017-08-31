/*! \file
 *
 * \brief Common data structures (shells, entries, etc)
 */

#pragma once

#include <array>
#include <vector>
#include <string>

namespace mirp {


/*! \brief A shell of cartesian gaussian functions (double precision) 
 *
 * This shell holds information as double precision
 */
struct gaussian_shell
{
    int Z;                      //!< Z-number of the center
    int am;                     //!< Angular momentum
    int nprim;                  //!< Number of primitives (segmented contraction)
    int ngeneral;               //!< Number of general contractions
    std::array<double, 3> xyz;  //!< Coordinates (in bohr)
    std::vector<double> alpha;  //!< Exponents of the gaussians
    std::vector<double> coeff;  //!< Contraction coefficients (unnormalized)
};


/*! \brief A shell of cartesian gaussian functions (strings) 
 *
 * This shell holds information as a string to allow for higher precision
 */
struct gaussian_shell_str
{
    int Z;                           //!< Z-number of the center
    int am;                          //!< Angular momentum
    int nprim;                       //!< Number of primitives (segmented contraction)
    int ngeneral;                    //!< Number of general contractions
    std::array<std::string, 3> xyz;  //!< Coordinates (in bohr)
    std::vector<std::string> alpha;  //!< Exponents of the gaussians
    std::vector<std::string> coeff;  //!< Contraction coefficients (unnormalized)
};


/*! \brief A single cartesian gaussian function */
struct gaussian_single
{
    std::array<int, 3> lmn;         //!< AM exponents
    std::array<std::string, 3> xyz; //!< Coordinates (in bohr)
    std::string alpha;              //!< Exponent of the gaussian
};


/*! \brief Data for a single cartesian integral
 *
 * The length of integral_single_data_entry::g represents the number of
 * centers for the integral.
 */
struct integral_single_data_entry
{
    std::vector<gaussian_single> g;  //!< Shell array for which the integrals calculated
    std::string integral;            //!< Actual data (integral) for this gaussian
};


/*! \brief Contents of a data file for a single cartesian integral
 *
 * The data file holds data for multiple integrals, plus a descriptive
 * header and some metadata.
 */
struct integral_single_data
{
    long ndigits;                                    //!< The number of digits of accuracy in the file
    std::string header;                              //!< Header or comments about the test
    std::vector<integral_single_data_entry> entries; //!< Actual data for the integrals
};



/*! \brief Data for a single contracted integral */
struct integral_data_entry
{
    std::vector<gaussian_shell_str> g;      //!< Shells for which the integrals are calculated
    std::vector<int> idx;               //!< Indices of the shells
    std::vector<std::string> integrals; //!< Computed integrals
};


/*! \brief Contents of a data file for contracted integrals
 *
 * The data file holds data for multiple integrals, plus a descriptive
 * header and some metadata.
 */
struct integral_data
{
    long ndigits;                             //!< The number of digits of accuracy in the file
    std::string header;                       //!< Header or comments about the test
    std::vector<integral_data_entry> entries; //!< Actual data for the integrals
};


/*! \brief Creates a basis from an XYZ file and a basis set file
 *
 * \throw std::runtime_error if there is a problem reading or parsing
 *        the files
 *
 * \param [in] xyzfile Path to a file containing atomic coordinates
 * \param [in] basfile Path to an NWChem-formatted basis set file
 * \return A vector of shells obtained from applying the basis to the coordinate file
 */
std::vector<gaussian_shell> read_construct_basis(const std::string & xyzfile,
                                                 const std::string & basfile);


} // close namespace mirp

