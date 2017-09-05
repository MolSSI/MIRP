/*! \file
 *
 * \brief Construction of a basis from coordinates and a basis set file
 */


#include "mirp_bin/reffile_io.hpp"
#include "mirp_bin/test_common.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>


/*! Factor for converting Angstrom to Bohr
 *
 * Source: NIST 2014 CODATA recommended values
 */
const double ang_to_bohr = 1.0/0.52917721067;


namespace mirp {

// Anonymous namespace for some helper functions and structs
namespace {

/*! \brief Helper struct representing an atom */
struct atom
{
    int Z;                     //!< Atomic Z number
    std::array<double, 3> xyz; //!< Coordinates of the atom 
};


/*! \brief Reads a simple XYZ file
 *
 * The first line is the number of atoms. The second is
 * an arbitrary comment line.
 *
 * The remaining lines contain the element symbol as a string,
 * followed by the x, y, and z coordinates in angstroms.
 */
std::vector<atom> read_xyz(const std::string & xyzfile)
{
    std::vector<atom> atom_vec;

    std::ifstream infile(xyzfile);

    if(!infile.good())
        throw std::runtime_error("Cannot open xyz file for reading");

    std::string line;
    std::getline(infile, line); // skip number of atoms
    std::getline(infile, line); // skip comment line

    while(std::getline(infile, line))
    {
        line = trim(line);
        atom a;
        std::stringstream ss(line);
        std::string element;
        ss >> element >> a.xyz[0] >> a.xyz[1] >> a.xyz[2];
        a.xyz[0] *= ang_to_bohr;
        a.xyz[1] *= ang_to_bohr;
        a.xyz[2] *= ang_to_bohr;
        a.Z = element_to_z(element);
        atom_vec.push_back(std::move(a));
    }

    return atom_vec;
}


/*! \brief Reads an NWChem-style basis set file and returns a map
 *
 * The map translates an atomic Z-number to a vector of shells for that element
 *
 * \note The shells returned do not have the xyz coordinates filled in
 *
 * \todo Needs hardening and more checking
 */
std::map<int, std::vector<gaussian_shell>>
read_basfile(const std::string & basfile)
{
    std::map<int, std::vector<gaussian_shell>> ret;

    std::ifstream infile(basfile);

    if(!infile.good())
        throw std::runtime_error("Cannot open basis set file for reading");

    std::string line;
    std::getline(infile, line);
    line = trim(line);

    while(infile)
    {
        // blank or comment line
        if(line.length() == 0 || line[0] == '#' ||
           (line.length() >= 5 && line.substr(0, 5) == "BASIS"))
        {
            std::getline(infile, line);
            line = trim(line);
            continue;
        }

        // end of the file
        if(line.length() >= 3 && line.substr(0, 3) == "END")
            break;

        if(!std::isalpha(line[0]))
            throw std::runtime_error("Malformed basis set file");

        gaussian_shell sh;
        sh.ngeneral = 0;

        std::vector<std::string> coeff_tmp;

        std::string element_str, am_str;
        std::stringstream ss(line);
        ss >> element_str >> am_str;

        int Z = element_to_z(element_str); 

        std::getline(infile, line);
        line = trim(line);

        while(infile && (line.length() && std::isdigit(line[0])))
        {
            auto c = split(line);
            sh.alpha.push_back(std::strtod(c[0].c_str(), nullptr)); 
            coeff_tmp.insert(coeff_tmp.end(), c.begin()+1, c.end());

            if(sh.ngeneral != 0 && sh.ngeneral != static_cast<int>(c.size()-1))
                throw std::runtime_error("Inconsistent number of general contractions");
            else
                sh.ngeneral = static_cast<int>(c.size())-1;

            std::getline(infile, line); 
            line = trim(line);
        } 

        sh.nprim = sh.alpha.size();

        // reorder the coefficients
        sh.coeff.resize(coeff_tmp.size());
        for(int i = 0; i < sh.nprim; i++)
        for(int j = 0; j < sh.ngeneral; j++)
            sh.coeff[j*sh.nprim+i] = std::strtod(coeff_tmp[i*sh.ngeneral+j].c_str(), nullptr);


        // handle special "combined" general contractions
        if(am_str.length() > 1)
        {
            for(int i = 0; i < sh.ngeneral; i++)
            {
                gaussian_shell sh2;
                sh2.am = amchar_to_int(am_str.at(i));
                sh2.nprim = sh.nprim;
                sh2.ngeneral = 1;
                sh2.alpha = sh.alpha;
                sh2.coeff.insert(sh2.coeff.begin(),
                                 sh.coeff.begin() + i*sh.nprim,
                                 sh.coeff.begin() + (i+1)*sh.nprim);

                ret[Z].push_back(sh2);
            }
        }
        else
        {
            sh.am = amchar_to_int(am_str.at(0));
            ret[Z].push_back(sh);
        }
    }

    return ret;
    
}

} // close anonymous namespace



/*! \brief Create a basis from an XYZ file and a basis set file
 *
 * \throw std::runtime_error if there is a problem reading or parsing
 *        the files
 */
std::vector<gaussian_shell> read_construct_basis(const std::string & xyzfile,
                                                 const std::string & basfile)
{
    const auto atoms = read_xyz(xyzfile);
    const auto basis_map = read_basfile(basfile);

    std::vector<gaussian_shell> basis;

    for(const auto & a : atoms)
    {
        const auto & shells = basis_map.at(a.Z);

        // notice the copy below
        for(auto s : shells)
        {
            s.xyz = a.xyz;
            basis.push_back(std::move(s));
        }
    }

    return basis;
}

} // close namespace mirp
