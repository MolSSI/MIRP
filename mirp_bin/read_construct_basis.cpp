/*! \file
 *
 * \brief Construction of a basis from coordinates and a basis set file
 */


#include "mirp_bin/reffile_io.hpp"
#include "mirp_bin/test_common.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <map>


/*! Factor for converting Angstrom to Bohr
 *
 * Source: NIST 2014 CODATA recommended values
 */
const double ang_to_bohr = 1.0/0.52917721067;


namespace mirp {

// Anonymous namespace for some helper functions and data
namespace {

/*! \brief Maps a character to an integer angular momentum ('s' = 0, etc)
 *
 * Angular momentum characters are stored lower case
 */
std::map<char, int> amchar_map
{
    {'s', 0},  {'p', 1},  {'d', 2},  {'f', 3},  {'g', 4},
    {'h', 5},  {'i', 6},  {'k', 7},  {'l', 8},  {'m', 9},
    {'n', 10}, {'o', 11}, {'q', 12}, {'r', 13}, {'t', 14},
    {'u', 15}, {'v', 16}, {'w', 17}, {'x', 18}, {'y', 19},
    {'z', 20}, {'a', 21}, {'b', 22}, {'c', 23}, {'e', 24}
};


/*! \brief Maps an elemental symbol to its atomic Z number ("he" = 2, etc)
 *
 * Element symbols are stored lower case
 */
std::map<std::string, int> z_map
{
    {"h",   1},   {"he", 2},   {"li",  3},   {"be", 4},   {"b",   5},   {"c",   6},   {"n",  7},
    {"o",   8},   {"f", 9},    {"ne",  10},  {"na", 11},  {"mg",  12},  {"al",  13},  {"si", 14},
    {"p",   15},  {"s", 16},   {"cl",  17},  {"ar", 18},  {"k",   19},  {"ca",  20},  {"sc", 21},
    {"ti",  22},  {"v", 23},   {"cr",  24},  {"mn", 25},  {"fe",  26},  {"co",  27},  {"ni", 28},
    {"cu",  29},  {"zn", 30},  {"ga",  31},  {"ge", 32},  {"as",  33},  {"se",  34},  {"br", 35},
    {"kr",  36},  {"rb", 37},  {"sr",  38},  {"y",  39},  {"zr",  40},  {"nb",  41},  {"mo", 42},
    {"tc",  43},  {"ru", 44},  {"rh",  45},  {"pd", 46},  {"ag",  47},  {"cd",  48},  {"in", 49},
    {"sn",  50},  {"sb", 51},  {"te",  52},  {"i",  53},  {"xe",  54},  {"cs",  55},  {"ba", 56},
    {"la",  57},  {"ce", 58},  {"pr",  59},  {"nd", 60},  {"pm",  61},  {"sm",  62},  {"eu", 63},
    {"gd",  64},  {"tb", 65},  {"dy",  66},  {"ho", 67},  {"er",  68},  {"tm",  69},  {"yb", 70},
    {"lu",  71},  {"hf", 72},  {"ta",  73},  {"w",  74},  {"re",  75},  {"os",  76},  {"ir", 77},
    {"pt",  78},  {"au", 79},  {"hg",  80},  {"tl", 81},  {"pb",  82},  {"bi",  83},  {"po", 84},
    {"at",  85},  {"rn", 86},  {"fr",  87},  {"ra", 88},  {"ac",  89},  {"th",  90},  {"pa", 91},
    {"u",   92},  {"np", 93},  {"pu",  94},  {"am", 95},  {"cm",  96},  {"bk",  97},  {"cf", 98},
    {"es",  99},  {"fm", 100}, {"md",  101}, {"no", 102}, {"lr",  103}, {"rf",  104}, {"db", 105},
    {"sg",  106}, {"bh", 107}, {"hs",  108}, {"mt", 109}, {"ds",  110}, {"rg",  111}, {"cn", 112},
    {"uut", 113}, {"fl", 114}, {"uup", 115}, {"lv", 116}, {"uus", 117}, {"uuo", 118}
};


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
        a.Z = z_map.at(to_lower(element));
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
        am_str = to_lower(am_str);

        int Z = z_map.at(to_lower(element_str)); 

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
                sh2.am = amchar_map.at(am_str.at(i));
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
            sh.am = amchar_map.at(am_str.at(0));
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
