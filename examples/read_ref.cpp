/*! \file
 *
 * \brief Example of reading a reference file for a 4-center integral
 */

#include <vector>
#include <array>
#include <iostream>
#include <locale> // tolower, isspace
#include <limits>
#include <fstream>

///////////////////////////////////////////
// This example is completely standalone
// and does not depend on MIRP in any way
///////////////////////////////////////////
    
// Maximum length of a line
static const std::streamsize max_length = std::numeric_limits<std::streamsize>::max();

/* Number of cartesian functions for a given angular momentum */
#define NCART(am) ((((am)+1)*((am)+2))/2)

/* Number of cartesian functions for 2 shells */
#define NCART2(am1, am2) \
        (NCART((am1)) * NCART((am2)))

/* Number of cartesian functions for 4 shells */
#define NCART4(am1, am2, am3, am4) \
        (NCART2((am1),(am2)) * NCART2((am3),(am4)))


/* A shell of cartesian gaussian functions (double precision) 
 *
 * This shell holds information as double precision
 */
struct gaussian_shell
{
    int Z;                      // Z-number of the center
    int am;                     // Angular momentum
    int nprim;                  // Number of primitives (segmented contraction)
    int ngeneral;               // Number of general contractions
    std::array<double, 3> xyz;  // Coordinates (in bohr)
    std::vector<double> alpha;  // Exponents of the gaussians
    std::vector<double> coeff;  // Contraction coefficients (unnormalized)
};



/* This reads a single double-precision number in hexfloat
 * format. This uses strtod rather than normal stream extraction
 * since some standard libraries don't support reading hexfloats
 * like that
 *
 * The stream is advanced past the hexfloat
 */
double read_hexdouble(std::istream & fs)
{

    std::string tmp;
    fs >> tmp;
    return std::strtod(tmp.c_str(), nullptr);
}


/* Skips blank lines and lines beginning with a comment character */
bool file_skip(std::istream & fs, char commentchar)
{
    while(true)
    {
        int c = fs.peek();

        if(!fs.good())                    // did we reach eof?
            break;
        else if(c == commentchar)
            fs.ignore(max_length, '\n');  // skip the rest of the line
        else if(std::isspace(c))          // is whitespace? (includes \n, \t, etc)
            fs.get();                     // discard it
        else
            break;                        // We found some data
    }

    return fs.good();
}


/* Reads in the basis information from the top of the file */
std::vector<gaussian_shell> reffile_read_basis(std::istream & fs)
{
    size_t nshell;

    file_skip(fs, '#');
    fs >> nshell;

    std::vector<gaussian_shell> shells(nshell);

    for(size_t i = 0; i < nshell; i++)
    {
        auto & s = shells[i];

        // angular momentum, # of primitives, # of general contractions
        file_skip(fs, '#');
        fs >> s.am >> s.nprim >> s.ngeneral;

        // coordinates
        file_skip(fs, '#');
        s.xyz[0] = read_hexdouble(fs);
        s.xyz[1] = read_hexdouble(fs);
        s.xyz[2] = read_hexdouble(fs);

        s.alpha.resize(s.nprim);
        s.coeff.resize(s.nprim*s.ngeneral);

        file_skip(fs, '#');
        for(int j = 0; j < s.nprim; j++)
            s.alpha[j] = read_hexdouble(fs);

        for(int n = 0; n < s.ngeneral; n++)
        {
            file_skip(fs, '#');
            for(int j = 0; j < s.nprim; j++)
                s.coeff[n*s.nprim+j] = read_hexdouble(fs);
        }
    }

    return shells;
}


int main(int argc, char ** argv)
{
    if(argc != 2)
    {
        std::cout << "Must give reference file as the only argument\n";
        return 1;
    }


    std::cout << "Opening file " << argv[1] << "\n";
    std::ifstream fs(argv[1]);
    if(!fs.is_open())
        throw std::runtime_error("Error opening input file");

    // read in the basis set information
    file_skip(fs, '#');
    std::vector<gaussian_shell> shells = reffile_read_basis(fs);

    std::cout << "Read " << shells.size() << " shells\n";
    size_t nintegrals_total = 0;

    while(fs.good())
    {
        // Shell quartet indices
        size_t p, q, r, s;
        fs >> p >> q >> r >> s;

        const gaussian_shell & s1 = shells.at(p);
        const gaussian_shell & s2 = shells.at(q);
        const gaussian_shell & s3 = shells.at(r);
        const gaussian_shell & s4 = shells.at(s);

        const size_t ncart = NCART4(s1.am, s2.am, s3.am, s4.am);
        const size_t ngen = s1.ngeneral * s2.ngeneral * s3.ngeneral * s4.ngeneral;
        const size_t nintegrals = ncart * ngen;

        std::vector<double> integrals_file(nintegrals);
        for(size_t i = 0; i < nintegrals; i++)
            integrals_file[i] = read_hexdouble(fs);

        // print what we read and the first/last integrals
        std::cout << "\nQuartet " << p << " " << q << " " << r << " " << s << " has " << nintegrals << " integrals\n";
        std::cout << "    " << integrals_file[0] << "\n";
        if(nintegrals > 1)
        {
            std::cout << "    ....\n";
            std::cout << "    " << integrals_file[nintegrals-1] << "\n";
        }
    
        nintegrals_total += nintegrals;

        if(!file_skip(fs, '#'))
            break;
    }

    std::cout << "\nRead in " << nintegrals_total << " integrals from the file\n";
    return 0;
}

