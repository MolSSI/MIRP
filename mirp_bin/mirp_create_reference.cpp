/*! \file
 *
 * \brief mirp_create_reference main function
 */

#include "mirp_bin/cmdline.hpp"
#include "mirp_bin/test_common.hpp"
#include "mirp_bin/ref_integral.hpp"

#include <mirp/kernels/all.h>

#include <sstream>
#include <iostream>
#include <stdexcept>

using namespace mirp;


static void print_help(void)
{
    std::cout << "\n"
              << "mirp_create_reference - Create a reference data file for use with other programs\n"
              << "\n"
              << "\n"
              << "Required arguments:\n"
              << "    --basis        Path to a basis set file (usually ends in .bas)\n"
              << "    --geometry     Path to an XYZ file(usually ends in .xyz)\n"
              << "    --outfile      Output file (usually ends in .dat). Existing data will\n"
              << "                       be overwritten\n"
              << "    --integral     The type of integral to compute. Possibilities are:\n"
              << "                       gtoeri\n"
              << "\n"
              << "\n"
              << "Optional arguments:\n"
              << "    --am           Comma-separated list of AM classes to calculate.\n"
              << "                   The AM should be represented by their letters.\n"
              << "                   (for example, for ERI: --am ssss,psps,dddd)\n"
              << "\n"
              << "\n"
              << "Other arguments:\n"
              << "    -h, --help     Display this help screen\n"
              << "\n";
}



/*! \brief Main function */
int main(int argc, char ** argv)
{
    std::string basfile, xyzfile, outfile;
    std::string integral;
    std::vector<std::vector<int>> amlist;

    try {
        auto cmdline = convert_cmdline(argc, argv);
        if(cmdline.size() == 0 || cmdline_get_switch(cmdline, "-h") || cmdline_get_switch(cmdline, "--help"))
        {
            print_help();
            return 0;
        }

        basfile = cmdline_get_arg_str(cmdline, "--basis");
        xyzfile = cmdline_get_arg_str(cmdline, "--geometry");
        outfile = cmdline_get_arg_str(cmdline, "--outfile");
        integral = cmdline_get_arg_str(cmdline, "--integral");

        if(cmdline_has_arg(cmdline, "--am"))
        {
            std::string amlist_str = cmdline_get_arg_str(cmdline, "--am");
            std::vector<std::string> amlist_split = split(amlist_str, ',');

            for(const auto & s : amlist_split)
            {
                std::vector<int> am_ntet;
                for(char c : s)
                    am_ntet.push_back(amchar_to_int(c));
                amlist.push_back(std::move(am_ntet));
            }

        }

        if(cmdline.size() != 0)
        {
            std::stringstream ss;
            ss << "Unknown command line arguments:\n";
            for(const auto & it : cmdline)
                ss << "  " << it << "\n";
            throw std::runtime_error(ss.str());
        }
    }
    catch(std::exception & ex)
    {
        std::cout << "\nError parsing command line: " << ex.what() << "\n\n";
        std::cout << "Run \"mirp_create_reference -h\" for help\n\n";
        return 1;
    }

    // Create a header from the command line
    std::string header("# Reference values for the ");
    header += integral;
    header += " integral generated with:\n";
    header += "#  ";
    for(int i = 0; i < argc; i++)
        header += " " + std::string(argv[i]);
    header += "\n#\n";

    try
    {
        if(integral == "gtoeri")
        {
            integral4_create_reference(xyzfile, basfile, outfile, header,
                                       amlist, mirp_gtoeri_exact);
        }
        else
        {
            std::cout << "Integral \"" << integral << "\" is not valid\n";
            return 3;
        }
    }
    catch(std::exception & ex)
    {
        std::cout << "Error while running tests: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
