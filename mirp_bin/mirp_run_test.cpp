/*! \file
 *
 * \brief mirp_run_test main function
 */

#include "mirp_bin/cmdline.hpp"
#include "mirp_bin/test_boys.hpp"
#include "mirp_bin/test_integral4.hpp"

#include <mirp/kernels/all.h>

#include <sstream>
#include <iostream>
#include <stdexcept>

using namespace mirp;

static void print_help(void)
{
    std::cout << "\n"
              << "mirp_run_test - Test MIRP functionality using test data\n"
              << "\n"
              << "\n"
              << "Required arguments:\n"
              << "    --file         File to test with\n"
              << "    --integral     The type of integral to compute. Possibilities are:\n"
              << "                       boys\n"
              << "                       eri\n"
              << "                       eri_single\n"
              << "    --float        Type of floating-point to test with. Possibilities are:\n"
              << "                       interval\n"
              << "                       double\n"
              << "                       exact\n"
              << "    --prec         Working precision in binary digits (bits) to test (required for --float interval)\n"
              << "\n"
              << "\n"
              << "Integral-dependent options:\n"
              << "\n"
              << "  Boys Function:\n"
              << "    --extra-m      Initially compute this many more m values (to test recursion)\n"
              << "\n"
              << "\n"
              << "Other arguments:\n"
              << "    -h, --help     Display this help screen\n"
              << "\n";
}

/*! \brief Main function */
int main(int argc, char ** argv)
{
    std::string file;
    std::string integral;
    std::string floattype;
    long working_prec = 0;
    int extra_m = 0;

    try {
        auto cmdline = convert_cmdline(argc, argv);
        if(cmdline.size() == 0 || cmdline_has_arg(cmdline, "-h") || cmdline_has_arg(cmdline, "--help"))
        {
            print_help();
            return 0;
        }

        file = cmdline_get_arg_str(cmdline, "--file");
        integral = cmdline_get_arg_str(cmdline, "--integral");
        floattype = cmdline_get_arg_str(cmdline, "--float");

        if(floattype != "double" && floattype != "exact")
            working_prec = cmdline_get_arg_long(cmdline, "--prec");
        else if(cmdline_has_arg(cmdline, "--prec"))
            throw std::runtime_error("--prec is not valid for this floating-point type");

        if(integral == "boys")
            extra_m = static_cast<int>(cmdline_get_arg_long(cmdline, "--extra-m", 0));
        else if(cmdline_has_arg(cmdline, "--extra-m"))
            throw std::runtime_error("--extra-m is not valid for this integral type");

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
        std::cout << "Run \"mirp_run_test -h\" for help\n\n";
        return 1;
    }


    try
    {
        long nfailed = -1;
        if(integral == "boys")
        {
            nfailed = boys_run_test_main(file, floattype, extra_m, working_prec);
        }
        else if(integral == "eri_single")
        {
            if(floattype == "interval")
            {
                nfailed = integral4_single_run_test(file, working_prec, mirp_eri_single_str);
            }
            else if(floattype == "double")
            {
                nfailed = integral4_single_run_test_d(file, mirp_eri_single_d);
            }
            else if(floattype == "exact")
            {
                nfailed = integral4_single_run_test_exact(file, mirp_eri_single_exact, mirp_eri_single);
            }
            else
            {
                std::cout << "Float type \"" << floattype << " not valid for integral \"" << integral << "\"\n";
                return 1;
            }
        }
        else if(integral == "eri")
        {
            if(floattype == "interval")
            {
                nfailed = integral4_run_test(file, working_prec, mirp_eri_str);
            }
            else if(floattype == "double")
            {
                nfailed = integral4_run_test_d(file, mirp_eri_d);
            }
            else if(floattype == "exact")
            {
                nfailed = integral4_run_test_exact(file, mirp_eri_exact, mirp_eri);
            }
            else
            {
                std::cout << "Float type \"" << floattype << " not valid for integral \"" << integral << "\"\n";
                return 1;
            }
        }
        else
        {
            std::cout << "Integral \"" << integral << "\" is not valid\n";
            return 3;
        }

        std::cout << nfailed << " tests failed\n";

        if(nfailed)
            return 1;
        else
            return 0;
    }
    catch(std::exception & ex)
    {
        std::cout << "Error while running tests: " << ex.what() << "\n";
        return 1;
    }
}
