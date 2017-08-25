/*! \file
 *
 * \brief mirp_run_test main function
 */

#include "mirp_bin/cmdline.hpp"
#include "mirp_bin/boys_test.hpp"
#include "mirp_bin/test_common.hpp"
#include <mirp/kernels/all.h>
#include <iostream>

using namespace mirp;

/*! \brief Main function */
int main(int argc, char ** argv)
{
    std::string file;
    std::string integral;
    std::string floattype;
    long target_prec = 0;
    long extra_m = 0;

    try {
        auto cmdline = convert_cmdline(argc, argv);

        file = cmdline_get_arg_str(cmdline, "--file");
        integral = cmdline_get_arg_str(cmdline, "--integral");
        floattype = cmdline_get_arg_str(cmdline, "--float");

        if(floattype != "double" && floattype != "exact")
            target_prec = cmdline_get_arg_long(cmdline, "--prec");

        if(integral == "boys")
            extra_m = cmdline_get_arg_long(cmdline, "--extra-m", 0);

    }
    catch(std::exception & ex)
    {
        std::cout << "Error parsing command line: " << ex.what() << "\n";
        return 1;
    }


    try
    {
        long nfailed = -1;
        if(integral == "boys")
        {
            nfailed = boys_run_test_main(file, floattype, extra_m, target_prec);
        }
        else if(integral == "eri_single")
        {
            if(floattype == "interval")
            {
                nfailed = detail::integral4_single_run_test(file, target_prec, mirp_eri_single_target_str);
            }
            else if(floattype == "double")
            {
                nfailed = detail::integral4_single_run_test_d(file, mirp_eri_single_d);
            }
            else if(floattype == "exact")
            {
                nfailed = detail::integral4_single_run_test_exact(file, mirp_eri_single_exact, mirp_eri_single_target);
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
                nfailed = detail::integral4_run_test(file, target_prec, mirp_eri_target_str);
            }
            else if(floattype == "double")
            {
                nfailed = detail::integral4_run_test_d(file, mirp_eri_d);
            }
            else if(floattype == "exact")
            {
                nfailed = detail::integral4_run_test_exact(file, mirp_eri_exact, mirp_eri_target);
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

    return 0;
}
