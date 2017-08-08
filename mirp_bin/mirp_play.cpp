/*! \file
 *
 * \brief Sandbox for testing features
 */

#include "mirp_bin/cmdline.hpp"
#include <mirp/mirp.h>
#include <iostream>

using namespace mirp;

int main(void)
{
    double alpha1[3] = { 130.7093200, 23.8088610, 6.4436083 };
    double coeff1[3] = { 0.15432897, 0.53532814, 0.44463454 };
    double A[3] = { 0.0, 0.0, 0.0 };
    double integral;

    mirp_integral4_exact(&integral,
                         0, A, 3, 1, alpha1, coeff1,
                         0, A, 3, 1, alpha1, coeff1,
                         0, A, 3, 1, alpha1, coeff1,
                         0, A, 3, 1, alpha1, coeff1,
                         mirp_eri_interval);

    printf("Value: %32.18e\n", integral);

    return 0;
}
