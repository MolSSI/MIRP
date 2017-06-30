#include "mirp_bin/shell.hpp"
#include <stdexcept>

namespace mirp {

int gaussian_shell::get_am(void) const noexcept
{
    return am_;
}

int gaussian_shell::get_nprim(void) const noexcept
{
    return nprim_;
}

int gaussian_shell::get_ngeneral(void) const noexcept
{
    return ngeneral_;
}


int gaussian_shell::get_nfunc(void) const noexcept
{
    int ncart = ((am_+1)*(am_+2))/2;
    return ncart * nprim_ * ngeneral_;
}


gaussian_shell::gaussian_shell(int am, int nprim, int ngeneral)
            : am_(am),
              nprim_(nprim),
              ngeneral_(ngeneral),
              alpha_(nprim),
              coeff_(nprim*ngeneral)
{ }


void
gaussian_shell::set_center(std::string x, std::string y, std::string z)
{
    center_ = std::array<std::string, 3>{{std::move(x),
                                          std::move(y),
                                          std::move(z)}};
}


std::array<std::string, 3>
gaussian_shell::get_center(void) const
{
    return center_;
}


void
gaussian_shell::set_alpha(int prim, std::string alpha)
{
    alpha_.at(prim) = std::move(alpha);
}


std::string
gaussian_shell::get_alpha(int prim) const
{
    return alpha_.at(prim);
}

void
gaussian_shell::set_coeff(int prim, int general, std::string coeff)
{
    coeff_.at(general * nprim_ + prim) = std::move(coeff);
}


std::string
gaussian_shell::get_coeff(int prim, int general) const
{
    return coeff_.at(general * nprim_ + prim);
}

} // close namespace mirp
