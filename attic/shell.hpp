#pragma once

#include <string>
#include <vector>
#include <array>

namespace mirp {

class gaussian_shell
{
        int am_;
        int nprim_;
        int ngeneral_;
        std::array<std::string, 3> center_;
        std::vector<std::string> alpha_;
        std::vector<std::string> coeff_;

    public:
        gaussian_shell(int am, int nprim, int ngeneral);

        ~gaussian_shell() = default;
        
        int get_am(void) const noexcept;

        int get_nprim(void) const noexcept;

        int get_ngeneral(void) const noexcept;

        int get_nfunc(void) const noexcept;

        void set_center(std::string x, std::string y, std::string z);

        std::array<std::string, 3> get_center(void) const;

        void set_alpha(int prim, std::string alpha);

        std::string get_alpha(int prim) const;

        void set_coeff(int prim, int general, std::string coeff);

        std::string get_coeff(int prim, int general) const;
};

} // close namespace mirp
