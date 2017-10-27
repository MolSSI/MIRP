
def print_integral_single_input(f, primitives):
    for prim in primitives:
        f.write("{} {} {} {} {} {} {} {}\n".format(*prim))  
    f.write("\n")


def print_integral_input(f, shells):
    for bf in shells:
        nprim = bf[2]
        ngen = bf[3]
        f.write("{} {} {} {} {} {} {}\n".format(*bf[:7]))  
        for j in range(0, nprim):
            f.write("{}".format(bf[7][j]))
            for k in range(0, ngen):
                f.write(" {}".format(bf[8][k][j]))
            f.write("\n")
    f.write("\n")
