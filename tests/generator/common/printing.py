
def print_integral_single_input(f, primitives):
    for prim in primitives:
        f.write("{} {} {} {} {} {} {}\n".format(*prim))  
    f.write("\n")


def print_integral_input(f, shells):
    for bf in shells:
        nprim = bf[1]
        ngen = bf[2]
        f.write("{} {} {} {} {} {}\n".format(*bf[:6]))
        for j in range(0, nprim):
            f.write("{}".format(bf[6][j]))
            for k in range(0, ngen):
                f.write(" {}".format(bf[7][k][j]))
            f.write("\n")
    f.write("\n")
