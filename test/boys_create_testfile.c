#include "mirp/boys.h"

#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <string.h> /* for strerror */
#include <errno.h>

#define MAXLINE 2048

static void print_usage(void)
{
    printf("Usage: boys_create_testfile [-h] [-p precision] [-n max_n] -i infile -o outfile outfile\n");
    printf("\n");
    printf("    -h             Print this help text\n");
    printf("    -p precision   Set the precision used in the calculation (default is 256 bits)\n");
    printf("    -n max_n       Maximum order of the Boys function to calculate\n");
    printf("    -i infile      File from which to take the x values for the Boys function\n");
    printf("    -o outfile     File to which to write the output\n");
    printf("\n");
}

static void print_line(int n, arb_t F, arb_t x, int ndigits)
{
    char * s_F = arb_get_str(F, ndigits, ARB_STR_NO_RADIUS); 
    char * s_x = arb_get_str(x, ndigits, ARB_STR_NO_RADIUS);

    printf("%6d   %s   %s\n", n, s_x, s_F);

    free(s_F);
    free(s_x);
}


int main(int argc, char **argv)
{
    /* Defaults */
    int prec = 256;
    int maxn = 50;
    int ndigits = 64;
    char * infile = NULL;
    char * outfile = NULL;

    /* Parse the arguments */
    int opt;
    while( (opt = getopt(argc, argv, "hp:n:i:o:")) != -1 )
    {
        switch(opt)
        {
            case 'h':
                print_usage();
                return 0;
            case 'p':
                prec = atoi(optarg);
                break;
            case 'n':
                maxn = atoi(optarg);
                break;
            case 'i':
                infile = optarg;
                break;
            case 'o':
                outfile = optarg;
                break;
            case '?':
                printf("Unknown option '%s'\n", optarg);
                return 1; 
        }
    }

    /* Check to see the correct options were set */
    if(infile == NULL)
    {
        printf("Input file (-i) not specified\n");
        return 2;
    }

    if(outfile == NULL)
    {
        printf("Output file (-o) not specified\n");
        return 2;
    }

    if(maxn < 0)
    {
        printf("Maximum order (-n) must be >= 0\n");
        return 2;
    }

    if(ndigits <= 0)
    {
        printf("Number of digits to print (-d) must be > 0\n");
        return 2;
    }

    if(prec <= 0)
    {
        printf("Number of bits of precision (-p) must be > 0\n");
        return 2;
    }

    /* Open the input and output files */
    FILE *fin = fopen(infile, "r");
    if(fin == NULL)
    {
        printf("Unable to open input file \"%s\" - %s\n", infile, strerror(errno));
        return 3;
    }

    FILE *fout = fopen(outfile, "w");
    if(fout == NULL)
    {
        printf("Unable to open output file \"%s\" - %s\n", outfile, strerror(errno));
        fclose(fin);
        return 3;
    }

    printf("-----------------------------------------------\n");
    printf("Summary\n");
    printf("-----------------------------------------------\n");
    printf("   Input file: %s\n", infile);
    printf("  Output file: %s\n", outfile);
    printf("    Precision: %d\n", prec);
    printf("Output digits: %d\n", ndigits);
    printf("    Max order: %d\n", maxn);
    printf("-----------------------------------------------\n");


    /* Read in the x values */
    ssize_t nread;
    char * line = NULL;
    size_t len = 0;
    while((nread = getline(&line, &len, fin)) != -1)
    {
        printf("Read line: %s\n", line); 
    }

    if(line)
        free(line);
    
/*
    double x,      F[maxn+1];
    arb_t  arb_x,  arb_F[maxn+1];
    mpfr_t mpfr_x, mpfr_F[maxn+1];

    arb_init(arb_x);
    mpfr_init2(mpfr_x, prec);

    for(int i = 0; i <= maxn; i++)
    {
        arb_init(arb_F[i]);
        mpfr_init2(mpfr_F[i], prec);
    }

    for(unsigned long i = 0; i < 5; i++)
    {
        x = (double)i / 10.0;

        arb_set_ui(arb_x, i);
        arb_div_ui(arb_x, arb_x, 10, prec);

        mpfr_set_ui(mpfr_x, i, MPFR_RNDN);
        mpfr_div_ui(mpfr_x, mpfr_x, 10, MPFR_RNDN);

        mirp_boys(F, maxn, x);
        mirp_boys_mp(prec, mpfr_F, maxn, mpfr_x);
        mirp_boys_interval(prec, arb_F, maxn, arb_x);

        for(int n = 0; n <= maxn; n++)
            print_line(n, arb_F[n], arb_x);
    }
*/
/*
    for(int i = 0; i < n; i++)
    {
        printf("-------------------------------------\n");
        printf("n = %d\n", i);
        printf(     "           double precision:  %-26.16e\n", F[i]); 
        mpfr_printf("    multi-precision [%5d]:  %Re\n", prec, mpfr_F[i]);
        printf(     "interval arithmetic [%5d]: %s\n", prec, arb_get_str(arb_F[i], 32, 0));
        printf("        bits of accuracy: %ld\n", arb_rel_accuracy_bits(arb_F[i]));
        printf("-------------------------------------\n");
    }
*/  
/*
    arb_clear(arb_x);    
    mpfr_clear(mpfr_x);

    for(int i = 0; i <= maxn; i++)
    {
        arb_clear(arb_F[i]);
        mpfr_clear(mpfr_F[i]);
    }
*/ 

    fclose(fin);
    fclose(fout);
    return 0;
}
