/*! \file
 *
 * \brief Functions for manipulating shell structures
 */

#include "mirp/shell.h"
#include <string.h>
#include <stdlib.h>

int mirp_shell_init(mirp_shell * shell, int am, int n_segment, int n_general)
{
    shell->am = am;
    shell->n_segment = n_segment;
    shell->n_general = n_general;

    shell->alpha = (char **)malloc(n_segment * sizeof(char*));
    if(shell->alpha == NULL)
        return -1;

    shell->coeff = (char **)malloc(n_segment * n_general * sizeof(char*));
    if(shell->coeff == NULL)
        return -1;

    for(int i = 0; i < n_segment; i++)
        shell->alpha[i] = NULL;
    for(int i = 0; i < (n_segment * n_general); i++)
        shell->coeff[i] = NULL;

    shell->center[0] = NULL;
    shell->center[1] = NULL;
    shell->center[2] = NULL;

    return 0;
}

void mirp_shell_free(mirp_shell * shell)
{
    for(int i = 0; i < shell->n_segment; i++)
    {
        if(shell->alpha[i] != NULL)
        {
            free(shell->alpha[i]);
            shell->alpha[i] = NULL;
        }
    }

    for(int i = 0; i < (shell->n_segment * shell->n_general); i++)
    {
        if(shell->coeff[i] != NULL)
        {
            free(shell->coeff[i]);
            shell->coeff[i] = NULL;
        }
    }

    if(shell->center[0] != NULL)
        free(shell->center[0]);
    if(shell->center[1] != NULL)
        free(shell->center[1]);
    if(shell->center[2] != NULL)
        free(shell->center[2]);

    shell->center[0] = NULL;
    shell->center[1] = NULL;
    shell->center[2] = NULL;

    shell->n_segment = 0;
    shell->n_general = 0;
    shell->am = 0;
}

int mirp_shell_set_alpha(mirp_shell * shell, int segment, const char * alpha)
{
    if(segment >= shell->n_segment)
        return -1;

    if(shell->alpha[segment] != NULL)
        free(shell->alpha[segment]);

    int size = strlen(alpha)+1;
    shell->alpha[segment] = (char *)malloc(size * sizeof(char));

    if(shell->alpha[segment] == NULL)
        return -1;

    strncpy(shell->alpha[segment], alpha, size); 
    return 0;
}

int mirp_shell_set_coeff(mirp_shell * shell, int segment, int general, const char * coeff)
{
    if(segment >= shell->n_segment)
        return -1;
    if(general >= shell->n_general)
        return -1;

    int idx = general * shell->n_segment + segment;
    if(shell->coeff[idx] != NULL)
        free(shell->coeff[idx]);

    int size = strlen(coeff)+1;
    shell->coeff[idx] = (char *)malloc(size * sizeof(char));

    if(shell->alpha[segment] == NULL)
        return -1;

    strncpy(shell->coeff[idx], coeff, size); 
    return 0;
}

int mirp_shell_get_alpha(char ** str, mirp_shell * shell, int segment)
{
    if(segment >= shell->n_segment)
        return -1;

    if(shell->alpha[segment] == NULL)
        return -2;

    int size = strlen(shell->alpha[segment])+1;
    (*str) = (char *)malloc(size * sizeof(char));

    if((*str) == NULL)
        return -1;

    strncpy((*str), shell->coeff[segment], size); 
    return 0;
}

int mirp_shell_get_coeff(char ** str, mirp_shell * shell, int segment, int general)
{
    if(segment >= shell->n_segment)
        return -1;
    if(general >= shell->n_general)
        return -1;

    int idx = general * shell->n_segment + segment;
    if(shell->coeff[idx] == NULL)
        return -2;

    int size = strlen(shell->coeff[idx])+1;
    (*str) = (char *)malloc(size * sizeof(char));

    if((*str) == NULL)
        return -1;

    strncpy((*str), shell->coeff[idx], size); 
    return 0;
}

int mirp_iterate_gaussian(int lmn[3])
{
    const int am = lmn[0] + lmn[1] + lmn[2];
    if(lmn[2] >= am)
        return 0;

    if(lmn[2] < (am - lmn[0]))
    {
        lmn[1]--;
        lmn[2]++;
    }
    else
    {
        lmn[0]--;
        lmn[1] = am-lmn[0];
        lmn[2] = 0;
    }
    return 1;
}

int mirp_shell_set_center(mirp_shell * shell, const char * xyz[3])
{
    for(int i = 0; i < 3; i++)
    {
        if(shell->center[i] != NULL)
            free(shell->center[i]);

        int size = strlen(xyz[i])+1;
        shell->center[i] = (char *)malloc(size * sizeof(char));
        if(shell->center[i] == NULL)
            return -1;

        strncpy(shell->center[i], xyz[i], size); 
    }

    return 0;
}

#ifdef __cplusplus
}
#endif

