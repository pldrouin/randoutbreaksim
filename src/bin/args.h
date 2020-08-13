/**
 * @file args.h
 * @brief Generic functions to process configuration parameters.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _ARGS_
#define _ARGS_

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>

/**
 * @brief Reads the next parameter.
 *
 * This function reads the next parameter. 
 *
 * @return the number of read characters.
 */
int getnextparam(FILE **fptra, int *fptri, const bool isarg, const int nargs, char const* const args[], int *parc, char *param);

/**
 * @brief Safely read an expected parameter.
 *
 * This function reads the next parameter. It calls exit(1) if the parameter is
 * missing
 */
void safegetnextparam(FILE **fptra, int *fptri, const bool isarg, const int nargs, char const* const args[], int *parc, char *param);

/**
 * @brief Verifies if two strings differ.
 *
 * This function verifies if two strings differ. It should be faster than strcmp
 * since it does not parse the whole strings unless necessary.
 *
 * @return true if the strings differ, and false otherwise.
 */
inline static bool argsdiffer(const char* arg1, const char* arg2) {int i=0; while(1) {if(arg1[i]!=arg2[i]) return true; if(!arg1[i]) return false; ++i;}}

#endif
