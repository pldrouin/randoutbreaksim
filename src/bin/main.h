/**
 * @file main.h
 * @brief Main function for the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include <stdio.h>
#include <unistd.h>

#include <math.h>

#include <gsl/gsl_rng.h>

#include "config.h"
#include "infindividual.h"
#include "simulation.h"
#include "standard_summary_stats.h"

/**
 * @brief Main function.
 */
int main(const int nargs, const char* args[]);
