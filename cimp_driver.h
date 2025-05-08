// Ansh Tarafder atarafd1
// Alex Ma ama12
// Amy Wang awang111
#include "ppm_io.h" // PPM I/O header
#include "img_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* printUsage
 * Prints correct command line usage for user to see.
 */
void printUsage(const char *prog_name);

/* runImgProcessing
 * Conducts input validation and calls the respective operations as entered in command line.
 * Returns 0 if successful, otherwise error occurred.
 */
int runImgProcessing(int argc, char **argv);