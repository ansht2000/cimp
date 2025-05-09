// Ansh Tarafder atarafd1
// Alex Ma ama12
// Amy Wang awang111
#ifndef IMG_PROCESSING_H
#define IMG_PROCESSING_H

#include "ppm_io.h"
#include <stdlib.h>
#include <math.h>

int grayscale(Image *img);

int binarize(Image *img, int threshold);

int crop(Image *img, int top_left_row, int top_left_col, int bottom_right_row, int bottom_right_col);

int transpose(Image *img);

int computeGrad(Image *img, int cols, int i);

int gradient(Image *img);

void createSeamsAndGradEnergies(Image *grad_img, int **seams, int *grad_energies);

void findLeastEnergySeamAndRemove(Image *img, Image *img_copy, int **seams, int *grad_energies);

int carveSeams(Image *img, int num_seams);

int seam(Image *img, float scale_factor_row, float scale_factor_col);

int max(int left, int right);

int blend(Image *img_one, Image *img_two, Image **img_blend, float alpha);

#endif // IMG_PROCESSING_H
