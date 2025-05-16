// Ansh Tarafder atarafd1
// Alex Ma ama12
// Amy Wang awang111
#ifndef IMG_PROCESSING_H
#define IMG_PROCESSING_H

#include "ppm_io.h"
#include <stdlib.h>

int grayscale(Image *img);

int binarize(Image *img, int threshold);

int crop(Image *img, int top_left_row, int top_left_col, int bottom_right_row, int bottom_right_col);

int flip(Image *img, unsigned char axis);

int rotate(Image *img, unsigned char direction);

int gradient(Image *img);

int seam(Image *img, float scale_factor_row, float scale_factor_col);

int blend(Image *img_one, Image *img_two, Image **img_blend, float alpha);

int pointilism(Image *img);

int dither(Image *img);

#endif // IMG_PROCESSING_H
