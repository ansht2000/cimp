#include "img_processing.h"
#include "ppm_io.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

static int clamp(int val, int min, int max) {
    if (val > max) {return max;}
    if (val < min) {return min;}
    return val;
}

/* grayscale
 * Convert the input image to grayscale.
 * Return 0 if successful.
 */
int grayscale(Image *img) {
    const size_t img_data_length = (size_t) img->rows * img->cols;
    uint8_t gray;
    Pixel *data = img->data;
    // assign the gray value for each pixel
    for (size_t i = 0; i < img_data_length; ++i, ++data) {
        gray = (77 * data->r + 150 * data->g + 29 * data->b) >> 8;
        data->r = data->g = data->b = gray;
    }
    return 0;
}

/* binarize
 * Convert the input image to black and white by thresholding.
 * Return 0 if successful.
 */
int binarize(Image *img, int threshold) {
    // grayscale image first
    grayscale(img);
    const size_t img_data_length = img->rows * img->cols;
    uint8_t gray;
    Pixel *data = img->data;
    // checks each pixel and converts its color to white/black depending on relation to threshold
    for (size_t i = 0; i < img_data_length; ++i, ++data) {
        // could pick any of the three channels
        gray = data->r;
        if (gray < threshold) {
            // set pixel color to black
            data->r = data->g = data->b = 0;
        } else {
            // set pixel color to white
            data->r = data->g = data->b = 255;
        }
    }
    return 0;
}

/* crop
 * Crop the input image given corner pixel locations.
 * Return 0 if successful.
 */
int crop(Image *img,
         int top_left_col, int top_left_row,
         int bottom_right_col, int bottom_right_row) {
    // set cropped image height and width
    size_t new_rows = bottom_right_row - top_left_row;
    size_t new_cols = bottom_right_col - top_left_col;
    Pixel *data = malloc(sizeof(Pixel) * new_rows * new_cols);
    // initialize cropped image data
    if (data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 8;
    }
    // iterate through the rows of the new data
    // and copy over the old data in one shot
    for (size_t i = 0; i < new_rows; ++i) {
        memcpy(data + new_cols * i,
               img->data + (top_left_row + i) * img->cols + top_left_col,
               sizeof(Pixel) * new_cols);
    }
    free(img->data);
    // update image with new rows, cols, data
    img->rows = new_rows;
    img->cols = new_cols;
    img->data = data;
    return 0;
}

unsigned char toLowerChar(unsigned char l) {
    if (l < 'a') {
        l ^= 32;
    }
    return l;
}

/* flip
 * Flip the image across the specified axis, 'v' for vertical axis, 'h' for horizontal axis
 * Return 0 if successful
 */
int flip(Image *img, unsigned char axis) {
    // if axis is uppercase, toggle it to lowercase
    axis = toLowerChar(axis);
    Pixel *data = malloc(img->cols * img->rows * sizeof(Pixel));
    if (data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 8;
    }
    switch (axis) {
    case 'v':
        for (size_t i = 0; i < (size_t)img->cols; ++i) {
            for (size_t j = 0; j < (size_t)img->rows; ++j) {
                data[j * img->cols + i] = img->data[j * img->cols + (img->cols - (i + 1))];
            }
        }
        break;
    case 'h':
        for (size_t i = 0; i < (size_t)img->rows; ++i) {
            memcpy(data + img->cols * i,
                   img->data + img->cols * (img->rows - (i + 1)),
                   sizeof(Pixel) * img->cols);
        }
        break;
    default:
        free(data);
        return 1;
    }
    free(img->data);
    img->data = data;
    return 0;
}

/* rotate
 * Rotate the image counter-clockwise or clockwise, 'l' for counter-clockwise, 'r' for clockwise
 * Return 0 if successful
 */
int rotate(Image *img, unsigned char direction) {
    direction = toLowerChar(direction);
    Pixel *data = malloc(img->cols * img->rows * sizeof(Pixel));
    if (data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 8;
    }
    size_t new_rows = img->cols;
    size_t new_cols = img->rows;
    switch (direction) {
    case 'l':
        // to rotate left, start copying og image data from the top right
        // and move pixel pointer columnwise
        for (size_t i = 0; i < new_rows; ++i) {
            for (size_t j = 0; j < new_cols; ++j) {
                data[i * new_cols + j] = img->data[j * img->cols + (img->cols - (i + 1))];
            }
        }
        break;
    case 'r':
        // to rotate right, start copying og image data from the bottom left
        // and move pixel pointer columnwise
        for (size_t i = 0; i < new_rows; ++i) {
            for (size_t j = 0; j < new_cols; ++j) {
                data[i * new_cols + j] = img->data[(img->rows - (j + 1)) * img->cols + i];
            }
        }
        break;
    default:
        free(data);
        return 1;
    }
    free(img->data);
    img->data = data;
    img->rows = new_rows;
    img->cols = new_cols;
    return 0;
}

/* transpose
 * Transpose the input image.
 * Return 0 if successful.
 */
int transpose(Image *img) {
    // flip image dimensions
    int trans_cols = img->rows;
    int trans_rows = img->cols;
    // allocate space for transpose image data
    Pixel *data = malloc(sizeof(Pixel) * (trans_rows) * trans_cols);
    if (data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 8;
    }

    // loop through rows and columns of transpose image data
    for (int i = 0; i < trans_rows; i++) {
        for (int j = 0; j < trans_cols; j++) {
            data[trans_cols * i + j] = img->data[trans_rows * j + i];
        }
    }

    free(img->data);
    // update image with new rows, cols, data
    img->rows = trans_rows;
    img->cols = trans_cols;
    img->data = data;
    return 0;
}

/* computeGrad
 * Helper function to compute the gradient.
 * Return the gradient.
 */
static int computeGrad(Image *img, size_t cols, size_t i) {
    // compute gradient from neighboring pixel values
    // by convolving with sobel 3 x 3 matrix:
    // g_x = [-1  0  1]   g_y = [ 1  2  1]
    //       [-2  0  2]       = [ 0  0  0]
    //       [-1  0  1]       = [-1 -2  1]
    // each directional gradient is scaled by 1/8 to get pixel accurate value
    int grad_x = (-img->data[i - (cols + 1)].r - 2 * img->data[i - 1].r - img->data[i + cols - 1].r
                 + img->data[i - (cols - 1)].r + 2 * img->data[i + 1].r + img->data[i + cols + 1].r) >> 3;
    int grad_y = (img->data[i - (cols + 1)].r + 2 * img->data[i - cols].r + img->data[i - (cols - 1)].r
                - img->data[i + cols - 1].r - 2 * img->data[i + cols].r - img->data[i + cols + 1].r) >> 3;
    // get norm of gradient
    uint8_t grad = sqrt(grad_x * grad_x + grad_y * grad_y);
    return grad;
}

/* gradient
 * Compute the gradient of an input image and save its magnitude as an image.
 * Return 0 if successful.
 */
int gradient(Image *img) {
    grayscale(img);
    const size_t cols = img->cols;
    const size_t rows = img->rows;
    const size_t img_data_length = rows * cols;
    uint8_t grad;
    // have to create new data because working on the same array
    // will mess up the gradient calculation
    Pixel *data = (Pixel *)malloc(sizeof(Pixel) * cols * rows);
    if (data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 8;
    }
    // set max_grad to 0 because all the grad values are positive
    uint8_t max_grad = 0;
    // iterate through image array and assign pixel values
    for (size_t i = 0; i < img_data_length; ++i, ++data) {
        if (i < cols || i > cols * (rows - 1) || i%cols == 0 || i%cols == cols - 1) {
            // set pixel on edges to black
            data->r = data->g = data->b = 0;
        } else {
            // set pixel to gradient value
            grad = computeGrad(img, cols, i);
            max_grad = imax(max_grad, grad);
            data->r = data->g = data->b = grad;
        }
    }
    // find the value to scale each pixel by so the max is 255
    float scale = 255.0 / max_grad;
    // reset data pointer to point to beginning of data
    data -= img_data_length;
    for (size_t i = 0; i < img_data_length; ++i, ++data) {
        data->r = data->g = data->b = clamp(round(scale * data->r), 0, 255);
    }
    free(img->data);
    img->data = data - img_data_length;
    return 0;
}

/* blend
 * Blend two given images together.
 * Return 0 if successful
*/
int blend(Image *img_one, Image *img_two, Image **img_blend, float alpha) {
    if (alpha < 0 || alpha > 1) {
        fprintf(stderr, "alpha value must be in the range [0, 1]\n");
        return 1;
    }

    int blend_img_cols = imax(img_one->cols, img_two->cols);
    int blend_img_rows = imax(img_one->rows, img_two->rows);
    *img_blend = NewImage(blend_img_rows, blend_img_cols);

    for (int i = 0; i < blend_img_rows; i++) {
        for (int j = 0; j < blend_img_cols; j++) {
            if (i > img_one->rows || j > img_one->cols) { 
                (*img_blend)->data[i * blend_img_cols + j] = img_two->data[i * img_two->cols + j];
            } else if (i > img_two->rows || j > img_two->cols) {
                (*img_blend)->data[i * blend_img_cols + j] = img_one->data[i * img_one->cols + j];
            } else {
                float pix_one_contrib_r = img_one->data[i * img_one->cols + j].r * alpha;
                float pix_one_contrib_g = img_one->data[i * img_one->cols + j].g * alpha;
                float pix_one_contrib_b = img_one->data[i * img_one->cols + j].b * alpha;
                float pix_two_contrib_r = img_two->data[i * img_two->cols + j].r * (1 - alpha);
                float pix_two_contrib_g = img_two->data[i * img_two->cols + j].g * (1 - alpha);
                float pix_two_contrib_b = img_two->data[i * img_two->cols + j].b * (1 - alpha);
                (*img_blend)->data[i * blend_img_cols + j].r = pix_one_contrib_r + pix_two_contrib_r;
                (*img_blend)->data[i * blend_img_cols + j].g = pix_one_contrib_g + pix_two_contrib_g;
                (*img_blend)->data[i * blend_img_cols + j].b = pix_one_contrib_b + pix_two_contrib_b;
            }
        }
    }

    return 0;
}

static int randInRange(int min, int max) {
    int rand_num = rand() % (max - min + 1) + min;
    return rand_num;
}

static void drawCircle(Point top_left, Point bottom_right, Point center, int radius, Image *og_img, Pixel *img_point_data) {
    int cols_len = og_img->cols;
    for (int i = top_left.x; i < bottom_right.x; i++) {
        for (int j = top_left.y; j < bottom_right.y; j++) {
            int x_len = abs(i - center.x);
            int y_len = abs(j - center.y);
            if (round(sqrt(x_len * x_len + y_len * y_len)) < radius) {
                img_point_data[j * cols_len + i] = og_img->data[center.y * cols_len + center.x];
            }
        }
    }
}

int pointilism(Image *img) {
    Pixel *img_point_data = malloc(img->cols * img->rows * sizeof(Pixel));
    // TODO: add support for custom background color
    // right now default to white for a canvas aesthetic
    for (int i = 0; i < img->cols * img->rows; i++) {
        img_point_data[i].r = 255;
        img_point_data[i].g = 255;
        img_point_data[i].b = 255;
    }
    int radius;
    for (int i = 0; i < img->rows; i++) {
        for (int j = 0; j < img->cols; j++) {
            int chance = randInRange(1, 100);
            if (chance == 69 || chance == 42 || chance == 28) {
                radius = randInRange(3, 10);
                Point top_left = {imax(j - radius, 0), imax(i - radius, 0)};
                struct _point bottom_right = {imin(j + radius, img->cols), imin(i + radius, img->rows)};
                struct _point center = {j, i};
                drawCircle(top_left, bottom_right, center, radius, img, img_point_data);
            }
        }
    }
    
    free(img->data);
    img->data = img_point_data;

    return 0;
}

static int findClosestPaletteColor(int old_pix_col) {
    return (old_pix_col < 128) ? 0 : 255;
}

int dither(Image *img) {
    grayscale(img);
    for (int i = 0; i < img->rows; i++) {
        for (int j = 0; j < img->cols; j++) {
            int idx = i * img->cols + j;
            int old_pix_col = img->data[idx].r;
            int new_pix_col = findClosestPaletteColor(old_pix_col);
            img->data[idx].r = new_pix_col;
            img->data[idx].g = new_pix_col;
            img->data[idx].b = new_pix_col;
            int quant_error = old_pix_col - new_pix_col;
            if (j + 1 < img->cols) {
                int idx = i * img->cols + (j + 1);
                img->data[idx].r = clamp(img->data[idx].r + round(quant_error * (7.0/16.0)), 0, 255);
                img->data[idx].g = clamp(img->data[idx].g + round(quant_error * (7.0/16.0)), 0, 255);
                img->data[idx].b = clamp(img->data[idx].b + round(quant_error * (7.0/16.0)), 0, 255);
            }
            if (i + 1 < img->rows) {
                int idx = (i + 1) * img->cols + j;
                img->data[idx].r = clamp(img->data[idx].r + round(quant_error * (5.0/16.0)), 0, 255);
                img->data[idx].g = clamp(img->data[idx].g + round(quant_error * (5.0/16.0)), 0, 255);
                img->data[idx].b = clamp(img->data[idx].b + round(quant_error * (5.0/16.0)), 0, 255);
                if (j - 1 >= 0) {
                    int idx = (i + 1) * img->cols + (j - 1);
                    img->data[idx].r = clamp(img->data[idx].r + round(quant_error * (3.0/16.0)), 0, 255);
                    img->data[idx].g = clamp(img->data[idx].g + round(quant_error * (3.0/16.0)), 0, 255);
                    img->data[idx].b = clamp(img->data[idx].b + round(quant_error * (3.0/16.0)), 0, 255);
                }
                if (j + 1 < img->cols) {
                    int idx = (i + 1) * img->cols + (j + 1);
                    img->data[idx].r = clamp(img->data[idx].r + round(quant_error * (1.0/16.0)), 0, 255);
                    img->data[idx].g = clamp(img->data[idx].g + round(quant_error * (1.0/16.0)), 0, 255);
                    img->data[idx].b = clamp(img->data[idx].b + round(quant_error * (1.0/16.0)), 0, 255);
                }
            }
        }
    }

    return 0;
}