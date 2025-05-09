// Ansh Tarafder atarafd1
// Alex Ma ama12
// Amy Wang awang111
#include "img_processing.h"
#include "ppm_io.h"
#include <math.h>
#include <stdlib.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/* grayscale
 * Convert the input image to grayscale.
 * Return 0 if successful.
 */
int grayscale(Image *img) {
    int img_data_length = img->rows * img->cols;
    unsigned char red;
    unsigned char green;
    unsigned char blue;
    double grayscale;
    Pixel buffer_pixel;
    // assign the grayscale value for each pixel
    for (int i = 0; i < img_data_length; i++) {
        buffer_pixel = img->data[i];
        red = buffer_pixel.r;
        green = buffer_pixel.g;
        blue = buffer_pixel.b;
        grayscale = (double) (.3 * red + .59 * green + .11 * blue);
        img->data[i].r = (unsigned char) floor(grayscale);
        img->data[i].g = (unsigned char) floor(grayscale);
        img->data[i].b = (unsigned char) floor(grayscale);
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
    int img_data_length = img->rows * img->cols;
    Pixel buffer_pixel;
    int grayscale;
    // checks each pixel and converts its color to white/black depending on relation to threshold
    for (int i = 0; i < img_data_length; i++) {
        buffer_pixel = img->data[i];
        grayscale = (int) buffer_pixel.r;
        if (grayscale < threshold) {
            // set pixel color to black
            img->data[i].r = 0;
            img->data[i].g = 0;
            img->data[i].b = 0;
        } else {
            // set pixel color to white
            img->data[i].r = 255;
            img->data[i].g = 255;
            img->data[i].b = 255;
        }
    }
    return 0;
}

/* crop
 * Crop the input image given corner pixel locations.
 * Return 0 if successful.
 */
int crop(Image *img, int top_left_col, int top_left_row, int bottom_right_col, int bottom_right_row) {
    // set cropped image height and width
    int im_crop_height = bottom_right_row - top_left_row;
    int im_crop_width = bottom_right_col - top_left_col;
    Image im_crop;
    im_crop.data = malloc(sizeof(Pixel) * im_crop_height * im_crop_width);
    // initialize cropped image data
    if (im_crop.data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 8;
    }
    // iterate through cropped image data
    for (int i = 0; i < im_crop_height; i++) {
        for (int j = 0; j < im_crop_width; j++) {
            int crop_index = im_crop_width * i + j;
            int full_index = img->cols * (top_left_row + i) + (top_left_col + j);
            im_crop.data[crop_index] = img->data[full_index];
        }
    }
    free(img->data);
    // update image with new rows, cols, data
    img->rows = im_crop_height;
    img->cols = im_crop_width;
    img->data = im_crop.data;
    return 0;
}

/* transpose
 * Transpose the input image.
 * Return 0 if successful.
 */
int transpose(Image *img) {
    Image im_trans;
    // flip image dimensions
    int trans_cols = img->rows;
    int trans_rows = img->cols;
    // initialize transpose image dimensions to use for initialization later
    im_trans.rows = trans_rows;
    im_trans.cols = trans_cols;
    // allocate space for transpose image data
    im_trans.data = malloc(sizeof(Pixel) * (trans_rows) * trans_cols);
    if (im_trans.data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 8;
    }
    // loop through rows and columns of transpose image data
    for (int i = 0; i < trans_rows; i++) {
        for (int j = 0; j < trans_cols; j++) {
            int trans_index = trans_cols * i + j;
            int mat_index = trans_rows * j + i;
            im_trans.data[trans_index] = img->data[mat_index];
        }
    }
    free(img->data);
    // update image with new rows, cols, data
    img->rows = trans_rows;
    img->cols = trans_cols;
    img->data = im_trans.data;
    return 0;
}

/* computeGrad
 * Helper function to compute the gradient.
 * Return the gradient.
 */
int computeGrad(Image *img, int cols, int i) {
    // get rgb values for neighboring pixels
    unsigned char gray_right = img->data[i + 1].r;
    unsigned char gray_left = img->data[i - 1].r;
    unsigned char gray_down = img->data[i + cols].r;
    unsigned char gray_up = img->data[i - cols].r;
    // compute gradient from neighboring pixel values
    int grad_x = (gray_right - gray_left) / 2;
    int grad_y = (gray_down - gray_up) / 2;
    int grad = floor(fabs((float)grad_x) + fabs((float)grad_y));
    if (grad > 255) grad = 255;
    return grad;
}

/* gradient
 * Compute the gradient of an input image and save its magnitude as an image.
 * Return 0 if successful.
 */
int gradient(Image *img) {
    grayscale(img);
    int cols = img->cols;
    int rows = img->rows;
    int grad;
    // have to create new data because working on the same array
    // will mess up the gradient calculation
    Pixel *grad_data = (Pixel *)malloc(sizeof(Pixel) * cols * rows);
    if (grad_data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 8;
    }
    int img_length = rows * cols;
    // iterate through image array and assign pixel values
    for (int i = 0; i < img_length; i++) {
        if (i < cols || i > cols * (rows - 1) || i%cols == 0 || i%cols == cols - 1) {
            // set pixel to black
            grad_data[i].r = 0;
            grad_data[i].g = 0;
            grad_data[i].b = 0;
        } else {
            // set pixel to grad value calculate by computeGrad
            grad = computeGrad(img, cols, i);
            grad_data[i].r = (unsigned char)grad;
            grad_data[i].g = (unsigned char)grad;
            grad_data[i].b = (unsigned char)grad;
        }
    }
    free(img->data);
    img->data = grad_data;
    return 0;
}

/* createSeamsAndGradEnergies
 * Create seams and find the gradient energy for each seam
 */
void createSeamsAndGradEnergies(Image *grad_img, int **seams, int *grad_energies) {
    // initialize variables for all the image dimensions
    int grad_img_rows = grad_img->rows;
    int grad_img_cols = grad_img->cols;
    // loop through gradient image in column-major order to create seams
    // ignores first and last columns because the gradient will be identically 0
    for (int i = 0; i < grad_img_cols; i++) {
        for (int j = 1; j < grad_img_rows; j++) {
            int next_col_index = seams[i][j - 1];
            // if at the first column of the image, set next_col_index to the second column
            if (next_col_index == 0) {next_col_index = 1;}
            // if at the last column of the image, set next_col_index to the second to last column
            if (next_col_index == grad_img_cols - 1) {next_col_index = grad_img_cols - 2;}
            int grad_value = (int) grad_img->data[j * grad_img_cols + next_col_index].r;
            // loop through pixel's bottom three neighbors
            for (int k = -1; k <= 1; ++k) {
                // set potential least neighbor column index and overall index in image data
                int potential_col_index = seams[i][j - 1] + k;
                int potential_index = j * grad_img_cols + potential_col_index;
                // if potential index is within the image bounds and grad value is less than neighbors
                // set column index to that pixel's column and set grad data to that pixel's gray value
                if (potential_col_index > 0 && potential_col_index < grad_img_cols - 1 && (int) grad_img->data[potential_index].r < grad_value) {
                    grad_value = (int) grad_img->data[potential_index].r;
                    next_col_index = potential_col_index;
                }
            }
            // store column index in seams and updated the gradient energy for that seam
            seams[i][j] = next_col_index;
            grad_energies[i] += grad_value;
        }
    }
    // set columns for the seams in the last two rows to the immediate bottom neighbor
    for (int i = 0, j = grad_img_rows - 1; i < grad_img_cols; i++) {
        seams[i][j] = seams[i][j - 1];
    }
}

/* findLeastEnergySeamAndRemove
 * Remove seam with least gradient energy from input image
 */
void findLeastEnergySeamAndRemove(Image *img, Image *img_copy, int **seams, int *grad_energies) {
    // initialize variables for image dimensions and least gradient energy
    int img_copy_cols = img_copy->cols;
    int img_rows = img->rows;
    int least_grad_energy = grad_energies[0];
    int least_seam_index = 0;
    // find the index of the seam with the least gradient energy in the seams array
    for (int i = 1; i < img_copy_cols; i++) {
        if (grad_energies[i] < least_grad_energy) {
            least_grad_energy = grad_energies[i];
            least_seam_index = i;
        }
    }
    // iterate through original image and copy data until the seam with the least gradient energy
    // then skip over the seam with the least gradient energy and continue copying data
    for (int i = 0; i < img_rows; i++) {
        int col_index_to_remove = seams[least_seam_index][i];
        for (int j = 0; j < col_index_to_remove; j++) {
            img->data[i * img->cols + j] = img_copy->data[i * img_copy_cols + j];
        }
        for (int j = col_index_to_remove + 1; j < img_copy_cols; j++) {
            img->data[i * img->cols + j - 1] = img_copy->data[i * img_copy_cols + j];
        }
    }
}

/* carveSeams
 * Operate seam carving instruction for the given number of seams
 * Return 0 if successful
 */
int carveSeams(Image *img, int num_seams) {
    // called in seam for carving column and row seams
    for (int i = 0; i < num_seams; i++) {
        // create copy of image to store the gradient data
        Image *grad_img = CopyImage(img);
        // create copy of image to store original rgb values
        Image *img_copy = CopyImage(img);
        // check that both images are valid
        if (!grad_img || !img_copy) {
            return 8;
        }
        img->cols--;
        gradient(grad_img);
        // allocate space for 2D seams array
        int **seams = malloc(sizeof(int *) * grad_img->cols);
        if (seams == NULL) {
            fprintf(stderr, "Memory allocation failed.\n");
            return 8;
        }
        for (int i = 0; i < grad_img->cols; i++) {
            seams[i] = malloc(sizeof(int) * grad_img->rows);
            if (seams[i] == NULL) {
                fprintf(stderr, "Memory allocation failed.\n");
                return 8;
            }
        }
        // allocate space for gradient energies array
        int *grad_energies = malloc(sizeof(int) * grad_img->cols);
        if (grad_energies == NULL) {
            fprintf(stderr, "Memory allocation failed.\n");
            return 8;
        }
        // loop through seams and gradient energies arrays to initialize with placeholder values
        for (int i = 0; i < grad_img->cols; i++) {
            for (int j = 0; j < grad_img->rows; j++) {
                seams[i][j] = i;
                grad_energies[i] = 0;
            }
        }
        // create and remove seams from input image
        createSeamsAndGradEnergies(grad_img, seams, grad_energies);
        findLeastEnergySeamAndRemove(img, img_copy, seams, grad_energies);
        // free allocated data for images and arrays
        for (int i = 0; i < grad_img->cols; i++) {
            free(seams[i]);
        }
        free(seams);
        free(grad_energies);
        FreeImage(grad_img);
        FreeImage(img_copy);
    }
    return 0;
}

/* seam
 * Simplified version of seam carving, which maintains the aspect ratio when rescaling an image.
 * Return 0 if successful.
 */
int seam(Image *img, float scale_factor_col, float scale_factor_row) {
    // instantiate integer status to check the status of helper functions
    int status;
    // initialize input and output image dimensions
    int input_rows = img->rows;
    int input_cols = img->cols;
    int output_rows = floor(scale_factor_row * input_rows);
    int output_cols = floor(scale_factor_col * input_cols);

    // set a baseline for output image dimensions
    if (output_rows < 2) { output_rows = 2; }
    if (output_cols < 2) { output_cols = 2; }

    // calculate number of seams to remove and remove them
    int num_row_seams = input_rows - output_rows;
    int num_col_seams = input_cols - output_cols;
    // return early after each step if operation is not successful
    // i.e. status is not 0
    status = carveSeams(img, num_col_seams);
    if (status != 0) {
        return status;
    }
    status = transpose(img);
    if (status != 0) {
        return status;
    }
    status = carveSeams(img, num_row_seams);
    if (status != 0) {
        return status;
    }
    status = transpose(img);
    if (status != 0) {
        return status;
    }
    return status;
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

    int blend_img_cols = max(img_one->cols, img_two->cols);
    int blend_img_rows = max(img_one->rows, img_two->rows);
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

int randInRange(int min, int max) {
    int rand_num = rand() % (max - min + 1) + min;
    return rand_num;
}

void drawCircle(Point top_left, Point bottom_right, Point center, int radius, Image *og_img, Pixel *img_point_data) {
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
                radius = randInRange(2, 10);
                struct _point top_left = {max(j - radius, 0), max(i - radius, 0)};
                struct _point bottom_right = {min(j + radius, img->cols), min(i + radius, img->rows)};
                struct _point center = {j, i};
                drawCircle(top_left, bottom_right, center, radius, img, img_point_data);
            }
        }
    }
    
    free(img->data);
    img->data = img_point_data;

    return 0;
}