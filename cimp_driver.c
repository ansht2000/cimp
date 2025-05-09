// Ansh Tarafder atarafd1
// Alex Ma ama12
// Amy Wang awang111
#include "cimp_driver.h"
#include "ppm_io.h"

/* printUsage
 * Prints correct command line usage for user to see.
 */
void printUsage(const char *prog_name) {
    fprintf(stdout, "Usage: %s <input_file> <output_file> <operation> [<args>]\n", prog_name);
    fprintf(stdout, "Operations:\n");
    fprintf(stdout, "  grayscale\n");
    fprintf(stdout, "  binarize <threshold>\n");
    fprintf(stdout, "  crop <x1> <y1> <x2> <y2>\n");
    fprintf(stdout, "  transpose\n");
    fprintf(stdout, "  gradient\n");
    fprintf(stdout, "  seam <x_scale> <y_scale>\n");
}

/* runImgProcessing
 * Conducts input validation and calls the respective operations as entered in command line.
 * Returns 0 if successful, otherwise error occurred.
 */
int runImgProcessing(int argc, char **argv) {
    int status = 0;
    if (argc < 4) {
        printUsage(argv[0]);
        return 1;  // Wrong usage
    }
    
    const char *operation = argv[3];
    FILE *input_file = fopen(argv[1], "r");;
    FILE *input_two;
    FILE *output_file;
    if (strcmp(operation, "blend") == 0) {
        input_two = fopen(argv[2], "r");
        output_file = fopen(argv[4], "w");
    } else {
        output_file = fopen(argv[2], "w");
    }
    

    Image *img = ReadPPM(input_file);
    if (img == NULL) {
        fprintf(stderr, "Error: Unable to open input file or read as PPM.\n");
        fclose(input_file);
        fclose(output_file);
        return 2;  // Input file I/O error or cannot read as PPM
    }

    if (strcmp(operation, "grayscale") == 0) {
        if (argc != 4) {
            fprintf(stderr, "Error: Incorrect number of arguments for grayscale.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 6;  // Incorrect number of arguments
        }
        status = grayscale(img);
    } else if (strcmp(operation, "binarize") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Error: Incorrect number of arguments for binarize.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 6;  // Incorrect number of arguments
        }
        int threshold = atoi(argv[4]);
        if (threshold < 0 || threshold > 255) {
            fprintf(stderr, "Error: Invalid threshold value for binarize. Must be between 0 and 255.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 7;  // Invalid arguments for operation
        }
        status = binarize(img, threshold);
    } else if (strcmp(operation, "crop") == 0) {
        if (argc != 8) {
            fprintf(stderr, "Error: Incorrect number of arguments for crop.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 6;  // Incorrect number of arguments
        }
        int x1 = atoi(argv[4]);
        int y1 = atoi(argv[5]);
        int x2 = atoi(argv[6]);
        int y2 = atoi(argv[7]);
        if (x1 < 0 || y1 < 0 || x2 < x1 || y2 < y1 || x2 > img->cols || y2 > img->rows) {
            fprintf(stderr, "Error: Invalid crop arguments.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 7;  // Invalid arguments for operation
        }
        status = crop(img, x1, y1, x2, y2);
    } else if (strcmp(operation, "transpose") == 0) {
        if (argc != 4) {
            fprintf(stderr, "Error: Incorrect number of arguments for transpose.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 6;  // Incorrect number of arguments
        }
        status = transpose(img);
    } else if (strcmp(operation, "gradient") == 0) {
        if (argc != 4) {
            fprintf(stderr, "Error: Incorrect number of arguments for gradient.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 6;  // Incorrect number of arguments
        }
        status = gradient(img);
    } else if (strcmp(operation, "seam") == 0) {
        if (argc != 6) {
            fprintf(stderr, "Error: Incorrect number of arguments for seam.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 6;  // Incorrect number of arguments
        }
        float x_scale = atof(argv[4]);
        float y_scale = atof(argv[5]);
        if (x_scale <= 0 || x_scale > 1 || y_scale <= 0 || y_scale > 1) {
            fprintf(stderr, "Error: Invalid seam scale factors. Must be between 0 and 1.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 7;  // Invalid arguments for operation
        }
        status = seam(img, x_scale, y_scale);
    } else if (strcmp(operation, "blend") == 0) {
        if (argc != 6) {
            fprintf(stderr, "Error: Incorrect number of arguments for blend.\n");
            FreeImage(img);
            fclose(input_file);
            fclose(output_file);
            return 6;  // Incorrect number of arguments
        }
        float alpha = atof(argv[5]);
        Image *img_two = ReadPPM(input_two);
        Image *img_blend;
        status = blend(img, img_two, &img_blend, alpha);
        FreeImage(img);
        img = img_blend;

    } else {
        fprintf(stderr, "Error: Unsupported image processing operation.\n");
        FreeImage(img);
        fclose(input_file);
        fclose(output_file);
        return 5;  // Unsupported image processing operations
    }

    if (WritePPM(output_file, img) < 0) {
        fprintf(stderr, "Error: Unable to write output file.\n");
        FreeImage(img);
        fclose(input_file);
        fclose(output_file);
        return 3;  // Output file I/O error
    }

    FreeImage(img);
    fclose(input_file);
    fclose(output_file);
    return status;  // No errors detected
}