#ifndef PPM_IO_H
#define PPM_IO_H

#include <stdio.h> // c file type: FILE
#include <stdint.h>

/* A struct to store a point (2D coordinate).
 */ 
typedef struct _point {
  int x;
  int y;
} Point;

/* A struct to store a single RGB pixel, one byte per color channel.
 */
typedef struct _pixel {
  uint8_t r;
  uint8_t g;
  uint8_t b;
} Pixel;

/* A struct to bundle together a pixel array with the other
 * image data we'll frequently want to pass around with it.
 * (This saves us from having to pass the same three 
 * variables to every function.) Note that no Pixels are
 * stored within this struct; the data field is a pointer.
 */
typedef struct _image {
  Pixel *data;  // pointer to array of Pixels
  int rows;     // number of rows of Pixels
  int cols;     // number of columns of Pixels
} Image;

/* ReadPPM
 * Read a PPM-formatted image from a file (assumes fp != NULL).
 * Returns the address of the heap-allocated Image struct it
 * creates and populates with the Image data.
 */
Image* ReadPPM(FILE *fp);

/* WritePPM
 * Write a PPM-formatted image to a file (assumes fp != NULL),
 * and return the number of pixels successfully written.
 */
int WritePPM(FILE *fp, const Image *img);

/* NewImage
 * Create a new image with r rows and c columns.
 * Returns the address of the new heap-allocated Image.
 */
Image* NewImage(int r, int c);

/* CopyImage
 * Create a copy of Image img.
 * Returns the address of the copied heap-allocated Image.
 */
Image* CopyImage (Image *img);

/* FreeImage
 * Free the memory allocation for Image img.
 */
void FreeImage (Image *img);

#endif // PPM_IO_H
