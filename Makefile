CC=gcc
CFLAGS=-std=c99 -pedantic -Wall -Wextra -g

## Below are commands to link and compile the checkerboard program
# Links together files needed to create the checkerboard executable
checkerboard : checkerboard.o ppm_io.o
	$(CC) -o $@ checkerboard.o ppm_io.o

cimp : cimp.o img_processing.o ppm_io.o cimp_driver.o
	$(CC) -o $@ cimp.o img_processing.o ppm_io.o cimp_driver.o -lm

# Compile the ppm i/o source code
ppm_io.o: ppm_io.c ppm_io.h
	$(CC) $(CFLAGS) -c -g ppm_io.c

img_processing.o: img_processing.c img_processing.h ppm_io.h
	$(CC) $(CFLAGS) -c -g img_processing.c

cimp.o: cimp.c ppm_io.h img_processing.h cimp_driver.h
	$(CC) $(CFLAGS) -c -g cimp.c

cimp_driver.o: cimp_driver.c cimp_driver.h ppm_io.h img_processing.h
	$(CC) $(CFLAGS) -c -g cimp_driver.c

checkerboard.o: checkerboard.c ppm_io.h
	$(CC) $(CFLAGS) -c -g checkerboard.c

# Removes all object files and the executable named cimp, so we can start fresh
clean:
	rm -f *.o *.ppm checkerboard cimp 

# Zips important files
zip: 
	zip midterm.zip img_processing.c img_processing.h Makefile ppm_io.c ppm_io.h cimp.c cimp_driver.c cimp_driver.h README gitlog.txt
