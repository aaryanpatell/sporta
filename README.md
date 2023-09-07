README.md

This code provides functions to decompress and save images in PGM format from HDF5 datasets using bitshuffle filter.

## Code Overview

The code consists of several functions and data structures. Here's an overview of each component:

### Data Structures

1. `image`: Represents an image with properties such as data, width, height, and sigma.
2. `Hdf5Header`: Stores information related to the HDF5 header, including distance, wavelength, beam center, and size.
3. `Spot`: Represents a connected component in the image, containing a label, number of pixels, and arrays to store the coordinates of the pixels.

### Functions

1. `saveImagePGM(float *data, int width, int height)`: Saves the image data in PGM format to a file named "output.pgm".
2. `decompressHDF5(const char* filename, const char* dataset_path, int index)`: Reads and decompresses a dataset from an HDF5 file. Returns the decompressed data as a float array.
3. `readHDF5(const char* filename)`: Reads specific information from an HDF5 file, such as the distance, wavelength, and beam center.
4. `create_spot(int label, int num_pixels)`: Creates a `Spot` structure with the given label and number of pixels.
5. `print_spot(const Spot* spot)`: Prints the details of a `Spot` structure to the console.
6. `print_spot_to_file(const Spot* spot, FILE* file)`: Writes the details of a `Spot` structure to a file.
7. `draw_circle_around_spot(float *image, int width, int height, Spot *spot, unsigned char circle_color)`: Draws a circle around a given spot in the image.
8. `destroy_spot(Spot* spot)`: Deallocates the memory used by a `Spot` structure.
9. `dfs(int x, int y, int current_label, int height, int width, float *image, int* num_pixels, int* pixels)`: Depth-first search algorithm to perform connected component analysis.
10. `connectedComponentAnalysis(int height, int width, float *image, int* num_spots)`: Performs connected component analysis on the image and returns an array of `Spot` structures.
11. `calculateAverageIntensity(float* image, int width, int height)`: Calculates the average intensity of the image.
12. `applySPORTA(float* image, int height, int width)`: Applies the SPORTA algorithm to the image, performing connected component analysis and drawing circles around the spots.
13. `applyThreshold(float* image, int width, int height)`: Applies a thresholding operation to the image based on the average intensity.
14. `applyGaussianBlur(float* image, int width, int height, float sigma)`: Applies a Gaussian blur to the image using a given sigma value.

## Requirements
* C compiler
* HDF5 library
* Bitshuffle HDF5 filter library [bshuf_h5filter files]

## Installation

1. Install the HDF5 library. For example, on Ubuntu, you can use the following command:

```css
sudo apt-get install libhdf5-dev
```
or

```css
sudo yum install hdf5-devel.x86_64
```

2. Download the `bshuf_h5filter.h` and `bshuf_h5filter.c` files and place them in the `src` directory. Ensure 'src' is placed within the same folder as the program.

3. To compile the program, use the 'Makefile' using the following command:

```css
make
```

## Usage


Note that `width` and `height` must be provided to `save_image_as_pgm`.

## License

TO BE DECIDED
