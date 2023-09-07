README.md

This code provides functions to decompress and save images in PGM format from HDF5 datasets using bitshuffle filter.

## Code Overview

The code consists of several functions and data structures. Here's an overview of each component:

### Data Structures

1. `image`: Represents an image with properties such as data, width, height, and sigma.
2. `Hdf5Header`: Stores information related to the HDF5 header, including distance, wavelength, beam center, and size.
3. `Spot`: Represents a connected component in the image, containing a label, number of pixels, and arrays to store the coordinates of the pixels.

## Requirements
* C compiler
* HDF5 library
* Bitshuffle HDF5 filter library [included]

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
