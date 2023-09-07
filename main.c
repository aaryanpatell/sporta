/**
    SPORTA: Spot Detection & Screening of X-Ray Diffraction Images​
    Copyright (c) 2023, Aaryan Patel

    MIT LICENSE

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    Acknowledgments:
    Canadian Light Source 
    Industrial Science
    Denis Spasyuk

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <unistd.h>
#include <bshuf_h5filter.h>
#include <math.h>
#include <dirent.h>
#include <regex.h>
#include <ctype.h>
#include <sys/wait.h>

#define UNCLASSIFIED -1
#define UNDEFINED -1
#define NOISE -2
#define CORE_POINT 0

typedef struct {
    double *data;
    int width;
    int height;
    double** pixels;
} Image;

typedef struct {
    double distance;
    double wavelength;
    double beam_center[2];
    double size[2];
    double psize;
    double threshold;
    int radius_90_percent;
    double resolution;
    int num_spots;
    double intensity;
    int ice;
    double max_distance;
    double snr;
} Hdf5Header;

typedef struct {
    double thresfactor; //intensity scaling factor
    double proximity; //pixel exclusion radius
    int index; // frame number
    bool output; // frame number
    char *filename; //data file name
    char *master; //master file name
    bool traversing; //traversing directories 
    char *data;
    float percent;
    char *searchdir;
    bool gnu;
    double gaussian;
    char path;
} config;

typedef struct {
    int label;
    int num_pixels;
    int* x_coords;
    int* y_coords;
    double intensity;
    int cluster_id;
    double reachability_distance;
    double core_distance;
} Spot;


typedef struct PositionSet {
    double startPosition;
    double endPosition;
    struct PositionSet* next;
} PositionSet;

PositionSet* positionSets = NULL;
Hdf5Header header;
Image img;
config cfg;

// directions for 4-connectivity (up, down, left, right)
int dx[4] = {-1, 1, 0, 0};
int dy[4] = {0, 0, -1, 1};

void printconfig(config cfg) {
    printf("\033[1;33mFilename: %s\033[0m\n", cfg.filename);
    printf("Master: %s\n", cfg.master);
    printf("Proximity: %f\n", cfg.proximity);
    printf("Index: %d\n", cfg.index);
    printf("Output: %d\n", cfg.output);
    printf("Traversing: %s\n", cfg.traversing ? "true" : "false");
    printf("Data: %s\n", cfg.data);
    printf("Search Directory: %s\n", cfg.searchdir);
}

void printheader(const Hdf5Header* header) {
    printf("=================================HEADER============================== \n");
    printf("Distance: %f\n", header->distance);
    printf("Wavelength: %f\n", header->wavelength);
    printf("Beam Center: (%f, %f)\n", header->beam_center[0], header->beam_center[1]);
    printf("Size: (%f, %f)\n", header->size[0], header->size[1]);
    printf("Pixel Size: %f\n", header->psize);
    printf("==================================INFO=============================== \n");
    printf("Threshold: %f\n", header->threshold);
    printf("SNR: %f\n", header->snr);
    printf("Max Radius: %f\n", header->max_distance);
    printf("\033[1;31mResolution: %f\033[0m\n", header->resolution);
    printf("\033[1;31mSpot Count: %d\033[0m\n", header->num_spots);  // assuming num_spots is pointing to a valid integer
    printf("Resolution: %f \n", header->resolution); 
    PositionSet* currentSet = positionSets;
    int groupNumber = 1;

    if (currentSet == NULL){
        printf("Ice Ring: NONE\n");
    }
    else{
        while (currentSet != NULL) {
            printf("Ice Ring %d: %f-%f\n", groupNumber, currentSet->startPosition, currentSet->endPosition);
            currentSet = currentSet->next;
            groupNumber++;
        }

        while (positionSets != NULL) {
            PositionSet* temp = positionSets;
            positionSets = positionSets->next;
            free(temp);
        }
    }
}

void mainprint(){
    printf("=================================DATASET============================= \n");
    printconfig(cfg);
    printheader(&header);
    printf("===================================================================== \n");
    
}

void noteprint(){
}

void savepgm(const char* filename, double* data, int width, int height) {
    char output_filename[64];
    snprintf(output_filename, sizeof(output_filename), "%s.pgm", filename);
    FILE* file = fopen(output_filename, "w");
    if (file == NULL) {
        printf("Error opening file for writing: %s\n", output_filename);
        return;
    }

    // Write the PGM header
    fprintf(file, "P2\n");
    fprintf(file, "%d %d\n", width, height);
    fprintf(file, "255\n");
    if (filename == "pre-output"){
    // Write the inverted image data
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int value = (int)(data[i * width + j]);
                if (value < 0 || value > 1000) {
                    value = 0;
                }
                fprintf(file, "%d ", value);
            }
            fprintf(file, "\n");
        }
    }
    else {
        // Write the inverted image data
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int value = (int)(data[i * width + j]);
                if (data[i * width + j] > 10000 && data[i * width + j] < 0) {
                    value = 0;
                }
                fprintf(file, "%d ", value);
            }
            fprintf(file, "\n");
        }
    }
    fclose(file);
}


double* decompresshdf5(const char* filename, const char* dataset_path, int index) {
    H5Eset_auto(H5E_DEFAULT, (H5E_auto_t)H5Eprint, stderr);
    herr_t status = bshuf_register_h5filter();
    if (status < 0) {
        fprintf(stderr, "Error registering Bitshuffle filter\n");
        return NULL;
    }

    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return NULL;
    }

    hid_t dataset_id = H5Dopen(file_id, dataset_path, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error opening dataset: %s\n", dataset_path);
        H5Fclose(file_id);
        return NULL;
    }

    // Check if the dataset is compressed with Bitshuffle
    hid_t filter_id;
    unsigned int flags;
    hid_t dataspace_id = H5Dget_space(dataset_id);
    hsize_t dims[3];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    int num_images = (int)dims[0];

    if (index >= num_images) {
        fprintf(stderr, "Index out of bounds\n");
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return NULL;
    }

    hsize_t start[3] = {index, 0, 0};
    hsize_t count[3] = {1, dims[1], dims[2]};
    hid_t mem_dataspace = H5Screate_simple(3, count, NULL);
    double* data_out = (double*)malloc(dims[1] * dims[2] * sizeof(double));

    if (data_out == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        H5Sclose(mem_dataspace);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return NULL;
    }

    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        fprintf(stderr, "Error selecting hyperslab\n");
        H5Sclose(mem_dataspace);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        free(data_out);
        return NULL;
    }

    if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, dataspace_id, H5P_DEFAULT, data_out) < 0) {
        fprintf(stderr, "Error reading data\n");
        H5Sclose(mem_dataspace);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        free(data_out);
        return NULL;
    }

    H5Sclose(mem_dataspace);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    return data_out;
}

void readhdf5(const char* filename){
    hid_t file_id, data_id, detector_id, dataset_id, detector_size_g, sample_id, gonio_id, beam_id, distance_id, detectorX_size,
        detectorY_size, wavelength_id, size_id, tector_x, detector_y, beamX_id, beamY_id, pixelX_id, detectorS_id, psize_id;
    
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    data_id = H5Gopen(file_id, "/entry/data", H5P_DEFAULT);
    detector_id = H5Gopen(file_id, "/entry/instrument/detector", H5P_DEFAULT);
    detectorS_id = H5Gopen(file_id, "/entry/instrument/detector/detectorSpecific", H5P_DEFAULT);
    beam_id = H5Gopen(file_id, "/entry/instrument/beam", H5P_DEFAULT);
    
    detectorX_size = H5Dopen(detectorS_id, "x_pixels_in_detector", H5P_DEFAULT);
    H5Dread(detectorX_size, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.size[0]);  
    H5Dclose(detectorX_size);
    detectorY_size = H5Dopen(detectorS_id, "y_pixels_in_detector", H5P_DEFAULT);
    H5Dread(detectorY_size, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.size[1]);  
    H5Dclose(detectorY_size);
    distance_id = H5Dopen2(detector_id, "detector_distance", H5P_DEFAULT);
    H5Dread(distance_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.distance);  
    H5Dclose(distance_id);
    wavelength_id = H5Dopen2(beam_id, "incident_wavelength", H5P_DEFAULT);
    H5Dread(wavelength_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.wavelength);
    H5Dclose(wavelength_id);
    beamX_id = H5Dopen2(detector_id, "beam_center_x", H5P_DEFAULT);
    H5Dread(beamX_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.beam_center[0]);
    H5Dclose(beamX_id);
    beamY_id = H5Dopen2(detector_id, "beam_center_y", H5P_DEFAULT);
    H5Dread(beamY_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.beam_center[1]);
    H5Dclose(beamY_id);
    psize_id = H5Dopen2(detector_id, "x_pixel_size", H5P_DEFAULT);
    H5Dread(psize_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.psize);
    H5Dclose(psize_id);
}

Spot* create_spot(int label, int num_pixels) {
    Spot* spot = malloc(sizeof(Spot));
    spot->label = label;
    spot->num_pixels = num_pixels;
    spot->x_coords = malloc(num_pixels * sizeof(int));
    spot->y_coords = malloc(num_pixels * sizeof(int));
    spot->core_distance = UNDEFINED; // Initialize core_distance
    spot->reachability_distance = UNDEFINED; // Initialize reachability_distance
    spot->cluster_id = UNDEFINED; // Initialize cluster_id
    return spot;
}

void filespot(const Spot* spot, int x, int y, FILE* file) {
    fprintf(file, "%4d  %5d     %5d,%5d %8.f     %d\n",
            spot->label, spot->num_pixels, y, x, spot->intensity, spot->cluster_id);
}

void drawcircle(double *image, int spot_radius, int width, int height, int center_x, int center_y, unsigned char circle_color) {
    int thickness = 3;
    for (int x = center_x - spot_radius; x <= center_x + spot_radius; x++) {
        for (int y = center_y - spot_radius; y <= center_y + spot_radius; y++) {
            // Calculate the distance from the center pixel
            int dx = x - center_x;
            int dy = y - center_y;
            int distance = dx * dx + dy * dy;

            bool inside_condition = distance <= spot_radius * spot_radius && distance > (spot_radius - thickness) * (spot_radius - thickness);

            // If the pixel is within the circle radius and on the circle's circumference, set it to the circle color
            if (inside_condition && x >= 0 && x < height && y >= 0 && y < width) {
                image[x * width + y] = circle_color;
            }
        }
    }
}

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

double calculateresolution(double radius) {
    return (header.wavelength / (sin(atan(radius * header.psize / header.distance) / 2) * 2)); 
}

void applysingulargaussian(double data[], int size, double sigma) {
    double* tempData = malloc(size * sizeof(double));
    memcpy(tempData, data, size * sizeof(double));

    double* kernel = malloc((6 * sigma + 1) * sizeof(double));
    double kernelSum = 0.0;

    int i;
    for (i = 0; i < 6 * sigma + 1; i++) {
        double x = i - 3 * sigma;
        kernel[i] = exp(-0.5 * x * x / (2 * (sigma * sigma)));
        kernelSum += kernel[i];
    }

    for (i = 0; i < 6 * sigma + 1; i++) {
        kernel[i] /= kernelSum;
    }

    for (i = 0; i < size; i++) {
        double blurredValue = 0.0;
        int j;
        for (j = 0; j < 6 * sigma + 1; j++) {
            int index = i + j - 3 * sigma;
            if (index >= 0 && index < size) {
                blurredValue += tempData[index] * kernel[j];
            }
        }
        data[i] = blurredValue;
    }

    free(tempData);
    free(kernel);
}


void addringset(PositionSet** head, double startPosition, double endPosition) {
    PositionSet* newPositionSet = (PositionSet*)malloc(sizeof(PositionSet));
    newPositionSet->startPosition = startPosition;
    newPositionSet->endPosition = endPosition;
    newPositionSet->next = *head;
    *head = newPositionSet;
}

double calculateavgintensity(double* image, int width, int height, double thresfactor) {
    int totalPixels = width * height;
    int sum = 0;
    for (int i = 0; i < totalPixels; i++) {
        unsigned char pixelValue = image[i];
        sum += pixelValue;
    }
    return (double) sum / totalPixels * thresfactor;
}


void intensityanalysis(double data[], int width, int height) {
    header.resolution = calculateresolution(header.max_distance);

    int newSize = (img.height + img.width) / 2;  // Modify the variable name to avoid confusion
    double max = 10; // set maximum resolution
    double mean = calculateavgintensity(img.data, img.width, img.height, 1);  // set threshold requirement for exclusion
    double ceiling = mean * 150;

    unsigned char circle_color = 255;

    freopen("/dev/null", "w", stderr);

    FILE* data_file = fopen("histogram_data.txt", "w");
    if (data_file == NULL) {
        fprintf(stderr, "Error opening temporary data file\n");
        return;
    }
    double sigma = 4;  // Adjust the sigma value as desired
    applysingulargaussian(data, newSize, sigma);

    if (cfg.gnu){
        for (int i = 0; i < newSize; i++) {
            double position = calculateresolution(i);
            int total = data[i];
            if (position < max) {
                fprintf(data_file, "%f %d\n", position, total);
            }
        }

        fclose(data_file);

        FILE* gnuplot_pipe = popen("gnuplot -persist", "w");
        if (gnuplot_pipe == NULL) {
            fprintf(stderr, "Error opening Gnuplot pipe\n");
            return;
        }
        fprintf(gnuplot_pipe, "set terminal qt title '%s'\n", cfg.filename);
        fprintf(gnuplot_pipe, "set title 'Powder X-Ray Diffraction Pattern'\n");
        fprintf(gnuplot_pipe, "set xlabel 'Resolution [Å]'\n");
        fprintf(gnuplot_pipe, "set ylabel 'Intensity [8-bit Greyscale]'\n");
        fprintf(gnuplot_pipe, "set xrange noextend\n");
        fprintf(gnuplot_pipe, "set style line 1 lt 1 lc rgb 'black' lw 1\n");
        fprintf(gnuplot_pipe, "plot 'histogram_data.txt' with filledcurves above y=0 ls 1 title 'Intensity'\n");
        fprintf(gnuplot_pipe, "exit\n");
        pclose(gnuplot_pipe);

        remove("histogram_data.txt");
    }
    double start = -1;
    double position;

    for (int i = 0; i < newSize; i++) {
        position = calculateresolution(i);
        int total = data[i];
        if (total > ceiling) {
            if (start == -1) {
                start = position;
            }
        } else {
            if (start != -1) {
                double end = calculateresolution(i - 1);
                addringset(&positionSets, start, end);
                start = -1;
            }
        }
    }

    if (start != -1) {
        double end = calculateresolution(newSize - 1);
        addringset(&positionSets, start, end);
    }
}


void dfs(int x, int y, int current_label, int height, int width, double *image, int* num_pixels, int* pixels) {
    if (x < 0 || x >= height || y < 0 || y >= width || image[x * width + y] != 1 || *num_pixels > 65) {
        return;
    }

    // Label the current pixel
    image[x * width + y] = current_label;

    // Store the pixel coordinates
    pixels[(*num_pixels) * 2] = x;
    pixels[(*num_pixels) * 2 + 1] = y;
    (*num_pixels)++;

    // Visit the neighbors
    for (int direction = 0; direction < 4; direction++) {
        dfs(x + dx[direction], y + dy[direction], current_label, height, width, image, num_pixels, pixels);
    }
}

Spot** connectedcomponentanalysis(int height, int width, double *image, double *data, double *preradial, double *postradial, int* num_spots) {
    int next_label = 0;
    int min_pixel = 4;
    int max_pixel = 60;
    Spot** spots = NULL;
    int max_spots = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (image[i * width + j] == 1) {
                int distance = (int) sqrt(pow(header.beam_center[0] - j, 2) + pow(header.beam_center[1] - i, 2));
                int current_label = next_label;
                int num_pixels = 0;
                double intensity = 0;
                int* pixels = malloc(height * width * 2 * sizeof(double));
                bool scdcondition = false;
                int encounter = 0;
                dfs(i, j, current_label, height, width, image, &num_pixels, pixels);
                for (int q = -cfg.proximity; q <= cfg.proximity; q++) {
                    for (int w = -cfg.proximity; w <= cfg.proximity; w++) {
                        int checkX = q + i;
                        int checkY = w + j;
                        if (checkX >= 0 && checkX < width && checkY >= 0 && checkY < height) {
                            if (image[checkX * width + checkY] != 0 && image[checkX * width + checkY] != current_label) {
                                encounter++;
                                if (encounter > 0) {
                                    scdcondition = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (scdcondition) {
                        break;
                    }
                }
                if (!scdcondition && num_pixels >= min_pixel && num_pixels <= max_pixel) {
                    if (postradial[distance] == 0){
                        postradial[distance] = data[i * width + j];
                    }
                    else {
                        postradial[distance] += data[i * width + j];
                    }
                    if (*num_spots >= max_spots) {
                        max_spots += 10;
                        spots = realloc(spots, max_spots * sizeof(Spot*));
                    }
                    spots[*num_spots] = create_spot(current_label, num_pixels);
                    spots[*num_spots]->x_coords = malloc(num_pixels * sizeof(double));
                    spots[*num_spots]->y_coords = malloc(num_pixels * sizeof(double));
                    for (int k = 0; k < num_pixels; k++) {
                        spots[*num_spots]->x_coords[k] = pixels[k * 2];
                        spots[*num_spots]->y_coords[k] = pixels[k * 2 + 1];
                        intensity += data[pixels[k * 2] * width + pixels[k * 2 + 1]];
                    }
                    spots[*num_spots]->intensity = intensity;
                    (*num_spots)++; 
                }
                
                free(pixels);
                next_label++;
            }
        }
    }
    return spots;
}

int comparedistances(const void* a, const void* b) {
    double distanceA = *(const double*)a;
    double distanceB = *(const double*)b;

    if (distanceA < distanceB) {
        return -1;
    } else if (distanceA > distanceB) {
        return 1;
    } else {
        return 0;
    }
}

double euclideandistance(Spot* a, Spot* b) {
    double dx = a->x_coords[0] - b->x_coords[0];
    double dy = a->y_coords[0] - b->y_coords[0];
    return sqrt(dx * dx + dy * dy);
}

void updatecoredistance(Spot* spot, Spot** neighbors, int num_neighbors, double epsilon, unsigned int minPts) {
    if (num_neighbors >= minPts) {
        double* distances = (double*)malloc(num_neighbors * sizeof(double));
        for (int i = 0; i < num_neighbors; ++i) {
            distances[i] = euclideandistance(spot, neighbors[i]);
        }
        qsort(distances, num_neighbors, sizeof(double), comparedistances);

        spot->core_distance = distances[minPts - 1] + epsilon;

        free(distances);
    } else {
        spot->core_distance = UNDEFINED;
    }
}

void optics(Spot** spots, int num_spots, double epsilon, unsigned int minPts) {
    for (int i = 0; i < num_spots; ++i) {
        if (spots[i]->cluster_id != UNCLASSIFIED) {
            continue;
        }

        Spot** neighbors = (Spot**)malloc(num_spots * sizeof(Spot*));
        int num_neighbors = 0;

        for (int j = 0; j < num_spots; ++j) {
            if (i == j) {
                continue;
            }

            double distance = euclideandistance(spots[i], spots[j]);
            if (distance <= epsilon) {
                neighbors[num_neighbors] = spots[j];
                ++num_neighbors;
            }
        }

        spots[i]->cluster_id = NOISE;
        spots[i]->reachability_distance = UNDEFINED;

        if (num_neighbors >= minPts) {
            updatecoredistance(spots[i], neighbors, num_neighbors, epsilon, minPts);

            for (int j = 0; j < num_neighbors; ++j) {
                Spot* neighbor = neighbors[j];
                if (neighbor->cluster_id == UNCLASSIFIED) {
                    double reachability_distance = euclideandistance(spots[i], neighbor);
                    if (reachability_distance <= epsilon) {
                        neighbor->reachability_distance = fmax(reachability_distance, spots[i]->core_distance);
                    }
                }
            }

            spots[i]->cluster_id = CORE_POINT;
        }

        free(neighbors);
    }
}

double calculate_snr(int* histogram, int num_bins, int threshold, int totalPixels) {
    // Calculate the mean intensity of the signal (μ_signal) within the object region
    double mean_signal = 0.0;
    int num_signal_pixels = 0;

    for (int i = threshold; i < num_bins; ++i) {
        mean_signal += i * histogram[i];
        num_signal_pixels += histogram[i];
    }

    mean_signal /= num_signal_pixels;

    // Calculate the standard deviation of the background noise (σ_noise) within the background region
    double variance_noise = 0.0;
    int num_noise_pixels = 0;

    for (int i = 0; i < threshold; ++i) {
        variance_noise += pow(i - mean_signal, 2) * histogram[i];
        num_noise_pixels += histogram[i];
    }

    variance_noise /= num_noise_pixels;
    double stddev_noise = sqrt(variance_noise);
    double snr = mean_signal / stddev_noise;

    return snr;
}


void applySPORTA(double* image, double* data, int height, int width){
    int num_spots = 0;
    double preradial[img.height + img.width], postradial[img.height + img.width];
    Spot** spots = connectedcomponentanalysis(height, width, image, data, preradial, postradial, &num_spots);
    header.num_spots = num_spots;

    double epsilon = 90;
    int minPts = 1;
    optics(spots, num_spots, epsilon, minPts);

    unsigned char circle_color = 255;
    FILE* file = fopen("output.txt", "w");
    fprintf(file, "Label  Pixels    Coordinate   Intensity   ClusterID\n");
    double pX = 0;
    double tightness = 100;
    double distance = 0;
    double max_distance = 0;

    Spot* furthest_spot = NULL;
    double max_distance_from_center = 0;

    long totalintensity = 0;

    for (int i = 0; i < num_spots; i++) {
        int totalX = 0;
        int totalY = 0;
        int num = spots[i]->num_pixels;
        for (int k = 0; k < num; k++) {
            totalintensity += spots[i]->intensity;
            totalX += spots[i]->x_coords[k];
            totalY += spots[i]->y_coords[k];
        }
        int aX = (int)(totalX / num);
        int aY = (int)(totalY / num);
        if (spots[i]->cluster_id == 0) {
            double spot_distance_from_center = sqrt(pow(header.beam_center[1] - aX, 2) + pow(header.beam_center[0] - aY, 2));
            if (spot_distance_from_center > max_distance_from_center) {
                max_distance_from_center = spot_distance_from_center;
                furthest_spot = spots[i];
            }
        }

        if (cfg.output){
            filespot(spots[i], aX, aY, file);
            drawcircle(image, 10, width, height, aX, aY, circle_color);
            drawcircle(data, 10, width, height, aX, aY, 255);
        }
    }
    header.max_distance = max_distance_from_center * 0.95;
    if (cfg.output){
        drawcircle(image, max_distance_from_center * 0.95, width, height, header.beam_center[1], header.beam_center[0], circle_color);
    }
    intensityanalysis(postradial, width, height);
    fclose(file); 
    free(spots);
}


void applythreshold(double* image, double* data, int width, int height) {
    int histogram[256] = {0};
    for (int i = 0; i < height * width; ++i) {
        int pixelValue = (int)image[i];
        if (pixelValue >= 0 && pixelValue <= 255) {
            histogram[pixelValue]++;
        }
    }

    int cumulative[256] = {0};
    cumulative[0] = histogram[0];
    for (int i = 1; i < 256; ++i) {
        cumulative[i] = cumulative[i - 1] + histogram[i];
    }
    // Calculate the total number of pixels
    int totalPixels = height * width;

    double threshold = 0.0;
    double w1 = (double)cumulative[1] / totalPixels;
    double w2 = 1.0 - w1;

    double sum1 = 0.0;
    for (int j = 0; j <= 1; ++j) {
        sum1 += j * histogram[j];
    }
    threshold = sum1 / (w1 * cumulative[1]);

    double sum2 = 0.0;
    for (int j = 2; j < 256; ++j) {
        sum2 += j * histogram[j];
    }
    threshold += sum2 / (w2 * (totalPixels - cumulative[1]));
    header.threshold = threshold;
    header.snr = calculate_snr(histogram, 256, threshold, width * height);


    // Apply the multi-otsu thresholds
    for (int i = 0; i < height * width; ++i) {
        int pixelValue = (int)image[i];

        // If the pixel intensity is greater than the threshold and less than 1000, set it to 1
        if (pixelValue > threshold && pixelValue < 1000) {
            image[i] = 1.0;
        }
        // Otherwise, set it to 0
        else {
            image[i] = 0.0;
        }
    }
}


void applygaussianblur(double* image, int width, int height, double sigma) {
    // Create the Gaussian kernel
    sigma = cfg.gaussian;
    int kernel_size = ceil(6 * sigma) + 1;
    int radius = kernel_size /2;
    float* kernel = (float*)malloc(kernel_size * sizeof(float));
    float kernel_sum = 0.0;

    for (int i = 0; i < kernel_size; i++) {
        int x = i - radius;
        kernel[i] = exp(-(x * x) / (3 * sigma * sigma));
        kernel_sum += kernel[i];
    }

    // Normalize the kernel
    for (int i = 0; i < kernel_size; i++) {
        kernel[i] /= kernel_sum;
    }

    // Create a temporary buffer to store the blurred image
    float* blurred_image = (float*)malloc(width * height * sizeof(float));

    // Convolve each pixel with the Gaussian kernel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float sum = 0.0;
            int count = 0;

            // Apply the kernel to the current pixel
            for (int i = 0; i < kernel_size; i++) {
                int image_x = x + i - radius;
                if (image_x >= 0 && image_x < width) {
                    sum += image[y * width + image_x] * kernel[i];
                    count++;
                }
            }

            // Store the blurred pixel in the temporary buffer
            blurred_image[y * width + x] = sum / count;
        }
    }

    // Copy the blurred image back to the original image buffer
    for (int i = 0; i < width * height; i++) {
        image[i] = blurred_image[i];
    }

    // Free the memory
    free(kernel);
    free(blurred_image);
}
char* generatemaster(const char* filename) {
    // Check if the input file exists
    if (access(filename, F_OK) == -1) {
        printf("File %s does not exist.\n", filename);
        return NULL;
    }

    // Copy the filename to a new variable that we can modify
    static char master[256];
    strncpy(master, filename, sizeof(master));
    master[sizeof(master) - 1] = '\0';  // Ensure null termination

    // Find the last occurrence of "_data_" in the filename
    char *substring = strstr(master, "_data_");
    if (substring == NULL) {
        printf("Invalid filename format.\n");
        return NULL;
    }

    // Replace "_data_" and everything after it with "_master.h5"
    strncpy(substring, "_master.h5", strlen("_master.h5"));
    // Make sure the new filename string is properly terminated
    substring[strlen("_master.h5")] = '\0';

    // Check if the new file exists
    if (access(master, F_OK) == -1) {
        printf("File %s does not exist.\n", master);
        return NULL;
    }

    return master;
}

double* deepcopyimage(const double* source, int width, int height) {
    int dataSize = width * height;
    double* target = (double*)malloc(dataSize * sizeof(double));

    for (int i = 0; i < dataSize; i++) {
        target[i] = source[i];
    }

    return target;
}

void setdefault(){
    header.size[0] = 3110;
    header.size[1] = 3269;
    header.distance = 0.2;
    header.wavelength = 0.953742;
    header.psize = 0.000075;
    header.beam_center[0] =1544;
    header.beam_center[1] = 1577;
}

int ends_with_data(const char *path, char *folder) {
    char *last_slash = strrchr(path, '/');
    if (last_slash != NULL && strcmp(last_slash + 1, folder) == 0) {
        return 1;
    }
    return 0;
}

int contains_proc(const char *name) {
    return strstr(name, "proc") != NULL;
}

int contains_native(const char *name) {
    return strstr(name, "native") != NULL;
}

int screener(const char *filename){
    cfg.master = generatemaster(cfg.filename);
    if (cfg.master != NULL){
        readhdf5(cfg.master);
    }
    else{
        setdefault();
    }
    
    double *data = decompresshdf5(cfg.filename, cfg.data, cfg.index);
    img.height = (int)header.size[0];
    img.width = (int)header.size[1];
    img.data = deepcopyimage(data, img.width, img.height);
    cfg.proximity = 2;
    if (cfg.gaussian > 0) applygaussianblur(img.data, img.width, img.height, cfg.gaussian);
    applythreshold(img.data, data, img.width, img.height);
    applySPORTA(img.data, data, img.width, img.height);
    if (cfg.output){
        savepgm("pre-output", data, img.height, img.width);
        savepgm("post-output", img.data, img.height, img.width);
    }
    mainprint();
    free(data);
    free(img.data);

    return 0;
}

int matches_pattern(const char *filename) {
    regex_t regex;
    int ret;
    if (regcomp(&regex, ".*_000001\\.h5$", REG_EXTENDED) != 0) {
        printf("Could not compile regex\n");
        return 0;
    }
    ret = regexec(&regex, filename, 0, NULL, 0);
    regfree(&regex);
    if (ret == 0) {
        return 1;
    } else {
        return 0;
    }
}

int getfile(const char *dirname) {
    DIR *dir;
    struct dirent *ent;

    if ((dir = opendir(dirname)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (ent->d_type == DT_REG) {
                char *filename = ent->d_name;
                if (matches_pattern(filename)) {
                    char filepath[1024];
                    snprintf(filepath, sizeof(filepath), "%s/%s", dirname, filename);
                    printf("Processing file: %s\n", filepath);
                    cfg.filename = filepath;
                    screener(cfg.filename);
                }
            }
        }
        closedir(dir);
        return 0;
    } else {
        perror("Could not open directory");
        return -1;
    }
}

char** getdirs(const char* currentdir, char* folder, int* count) {
    struct dirent* dir_entry;
    DIR* dir = opendir(currentdir);

    if (dir == NULL) {
        perror("Unable to read directory");
        return NULL;
    }

    char** directories = NULL;
    *count = 0;

    while ((dir_entry = readdir(dir)) != NULL) {
        if (dir_entry->d_type == DT_DIR) {
            if (strcmp(dir_entry->d_name, ".") == 0 || strcmp(dir_entry->d_name, "..") == 0)
                continue;
            char path[1024];
            snprintf(path, sizeof(path), "%s/%s", currentdir, dir_entry->d_name);
            if (ends_with_data(path, folder) && !contains_proc(path) && !contains_native(path)) {
                (*count)++;
                directories = realloc(directories, sizeof(char*) * (*count));
                directories[*count - 1] = strdup(path);
            }

            int child_count = 0;
            char** child_directories = getdirs(path, folder, &child_count);
            if (child_directories != NULL) {
                directories = realloc(directories, sizeof(char*) * (*count + child_count));
                memcpy(directories + *count, child_directories, sizeof(char*) * child_count);
                *count += child_count;
                free(child_directories);
            }
        }
    }
    closedir(dir);
    return directories;
}

int argparse(int argc, char **argv) {
    cfg.master = NULL;
    cfg.output = false;
    cfg.index = 1;
    cfg.data = "/entry/data/data";
    cfg.thresfactor = 6;
    cfg.proximity = 2;
    cfg.percent = 0.95;
    cfg.traversing = 0;  
    cfg.searchdir = "data";
    cfg.gnu = -1;
    header.ice = 0;
    int c;

    while ((c = getopt(argc, argv, "hof:i:d:t:s:rp:g:")) != -1) {
        switch (c) {
            case 'h':
                printf("Usage: ./imager -f <filepath> -o\n");
                printf("Options:\n");
                printf("  -h          Display this help message\n");
                printf("  -f <file>   Specify the input file\n");
                printf("  -o          Enable output\n");
                printf("  -i <index>  Set the index to a specified value\n");
                printf("  -d <data>   Set the data to a specified value\n");
                printf("  -r <dir>    Enable traversing on path <dir>\n");
                printf("  -t <value>  Set the threshold scale factor to a specified value\n");
                printf("  -s <dir>    Set the search directory to a specified value\n");
                printf("  -e <value>  Set the ring exclusion proximity radius\n");
                printf("  -p          Enable histogram through GNUPLOT\n");
                printf("  -g <sigma>  Set gaussian filtering sigma value\n");
                printf("  -n          Display the notes, tips, and comments");
                return(0);
            case 'd':
                cfg.data = optarg;
                break;
            case 'o':
                cfg.output = true;
                break;
            case 'i':
                cfg.index = atoi(optarg);
                break;
            case 'r':
                cfg.traversing =1;
                break;   
            case 't':
                cfg.thresfactor = atoi(optarg);
                break;
            case 's':
                cfg.searchdir =optarg;
                break;
            case 'e':
                cfg.percent =atof(optarg);
                break;
            case 'p':
                cfg.gnu = true;
                break;
            case 'g':
                cfg.gaussian = atof(optarg);
                break;
            case 'f':
                cfg.filename = optarg;
                break;  
            case 'n':
                noteprint();
                return 0;
            case '?':
                if (optopt == 'f' || optopt == 'i' || optopt == 't' || optopt == 's' || optopt == 'd') {
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                } else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                } else {
                    printf("Help: use -f for file input, -o for output\n");

                }
                return 1;
            default:
                abort();
        }
    }
    if (!cfg.traversing) {
        screener(cfg.filename);
    } else {
        char cwd[PATH_MAX];
        getcwd(cwd, sizeof(cwd));
        int count;
        char** directories = getdirs(cwd, cfg.searchdir, &count);
        if (directories != NULL) {
            printf("Traversal Found %d Directories:\n", count);
            for (int i = 0; i < count; i++) {
                pid_t pid = fork();
                if (pid < 0) {
                    fprintf(stderr, "Fork failed\n");
                    exit(1);
                } else if (pid == 0) {
                    getfile(directories[i]);
                    free(directories[i]);
                    exit(0);
                }
            }
            int status;
            pid_t child_pid;
            while ((child_pid = wait(&status)) > 0);

            for (int i = 0; i < count; i++) {
                free(directories[i]);
            }
            free(directories);
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    argparse(argc, argv);
    return 0;
}