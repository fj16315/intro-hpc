
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"

// Define output file name
#define OUTPUT_FILE "stencil.pgm"
#define MASTER 0

void stencil(const int nx, const int ny, float *  image, float *  tmp_image, const int rank, const int size);
void init_image(const int nx, const int ny, float *  image);
void init_tmp(const int nx, const int ny, float * tmp_image);
void output_image(const char * file_name, const int nx, const int ny, float *image);
double wtime(void);

int main(int argc, char *argv[]) {

  int rank;               /* 'rank' of process among it's cohort */
  int size;               /* size of cohort, i.e. num processes started */
  int flag;               /* for checking whether MPI_Init() has been called */
  enum bool {FALSE,TRUE}; /* enumerated type: false = 0, true = 1 */
  int tag = 0;           /* scope for adding extra information to a message */
  MPI_Status status;     /* struct used by MPI_Recv */

  MPI_Init(&argc, &argv);

  MPI_Initialized(&flag);
  if (flag != TRUE) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  flag = FALSE;

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int portionHeight = ny / size;
  int masterPortionHeight = ny-(portionHeight*(size-1));
  int niters = atoi(argv[3]);
  float *image;
  float *tmp_image;
  float *imagePortion;

  if (rank == MASTER) {
    // Allocate the image
    image = malloc(sizeof(float)*nx*ny);
    tmp_image = malloc(sizeof(float)*nx*masterPortionHeight);
    imagePortion = malloc(sizeof(float)*nx*masterPortionHeight);
    // Set the input image
    init_image(nx, ny, image);
    init_tmp(nx, masterPortionHeight, tmp_image);
    printf("Master about to scatter\n");
    for (int i = 1; i < size; i++) {
      for (int j = 0; j < nx*portionHeight; j++) {
        imagePortion[j] = image[(nx*masterPortionHeight)+((i-1)*nx*portionHeight)+j];
      }
      MPI_Send(imagePortion, nx*portionHeight, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
    }
    printf("Master scattered\n");
    for (int j = 0; j < nx*masterPortionHeight; j++) {
      imagePortion[j] = image[j];
    }
  }
  else {
    imagePortion = malloc(sizeof(float)*nx*portionHeight);
    tmp_image = malloc(sizeof(float)*nx*portionHeight);
    init_tmp(nx, portionHeight, tmp_image);
    MPI_Recv(imagePortion, nx*portionHeight, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double tic = wtime();
  if (rank == MASTER) {
    printf("Master running stencil!\n");
    for (int t = 0; t < niters; ++t) {
      stencil(nx, masterPortionHeight, imagePortion, tmp_image, rank, size);
      stencil(nx, masterPortionHeight, tmp_image, imagePortion, rank, size);
    }
  }
  else {
    for (int t = 0; t < niters; ++t) {
      stencil(nx, portionHeight, imagePortion, tmp_image, rank, size);
      stencil(nx, portionHeight, tmp_image, imagePortion, rank, size);
    }
  }
  double toc = wtime();
  double finalTime = toc-tic;

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank != MASTER) {
    MPI_Send(imagePortion, nx*portionHeight, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD);
    MPI_Send(&finalTime, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
  }
  else {
    double maxTime = finalTime;
    printf("Master ran with time: %lf s\n", finalTime);
    for (int j = 0; j < nx*masterPortionHeight; j++) {
      image[j] = imagePortion[j];
    }
    printf("Before gather\n");

    for (int i = 1; i < size; i++) {
      MPI_Recv(imagePortion, nx*portionHeight, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&finalTime, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      printf("Rank %d ran with time: %lf s\n", i, finalTime);
      if (finalTime > maxTime) {
        maxTime = finalTime;
      }
      for (int j = 0; j < nx*portionHeight; j++) {
        image[(nx*masterPortionHeight)+((i-1)*nx*portionHeight)+j] = imagePortion[j];
      }
    }

    printf("After gather\n");
    // Output
    printf("------------------------------------\n");
    printf(" runtime from rank %d: %lf s\n", rank, maxTime);
    printf("------------------------------------\n");

    output_image(OUTPUT_FILE, nx, ny, image);
    free(image);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}

void stencil(const int nx, const int ny, float *  restrict image, float *  restrict tmp_image, const int rank, const int size) {
  float *currentTopLine = malloc(sizeof(float)*nx);
  float *currentBottomLine = malloc(sizeof(float)*nx);
  float *bottomTopLine = malloc(sizeof(float)*nx);
  float *topBottomLine = malloc(sizeof(float)*nx);
  MPI_Status status;
  MPI_Request req;
  int flag = 0;

  if (rank != MASTER && rank != size-1) {
    for (int i = 0; i < nx; i++) {
      currentTopLine[i] = image[i];
      currentBottomLine[i] = image[((ny-1)*nx)+i];
    }
  }
  else if (rank == MASTER) {
    for (int i = 0; i < nx; i++) {
      currentBottomLine[i] = image[((ny-1)*nx)+i];
    }
  }
  else {
    for (int i = 0; i < nx; i++) {
      currentTopLine[i] = image[i];
    }
  }

  if (rank != size-1 && size > 1) {
    MPI_Isend(currentBottomLine, nx, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &req);
  }
  if (rank != MASTER) {
    MPI_Irecv(topBottomLine, nx, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &req);
  }
  while (!flag && size > 1) {
    MPI_Test(&req, &flag, &status);
  }
  flag = 0;

  if (rank != MASTER) {
    MPI_Isend(currentTopLine, nx, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &req);
  }
  if (rank != size-1 && size > 1) {
    MPI_Irecv(bottomTopLine, nx, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &req);
  }
  while (!flag && size > 1) {
    MPI_Test(&req, &flag, &status);
  }
  flag = 0;
  //First rank top and bottom lines and corners
  if (rank == MASTER) {
    for (int i = 1; i < nx-1; i++) {
      tmp_image[i] = image[i] * 0.6f; //Weight current pixel
      tmp_image[i] += (image[i-1] + image[i+1] + image[i+nx]) * 0.1f;
      tmp_image[((ny-1)*nx)+i] = image[((ny-1)*nx)+i] * 0.6f;
      tmp_image[((ny-1)*nx)+i] += (image[((ny-1)*nx)+i-1] + image[((ny-1)*nx)+i+1] + image[((ny-2)*nx)+i] + bottomTopLine[i]) * 0.1f;
    }
    tmp_image[0] = image[0] * 0.6f;
    tmp_image[0] += (image[1] + image[nx]) * 0.1f;
    tmp_image[nx-1] = image[nx-1] * 0.6f;
    tmp_image[nx-1] += (image[nx-2] + image[(2*nx)-1]) * 0.1f;
    tmp_image[(ny-1)*nx] = image[(ny-1)*nx] * 0.6f;
    tmp_image[(ny-1)*nx] += (image[(ny-2)*nx] + image[((ny-1)*nx)+1] + bottomTopLine[0]) * 0.1f;
    tmp_image[(ny*nx)-1] = image[(ny*nx)-1] * 0.6f;
    tmp_image[(ny*nx)-1] += (image[((ny-1)*nx)-1] + image[(ny*nx)-2] + bottomTopLine[nx-1]) * 0.1f;
  }
  //Last rank top and bottom lines and corners
  else if (rank == size-1) {
    for (int i = 1; i < nx-1; i++) {
      tmp_image[i] = image[i] * 0.6f; //Weight current pixel
      tmp_image[i] += (image[i-1] + image[i+1] + image[i+nx] + topBottomLine[i]) * 0.1f;
      tmp_image[((ny-1)*nx)+i] = image[((ny-1)*nx)+i] * 0.6f; //Weight current pixel
      tmp_image[((ny-1)*nx)+i] += (image[((ny-1)*nx)+i-1] + image[((ny-1)*nx)+i+1] + image[((ny-2)*nx)+i]) * 0.1f;
    }
    tmp_image[0] = image[0] * 0.6f;
    tmp_image[0] += (image[1] + image[nx] + topBottomLine[0]) * 0.1f;
    tmp_image[nx-1] = image[nx-1] * 0.6f;
    tmp_image[nx-1] += (image[nx-2] + image[(2*nx)-1] + topBottomLine[nx-1]) * 0.1f;
    tmp_image[(ny-1)*nx] = image[(ny-1)*nx] * 0.6f;
    tmp_image[(ny-1)*nx] += (image[(ny-2)*nx] + image[((ny-1)*nx)+1]) * 0.1f;
    tmp_image[(ny*nx)-1] = image[(ny*nx)-1] * 0.6f;
    tmp_image[(ny*nx)-1] += (image[((ny-1)*nx)-1] + image[(ny*nx)-2]) * 0.1f;
  }
  //Other ranks top and bottom lines and corners
  else {
    for (int i = 1; i < nx-1; i++) {
      tmp_image[i] = image[i] * 0.6f; //Weight current pixel
      tmp_image[i] += (image[i-1] + image[i+1] + image[i+nx] + topBottomLine[i]) * 0.1f;
      tmp_image[((ny-1)*nx)+i] = image[((ny-1)*nx)+i] * 0.6f; //Weight current pixel
      tmp_image[((ny-1)*nx)+i] += (image[((ny-1)*nx)+i-1] + image[((ny-1)*nx)+i+1] + image[((ny-2)*nx)+i] + bottomTopLine[i]) * 0.1f;
    }
    tmp_image[0] = image[0] * 0.6f;
    tmp_image[0] += (image[1] + image[nx] + topBottomLine[0]) * 0.1f;
    tmp_image[nx-1] = image[nx-1] * 0.6f;
    tmp_image[nx-1] += (image[nx-2] + image[(2*nx)-1] + topBottomLine[nx-1]) * 0.1f;
    tmp_image[(ny-1)*nx] = image[(ny-1)*nx] * 0.6f;
    tmp_image[(ny-1)*nx] += (image[(ny-2)*nx] + image[((ny-1)*nx)+1] + bottomTopLine[0]) * 0.1f;
    tmp_image[(ny*nx)-1] = image[(ny*nx)-1] * 0.6f;
    tmp_image[(ny*nx)-1] += (image[((ny-1)*nx)-1] + image[(ny*nx)-2] + bottomTopLine[nx-1]) * 0.1f;
  }

  //#pragma omp simd
  //Left and rightmost sides
  for (int j = 1; j < ny-1; j++) {
    tmp_image[j*nx] = image[j*nx] * 0.6f; //Weight current pixel
    tmp_image[j*nx] += (image[(j-1)*nx] + image[(j+1)*nx] + image[(j*nx)+1]) * 0.1f;
    tmp_image[((j+1)*nx)-1] = image[((j+1)*nx)-1] * 0.6f;
    tmp_image[((j+1)*nx)-1] += (image[(j*nx)-1] + image[((j+2)*nx)-1] + image[((j+1)*nx)-2]) * 0.1f;
  }

  //#pragma omp simd
  //Other pixels
  for (int i = 1; i < ny-1; ++i) {
    for (int j = 1; j < nx-1; ++j) { //Image stored as arrayof Doubles, column by column
      tmp_image[j+i*nx] = image[j+i*nx] * 0.6f; //Weight current pixel
      tmp_image[j+i*nx] += (image[j  +(i-1)*nx] + image[j  +(i+1)*nx] + image[j-1+i*nx] + image[j+1+i*nx]) * 0.1f;
    }
  }
}

// Create the input image
void init_image(const int nx, const int ny, float *  image) {
  // Zero everything
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      image[0+j+i*ny] = 0.0f;
    }
  }

  // Checkerboard
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int ii = (i*nx/8); ii < ((i+1)*nx/8); ++ii) {
        for (int jj = (j*ny/8); jj < ((j+1)*ny/8); ++jj) {
          if ((i+j)%2)
          image[jj+ii*ny] = 100.0f;
        }
      }
    }
  }
}

void init_tmp(const int nx, const int ny, float * tmp_image) {
  // Zero everything
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      tmp_image[j+i*nx] = 0.0f;
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const int nx, const int ny, float *image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  double maximum = 0.0;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      fputc((char)(255.0*image[j+i*ny]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
