
#include "gdal.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <errno.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/core/cuda/common.hpp>
#include <opencv2/cudev.hpp>

#include <cuda.h>
#include <cuda_runtime.h>

using namespace cv;

const char* FILENAME = "../data/cdnh44o.tif";  // File Location
const int CSZ = 15;                            // Cell Size
const bool DISPLAY = 0;                        // Display Contours
const bool ANIMATE = 0;                        // Animate Contours
const bool SAVE = 0;                           // Save Contour Image
const float ISOVALUE = 0;                      // Minimum Elevation of contour visible (0 - everything is visible, 50 - ) TO-DO
const float numContours = 1000;                // Num of different contour levels to show in output image

float adfMinMax[2];
int N = 2000, M = 2000;                                                                              
Mat img = Mat::zeros((N - 1) * CSZ, (M - 1) * CSZ, CV_8UC1);

void readDem(float* data)
{
    GDALDatasetH hDataset;
    GDALAllRegister();
    const GDALAccess eAccess = GA_ReadOnly;
    hDataset = GDALOpen(FILENAME, eAccess);
    if (hDataset == NULL)
    {
        printf("Error in loading Raster File\n");
        return;
    }

    GDALRasterBandH hBand;
    hBand = GDALGetRasterBand(hDataset, 1);
    int bGotMin, bGotMax;
    adfMinMax[0] = GDALGetRasterMinimum(hBand, &bGotMin);
    adfMinMax[1] = GDALGetRasterMaximum(hBand, &bGotMax);
    if (!(bGotMin && bGotMax))
        GDALComputeRasterMinMax(hBand, TRUE, (double*)adfMinMax);

    float *pafScanline;
    int nXSize = GDALGetRasterBandXSize(hBand);
    int nYSize = GDALGetRasterBandYSize(hBand);

    if(N > nYSize || M > nXSize)
    {
        N = nYSize;
        M = nXSize;
    }
    pafScanline = (float *)CPLMalloc(sizeof(float) * nXSize);

    for (int i = 0; i < N; i++)
    {
        CPLErr error = GDALRasterIO(hBand, GF_Read, 0, i, nXSize, 1,
                    pafScanline, nXSize, 1, GDT_Float32,
                    0, 0);
        for (int j = 0; j < M; j++)
        {
            data[i * M + j] = pafScanline[j];
            if(i == 0 && j == 0)
            {
                adfMinMax[0] = data[i * M + j];
                adfMinMax[1] = data[i * M + j];
            }
            adfMinMax[0] = std::min(data[i * M + j], adfMinMax[0]);
            adfMinMax[1] = std::max(data[i * M + j], adfMinMax[1]);
        }
    }

    CPLFree(pafScanline);
    GDALClose(hDataset);

}

__global__ void drawContours(short* contourGrid, cv::cudev::PtrStepSz<uchar> img, int N, int M)
{

    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int chunkRow = N / blockDim.y, chunkCol = M / blockDim.x;
    for(int i = row; i < row + chunkRow && i < N - 1; i++)
    {
        for(int j = col; j < col + chunkCol && j < M - 1; j++)
        {
            if(contourGrid[i * (M - 1) + j] == 1 || contourGrid[i * (M - 1) + j] == 14)
            {
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ/2, CSZ), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 2 || contourGrid[i * (M - 1) + j] == 13)
            {
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, CSZ), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 3 || contourGrid[i * (M - 1) + j] == 12)
            {
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 4 || contourGrid[i * (M - 1) + j] == 11)
            {
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, 0), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 6 || contourGrid[i * (M - 1) + j] == 9)
            {
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, 0), Point(CSZ/2, CSZ), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 7 || contourGrid[i * (M - 1) + j] == 8)
            {
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ/2, 0), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 5)
            {
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ/2, 0), Scalar(255), 2, LINE_8);
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, CSZ), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 10)
            {
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, 0), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ/2, CSZ), Scalar(255), 2, LINE_8);
            }
        }
    }
}

__global__ void makeContourGrid(short* contourGrid, const float* data, int N, int M, int isovalue = ISOVALUE)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int chunkRow = N / blockDim.y, chunkCol = M / blockDim.x;
    for(int i = row; i < row + chunkRow && i < N - 1; i++)
        for(int j = col; j < col + chunkCol && j < M - 1; j++)
                contourGrid[i * (M - 1) + j] = ((data[i * M + j] > isovalue) * 8) + ((data[i * M + j + 1] > isovalue) * 4) + ((data[(i + 1) * M + j + 1] > isovalue) * 2) + ((data[(i + 1) * M + j] > isovalue) * 1);
}

int main(int argc, char* argv[])
{

    // loading data in host memory

    float *h_A = (float*)malloc(sizeof(float) * N * M);
    readDem(h_A);
    float stepSz = (adfMinMax[1] - std::max(ISOVALUE, adfMinMax[0])) / numContours;
    printf("Max: %f, Min: %f\n", adfMinMax[1], adfMinMax[0]);
    printf("Step Size: %f\n", stepSz);
    short *contourGrid = (short*)malloc(sizeof(short) * (N - 1) * (M - 1));


    cudaError_t err = cudaSuccess;

    // device memory allocation
    float* d_A = NULL;
    err = cudaMalloc((void**)&d_A, sizeof(float) * N * M);
    if(err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device DEM array (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    short* d_cg = NULL;
    err = cudaMalloc((void**)&d_cg, sizeof(short) * (N - 1) * (M - 1));
    if(err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device Contour Grid (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // uint8_t *d_img = NULL;
    // err = cudaMalloc((void**)&d_img, sizeof(uint8_t) * ((N - 1) * CSZ) * ((M - 1) * CSZ));
    // if(err != cudaSuccess)
    // {
    //     fprintf(stderr, "Failed to allocate device Contour image (error code %s)!\n", cudaGetErrorString(err));
    //     exit(EXIT_FAILURE);
    // }
    cuda::GpuMat d_img;

    // host to device data transfer
    err = cudaMemcpy(h_A, d_A, sizeof(float) * N * M, cudaMemcpyHostToDevice);
    // err = cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);
    if(err != cudaSuccess)
    {
        fprintf(stderr, "Failed to transfer data from h_A to d_A (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(contourGrid, d_cg, sizeof(short) * (N - 1) * (M - 1), cudaMemcpyHostToDevice);
    if(err != cudaSuccess)
    {
        fprintf(stderr, "Failed to transfer data from contourGrid to d_cg (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    d_img.upload(img);
    // err = cudaMemcpy(img, d_img, sizeof(img), cudaMemcpyHostToDevice);
    // if(err != cudaSuccess)
    // {
    //     fprintf(stderr, "Failed to transfer data from img to d_img (error code %s)!\n", cudaGetErrorString(err));
    //     exit(EXIT_FAILURE);
    // }

    // kernel launch
    dim3 block(32, 32);
    dim3 grid(16, 16);

    for (float i = std::max(ISOVALUE, adfMinMax[0]); i < adfMinMax[1]; i+=stepSz)
    {
        makeContourGrid<<<block, grid>>>(d_cg, d_A, N, M, i);
        drawContours<<<block, grid>>>(d_cg, d_img, N, M);
        // if(ANIMATE)
        // {
        //     namedWindow("Contours", 0);
        //     imshow("Contours", img);
        //     waitKey(1);
        // }
    }

    img = d_img;
    free(h_A);
    free(contourGrid);

    if(DISPLAY)
    {
        namedWindow("Contours", 0);
        imshow("Contours", img);
        waitKey(0);
    }

    if(SAVE)
    {
      imwrite("Contour.png", img);
    }

    
    return 0;
}
