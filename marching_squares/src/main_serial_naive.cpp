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

using namespace cv;

const char* FILENAME = "/home/rohit/Acads/HP3/term_proj/HP3_GDAL/marching_squares/data/cdnh44o.tif";  // File Location
const int CSZ = 11;                                                                                   // Cell Size
const bool INPISTXT = false;                                                                          // Input Type: True if input is .txt file
const bool DISPLAY = false;                                                                           // Display Contours

float ISOVALUE = 0, adfMinMax[2];
int N = 1024, M = 1024;                                                                              
Mat img = Mat::zeros((N - 1) * CSZ, (M - 1) * CSZ, CV_8UC1);


void readDem(float* data, bool inpIsTxt = INPISTXT)
{
    if(!inpIsTxt)
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

        GDALDriverH hDriver;
        hDriver = GDALGetDatasetDriver(hDataset);

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

        if(N > nXSize || M > nYSize)
        {
            N = nXSize;
            M = nYSize;
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
    else
    {
        FILE *fptr;
        fptr = fopen("../data/cmu.txt", "r");
        if(fptr == NULL)
        {
            printf("Error reading file!\n");
            exit(1);
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                fscanf(fptr, "%f", &data[i * M + j]);
                if(i == 0 && j == 0)
                {
                    adfMinMax[0] = data[i * M + j];
                    adfMinMax[1] = data[i * M + j];
                }
                adfMinMax[0] = std::min(data[i * M + j], adfMinMax[0]);
                adfMinMax[1] = std::max(data[i * M + j], adfMinMax[1]);
            }
        }

        fclose(fptr);

    }

}

void drawContours(short* contourGrid)
{
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < M - 1; j++)
        {
            if(contourGrid[i * (M - 1) + j] == 1 || contourGrid[i * (M - 1) + j] == 14)
            {
                // contours[1].copyTo(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)));
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ/2, CSZ), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 2 || contourGrid[i * (M - 1) + j] == 13)
            {
                // contours[2].copyTo(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)));
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, CSZ), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 3 || contourGrid[i * (M - 1) + j] == 12)
            {
                // contours[3].copyTo(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)));
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 4 || contourGrid[i * (M - 1) + j] == 11)
            {
                // contours[4].copyTo(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)));
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, 0), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 6 || contourGrid[i * (M - 1) + j] == 9)
            {
                // contours[6].copyTo(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)));
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, 0), Point(CSZ/2, CSZ), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 7 || contourGrid[i * (M - 1) + j] == 8)
            {
                // contours[7].copyTo(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)));
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ/2, 0), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 5)
            {
                // contours[5].copyTo(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)));
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ/2, 0), Scalar(255), 2, LINE_8);
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, CSZ), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
            }
            else if(contourGrid[i * (M - 1) + j] == 10)
            {
                // contours[8].copyTo(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)));

                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(CSZ/2, 0), Point(CSZ, CSZ/2), Scalar(255), 2, LINE_8);
                line(img(Rect(j * CSZ, i * CSZ, CSZ, CSZ)), Point(0, CSZ/2), Point(CSZ/2, CSZ), Scalar(255), 2, LINE_8);
            }
        }
    }

    // namedWindow("Contour Image", 0);
    // imshow("Contour Image", img);
    // waitKey();
}

void applyThresh(bool* binData, float* data, int isovalue = ISOVALUE)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            binData[i * M + j] = (data[i * M + j] > isovalue);
        }
    }
    
}

int makeContourUtil(int a, int b, int c, int d)
{
    return (a * 8) + (b * 4) + (c * 2) + (d * 1);
}

void makeContourGrid(short* contourGrid, bool* binData)
{
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < M - 1; j++)
        {
            contourGrid[i * (M - 1) + j] = makeContourUtil(binData[i * M + j], 
                                                binData[i * M + j + 1],
                                                binData[(i + 1) * M + j + 1],
                                                binData[(i + 1) * M + j]);
        }
    } 
}

void writeContourFile(short *contourPlot)
{
    FILE *fptr;
    fptr = fopen("../data/contour_plot.txt", "w");
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < M - 1; j++)
        {
            fprintf(fptr, "%d", contourPlot[i * (M - 1) + j]);
            if(j != M - 1)
                fprintf(fptr, " ");

        }
        fprintf(fptr, "\n");
    }
}

int main(int argc, char* argv[])
{

    // preprocess();
    
    float *data = (float*)malloc(sizeof(float) * N * M);
    float numContours = 1000;
    readDem(data);
    printf("Max: %f, Min: %f\n", adfMinMax[1], adfMinMax[0]);
    float stepSz = (adfMinMax[1] - std::max(ISOVALUE, adfMinMax[0])) / numContours;
    printf("Step Size: %f\n", stepSz);
    for (float i = std::max(ISOVALUE, adfMinMax[0]); i < adfMinMax[1]; i+=stepSz)
    {
        bool *binData = (bool*)malloc(sizeof(bool) * N * M);
        applyThresh(binData, data, i);
        
        short *contourGrid = (short*)malloc(sizeof(short) * (N - 1) * (M - 1));
        makeContourGrid(contourGrid, binData);
        drawContours(contourGrid);
        
        if(DISPLAY)
        {
            namedWindow("Contours", 0);
            imshow("Contours", img);
            waitKey(1);
        }
    }

    if(DISPLAY)
    {
        namedWindow("Contours", 1);
        imshow("Contours", img);
        waitKey(0);
    }

    
    return 0;
}