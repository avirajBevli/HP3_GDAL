#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <map>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <cmath>
#include "device_launch_parameters.h"
#include "helper_cuda.h"
#include "kernel.h"
#include "npp.h"
#include "EasyBMP.h"
#include "EasyBMP_DataStructures.h"
#include "EasyBMP_BMP.h"


using namespace std;

//uble PCFreq = 0.0;

int radiusGlob = 0;
int currenIteration = 0;
int iterationCount = 1;
ofstream resultFile;

#define MEMORYMETRICS
#define BLOCK_DIM 512


__global__ void cudaR2_VSA(vs_t* viewshed, elev_t* elev, elev_t observer_elev, int minX, int maxX, int minY, int maxY, int observerX, int observerY, int ncols)
{

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int width = (maxX - minX) + 1;
	int height = (maxY - minY) + 1;

	int totalCell = ((width + height) * 2) - 4;

	if (idx >= totalCell)
	{
		return;
	}

	int x, y;
	if (idx < width)
	{
		x = minX + idx;
		y = minY;
	}
	else if (idx >= width && idx < width + height - 1)
	{
		x = maxX;
		y = minY + (idx + 1) - width;
	}
	else if (idx >= width + height - 1 && idx < width + height + width - 2)
	{
		x = maxX - ((idx + 1) - (width + height - 1));
		y = maxY;
	}
	else if (idx >= width + height + width - 2 && idx < totalCell)
	{
		x = minX;
		y = maxY - ((idx + 1) - (width + height + width - 2));
	}


	int x1 = observerX;
	int y1 = observerY;
	int x2 = x;
	int y2 = y;

	int delta_x(x2 - x1);
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;


	float maxGradient = -10000;

	if (delta_x >= delta_y)
	{
		int error(delta_y - (delta_x >> 1));

		while (x1 != x2)
		{
			if ((error >= 0) && (error || (ix > 0)))
			{
				error -= delta_x;
				y1 += iy;
			}
			error += delta_y;
			x1 += ix;


			int currentIndex = (x1 + y1 * ncols);
			int deltaY = y1 - observerY;
			int deltaX = x1 - observerX;
			float dist2 = deltaX * deltaX + deltaY * deltaY;

			double diff_elev = elev[currentIndex] - observer_elev;
			float gradient = (diff_elev * diff_elev) / dist2;
			if (diff_elev < 0) gradient *= -1;

			if (gradient > maxGradient)
			{
				maxGradient = gradient;
				viewshed[currentIndex] = 1;
			}
			else
			{
				viewshed[currentIndex] = 0;
			}
		}
	}
	else
	{
		int error(delta_x - (delta_y >> 1));

		while (y1 != y2)
		{
			if ((error >= 0) && (error || (iy > 0)))
			{
				error -= delta_y;
				x1 += ix;
			}
			error += delta_x;
			y1 += iy;

			int currentIndex = (x1 + y1 * ncols);

			int deltaY = y1 - observerY;
			int deltaX = x1 - observerX;
			float dist2 = deltaX * deltaX + deltaY * deltaY;

			double diff_elev = elev[currentIndex] - observer_elev;
			float gradient = (diff_elev * diff_elev) / dist2;
			if (diff_elev < 0) gradient *= -1;

			if (gradient > maxGradient)
			{
				maxGradient = gradient;
				viewshed[currentIndex] = 1;
			}
			else
			{
				viewshed[currentIndex] = 0;
			}

		}
	}
}

void cudaR2_caller(vs_t* viewshed, elev_t* elev, elev_t observer_elev, int minX, int maxX, int minY, int maxY, int observerX, int observerY, int ncols)
{
	int width = (maxX - minX) + 1;
	int height = (maxY - minY) + 1;
	int totalCell = ((width + height) * 2) - 4;
	int size = (int)std::ceil((float)((float)totalCell / (float)1024));
	cudaR2_VSA<<< size, 1024 >>> (viewshed, elev, observer_elev, minX, maxX, minY, maxY, observerX, observerY, ncols);
}

unsigned long int SIZECONV = 1024 * 1024;

void write_to_picture(vs_t** viewshed, int ncols, int nrows, char* fileName)
{
	RGBApixel visibleColor;
	visibleColor.Red = 255;
	visibleColor.Green = 255;
	visibleColor.Blue = 255;     //color: white

	RGBApixel notVisibleColor;
	notVisibleColor.Red = 0;
	notVisibleColor.Green = 0;
	notVisibleColor.Blue = 0;    //color: black
	BMP img;
	img.SetSize(ncols, nrows);
	
	for (int i = 0; i < nrows; i++)
	{
		for (int j = 0; j < ncols; j++)
		{
			if (viewshed[i][j] == 1)
			{
				img.SetPixel(j, i, visibleColor);
			}
			else
			{
				img.SetPixel(j, i, notVisibleColor);
			}
		}
	}
	img.WriteToFile(fileName);
}


void write_to_csv(vs_t** viewshed, int ncols, int nrows, char* fileName)
{
	std::ofstream myfile;
	myfile.open(fileName);
	for (int i = 0; i < nrows; i++)
	{

		for (int j = 0; j < ncols; j++)
		{
			if (viewshed[i][j] == 1)
			{
				myfile << "1,";
			}
			else
			{
				myfile << "0,";
			}
		}
		myfile << "\n";
	}
	myfile.close();
}



void write_to_picture(vs_t* viewshed, int ncols, int nrows, char* fileName)
{
	RGBApixel visibleColor;
	visibleColor.Red = 255;
	visibleColor.Green = 255;
	visibleColor.Blue = 255;  //color: white

	RGBApixel notVisibleColor;
	notVisibleColor.Red = 0;
	notVisibleColor.Green = 0;
	notVisibleColor.Blue = 0;  //color: black
	BMP img;
	img.SetSize(ncols, nrows);
	for (int i = 0; i < nrows; i++)
	{
		for (int j = 0; j < ncols; j++)
		{
			int index = i * ncols + j;
			if (viewshed[index] == 1)
			{
				img.SetPixel(j, i, visibleColor);
			}
			else
			{
				img.SetPixel(j, i, notVisibleColor);
			}
		}
	}
	img.WriteToFile(fileName);
}
void write_to_csv(vs_t* viewshed, int ncols, int nrows, char* fileName) {
	std::ofstream myfile;
	myfile.open(fileName);

	for (int i = 0; i < nrows; i++)
	{
		for (int j = 0; j < ncols; j++)
		{
			int index = i * ncols + j;
			if (viewshed[index] == 1)
			{
				myfile << "1,";
			}
			else
			{
				myfile << "0,";
			}
		}
		myfile << "\n";
	}
	myfile.close();

}

void draw_lines(int x, int y, int observerX, int observerY, elev_t** elev, vs_t** viewshed, elev_t observer_elev)
{
	int x1 = observerX;
	int y1 = observerY;
	int x2 = x;
	int y2 = y;

	int delta_x(x2 - x1);
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;


	float maxGradient = -10000;

	if (delta_x >= delta_y)
	{
		int error(delta_y - (delta_x >> 1));

		while (x1 != x2)
		{
			if ((error >= 0) && (error || (ix > 0)))
			{
				error -= delta_x;
				y1 += iy;
			}
			error += delta_y;
			x1 += ix;

			int deltaY = y1 - observerY;
			int deltaX = x1 - observerX;
			float dist2 = deltaX * deltaX + deltaY * deltaY;

			double diff_elev = elev[y1][x1] - observer_elev;
			float gradient = (diff_elev * diff_elev) / dist2;
			if (diff_elev < 0) gradient *= -1;

			if (gradient > maxGradient)
			{
				maxGradient = gradient;
				viewshed[y1][x1] = 1;
			}
			else
			{
				viewshed[y1][x1] = 0;
			}
		}
	}
	else
	{
		int error(delta_x - (delta_y >> 1));
		while (y1 != y2)
		{
			if ((error >= 0) && (error || (iy > 0)))
			{
				error -= delta_y;
				x1 += ix;
			}
			error += delta_x;
			y1 += iy;

			int deltaY = y1 - observerY;
			int deltaX = x1 - observerX;
			float dist2 = deltaX * deltaX + deltaY * deltaY;

			double diff_elev = elev[y1][x1] - observer_elev;
			float gradient = (diff_elev * diff_elev) / dist2;
			if (diff_elev < 0) gradient *= -1;

			if (gradient > maxGradient)
			{
				maxGradient = gradient;
				viewshed[y1][x1] = 1;
			}
			else
			{
				viewshed[y1][x1] = 0;
			}
		}
	}
}


void R2_algo_CPU(int nrows, int ncols, int radius, int observerX, int observerY, int observer_ht, elev_t** elev)
{

	vs_t** viewshed;
	elev_t observer_elev;

	observer_elev = elev[observerY][observerX] + observer_ht;

	//area: Point of interest
	int minY = max(0, observerY - radius);
	int maxY = min(nrows - 1, observerY + radius);
	int minX = max(0, observerX - radius);
	int maxX = min(ncols - 1, observerX + radius);


	int width = (maxX - minX) + 1;
	int height = (maxY - minY) + 1;

	int x = minX;
	int y = minY;
	
	viewshed = new vs_t * [nrows];
	for (int i = 0; i < nrows; i++)
		viewshed[i] = new vs_t[ncols];

	viewshed[observerY][observerX] = 1;

	for (int y = maxY, x = maxX; y >= minY; y--) { 
		draw_lines(x, y, observerX, observerY, elev, viewshed, observer_elev);
	}
	for (int x = maxX, y = minY; x >= minX; x--) { 
		draw_lines(x, y, observerX, observerY, elev, viewshed, observer_elev);
	}
	for (int y = minY, x = minX; y <= maxY; y++) { 
		draw_lines(x, y, observerX, observerY, elev, viewshed, observer_elev);
	}
	for (int x = minX, y = maxY; x <= maxX; x++) { 
		draw_lines(x, y, observerX, observerY, elev, viewshed, observer_elev);
	}

	if (currenIteration == iterationCount)
	{
		char buffer[100];
		sprintf(buffer, "CPU-R2-%d-%d.bmp", observerX, observerY);
		write_to_picture(viewshed, ncols, nrows, buffer);
	}
	if (currenIteration == iterationCount)
	{
		char buffer[100];
		sprintf(buffer, "CPU-R2-%d-%d.csv", observerX, observerY);
		write_to_csv(viewshed, ncols, nrows, buffer);

	}

	for (int i = 0; i < nrows; i++) {
		delete[] viewshed[i];
	}

	delete[] viewshed;
}


double R2_algo_GPU(int nrows, int ncols, int radius, int observerX, int observerY, int observer_ht, elev_t* elev)
{
	elev_t observer_elev;
	observer_elev = elev[observerY * ncols + observerX] + observer_ht;

	cudaEvent_t startEvent, stopEvent;
	cudaEventCreate(&startEvent);
	cudaEventCreate(&stopEvent);

	//area: Point of Interest
	int minY = max(0, observerY - radius);
	int maxY = min(nrows - 1, observerY + radius);
	int minX = max(0, observerX - radius);

	cudaEventRecord(startEvent);

	vs_t* viewshed = new vs_t[nrows*ncols];
	vs_t* d_viewshed;
	checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_viewshed), sizeof(vs_t) * nrows * ncols));
	elev_t* d_elev;
	checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_elev), sizeof(elev_t) * nrows * ncols));
	checkCudaErrors(cudaMemcpy(d_elev, elev, sizeof(elev_t) * nrows * ncols, cudaMemcpyHostToDevice));
	cudaR2_caller(d_viewshed, d_elev, observer_elev, minX, maxX, minY, maxY, observerX, observerY, ncols);
	checkCudaErrors(cudaMemcpy(viewshed, d_viewshed, sizeof(vs_t) * nrows * ncols, cudaMemcpyDeviceToHost));
	
	cudaEventRecord(stopEvent);
	cudaEventSynchronize(stopEvent);
	
	float result = 0.0;

	cudaEventElapsedTime(&result, startEvent, stopEvent);


	if (currenIteration == iterationCount)
	{
		char buffer[100];
		sprintf(buffer, "GPU-R2-%d-%d.bmp", observerX, observerY);
		write_to_picture(viewshed, ncols, nrows, buffer);
	}
	if (currenIteration == iterationCount)
	{
		char buffer[100];
		sprintf(buffer, "GPU-R2-%d-%d.csv", observerX, observerY);
		write_to_csv(viewshed, ncols, nrows, buffer);

	}

	delete[] viewshed;
	checkCudaErrors(cudaFree(d_elev));
	checkCudaErrors(cudaFree(d_viewshed));

	return result;
}

int main(int argc, char** argv)
{
	int observer_ht=15;//observer height
	int nrowCounts;   
	int nColumnCounts;
	int radius;
	char elevPaths[20] = "heights1.txt";
	nrowCounts = 3600; nColumnCounts = 3600;
	int n = 10;
	double result;
	cout << "\n";
	cout << "Enter viewshed analysis radius ";
	cin >> radius;
	//dius = 1000;
	int observerXs[1] = {};     
	int observerYs[1] = {};
	int k = 0;
	cout << "Enter Cordinate of Observer ";
	cin >> observerXs[0] >> observerYs[0];
	//serverXs[0] = 2000; observerYs[0] =2000;
	float result_GPUR2 = 0;

	resultFile.open("results.txt");
	for (int i = 0; i < 1; i++)
	{
		radiusGlob = radius;
		elev_t** elev;
		//Read elevation
		FILE* f = fopen(elevPaths, "rb");
		elev = new elev_t * [nrowCounts];
		for (int k = 0; k < nrowCounts; k++) {
			elev[k] = new elev_t[nColumnCounts];
			fread(reinterpret_cast<char*>(elev[k]), sizeof(elev_t), nColumnCounts, f);
		}
		fclose(f);
		elev_t* elev1D = new elev_t[nrowCounts * nColumnCounts];
		//Read elevation
		f = fopen(elevPaths, "rb");
		int readCount = 0;
		for (int k = 0; k < nrowCounts; k++) {
			fread(reinterpret_cast<char*>(elev1D + readCount), sizeof(elev_t), nColumnCounts, f);
			readCount += nColumnCounts;
		}

		fclose(f);
		result_GPUR2 = 0;

		event_t* d_events;
		checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_events), sizeof(event_t) * 3601 * 3601 * 3));

		for (int j = 0; j < iterationCount; j++)
		{
			currenIteration = j + 1;
			float result;
			R2_algo_CPU(nrowCounts, nColumnCounts, radius, observerXs[i], observerYs[i], observer_ht, elev);
			
			std::cout << "Running time:GPU R2 with Raster Image " <<  std::endl;
			resultFile << "Running time: GPU R2 with Raster Image " << std::endl;
			result = R2_algo_GPU(nrowCounts, nColumnCounts, radius, observerXs[0], observerYs[0], observer_ht, elev1D);
			std::cout << "GPU R2 running time is " << result << std::endl << std::endl << std::endl;
			resultFile << "GPU R2 running time is " << result << std::endl << std::endl << std::endl;
			result_GPUR2 += result;
		}

		delete elev1D;

		for (int k = 0; k < nrowCounts; k++) {
			delete[] elev[k];
		}
		delete[] elev;

	}

	resultFile.close();
	std::cout << "Finished" << endl;
	cudaDeviceReset();
	return 0;
}
