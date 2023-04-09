#include "gdal.h"
#include "cpl_conv.h" // for CPLMalloc()

#include <errno.h>

const bool DEBUG = 0;
const bool GET_GEOTRANSFORM = 0;
const char* FILENAME = "/home/rohit/Acads/HP3/term_proj/data/cdnh44o_v3r1/cdnh44o.tif";

struct Point
{
    Point 
    Point(double _x, double _y) : x(_x), y(_y)
    {
    }

    double x;
    double y;
};

void rasterInfo(GDALDatasetH hDataset, GDALDriverH hDriver, double* adfGeoTransform)
{
    if(DEBUG)
    {
        printf("Driver: %s/%s\n",
            GDALGetDriverShortName(hDriver),
            GDALGetDriverLongName(hDriver));
        printf("Size is %dx%dx%d\n",
            GDALGetRasterXSize(hDataset),
            GDALGetRasterYSize(hDataset),
            GDALGetRasterCount(hDataset));
        if (GDALGetProjectionRef(hDataset) != NULL)
            printf("Projection is `%s'\n", GDALGetProjectionRef(hDataset));
    }
    if(GET_GEOTRANSFORM)
    {
        if (GDALGetGeoTransform(hDataset, adfGeoTransform) == CE_None)
        {
            if(DEBUG)
            {
                printf("Origin = (%.6f,%.6f)\n",
                    adfGeoTransform[0], adfGeoTransform[3]);
                printf("Pixel Size = (%.6f,%.6f)\n",
                    adfGeoTransform[1], adfGeoTransform[5]);
            }
        }
    }
}

void rasterBandInfo(GDALRasterBandH hBand)
{
    if(DEBUG)
    {
        int nBlockXSize, nBlockYSize;
        int bGotMin, bGotMax;
        double adfMinMax[2];
        GDALGetBlockSize(hBand, &nBlockXSize, &nBlockYSize);
        printf("Block=%dx%d Type=%s, ColorInterp=%s\n",
            nBlockXSize, nBlockYSize,
            GDALGetDataTypeName(GDALGetRasterDataType(hBand)),
            GDALGetColorInterpretationName(
                GDALGetRasterColorInterpretation(hBand)));
        adfMinMax[0] = GDALGetRasterMinimum(hBand, &bGotMin);
        adfMinMax[1] = GDALGetRasterMaximum(hBand, &bGotMax);
        if (!(bGotMin && bGotMax))
            GDALComputeRasterMinMax(hBand, TRUE, adfMinMax);
        printf("Min=%.3fd, Max=%.3f\n", adfMinMax[0], adfMinMax[1]);
    }
}



int run(const char *pszFilename)
{

    // Loading the Raster Dataset

    GDALDatasetH hDataset;
    GDALAllRegister();
    const GDALAccess eAccess = GA_ReadOnly;
    hDataset = GDALOpen(pszFilename, eAccess);
    if (hDataset == NULL)
    {
        printf("Error in loading Raster File\n");
        return -1;
    }

    printf("Loaded Raster File Successfully\n");

    // Raster Data Information

    GDALDriverH hDriver;
    double adfGeoTransform[6];
    hDriver = GDALGetDatasetDriver(hDataset);
    rasterInfo(hDataset, hDriver, adfGeoTransform);


    // Accessing Band Specific Information

    GDALRasterBandH hBand;
    hBand = GDALGetRasterBand(hDataset, 1);
    rasterBandInfo(hBand);



    // Accessing Values in the raster file

    float *pafScanline;
    int nXSize = GDALGetRasterBandXSize(hBand);
    int nYSize = GDALGetRasterBandYSize(hBand);
    pafScanline = (float *)CPLMalloc(sizeof(float) * nXSize);

    for (int i = 0; i < 10; i++)
    {
        GDALRasterIO(hBand, GF_Read, 0, i, nXSize, 1,
                     pafScanline, nXSize, 1, GDT_Float32,
                     0, 0);
        for (int j = 0; j < 10; j++)
            printf("%f ", pafScanline[j]);

        printf("\n");
    }

    CPLFree(pafScanline);
    GDALClose(hDataset);

    // vector processing

    return 0;
}

void help()
{
    printf("HELP Prompt to aid user when running the program!")
}

int main(int argc, const char *argv[])
{

    if(DEBUG)
        help();
    
    const char *pszFilename;
    if(argc == 1)
        pszFilename = FILENAME;
    else if(argc == 2)
        pszFilename = argv[1];

    run(pszFilename);
    return 0;
}