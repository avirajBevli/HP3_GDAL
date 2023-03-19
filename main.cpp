#include <iostream>
#include "gdal_priv.h"
#include "gdal_alg.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <errno.h>
using namespace std;

int main(int argc, const char* argv[]){
    const char* dem_filename = "/Users/avirajbevli/Desktop/Sem10/HPPP/cuda_term project/data/cdnf43w/cdnf43w.tif";
    GDALDatasetUniquePtr dem_dataset;
    GDALAllRegister();
    const GDALAccess eAccess = GA_ReadOnly;
    dem_dataset = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen( dem_filename, eAccess )));
    if( !dem_dataset ){
        cout<<"DEM dataset loading mei ERROR ho gaya broo!"<<endl;
    }
    GDALRasterBand *band = dem_dataset->GetRasterBand(1);


    GDALDatasetH dem_dataset_simple = GDALOpen(dem_filename, GA_ReadOnly);
    GDALRasterBandH dem_band_simple = GDALGetRasterBand(dem_dataset_simple, 1);



    const char* inputRasterPath = "/Users/avirajbevli/Desktop/Sem10/HPPP/cuda_term project/data/cdnf43w/cdnf43w.tif";

    // Open the input raster dataset
    GDALDataset* input_ds = static_cast<GDALDataset*>(GDALOpen(inputRasterPath, GA_ReadOnly));

    // Get the geotransform of the input raster
    double geoTransform[6];
    dem_dataset->GetGeoTransform(geoTransform);

    // Get the raster size of the input raster
    int rasterXSize = input_ds->GetRasterXSize();
    int rasterYSize = input_ds->GetRasterYSize();

    // Get the data type of the input raster
    GDALDataType dataType = input_ds->GetRasterBand(1)->GetRasterDataType();

    // Define the output raster file path
    const char *pszDstFilename = "/Users/avirajbevli/Desktop/Sem10/HPPP/cuda_term project/viewshed_data/viewshed_file.tif";

    // Create the output raster dataset
    GDALDataset* output_ds = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(
        pszDstFilename, rasterXSize, rasterYSize, 1, dataType, nullptr);

    // Set the geotransform of the output raster
    output_ds->SetGeoTransform(geoTransform);
    char* wkt = NULL;
    OGRSpatialReference proj;
    proj.SetWellKnownGeogCS("WGS84");
    proj.exportToWkt(&wkt); // Export the projection in WKT format
    output_ds->SetProjection(wkt); // Set the projection



    // Open Memory Dataset
    // GDALDatasetH dem_dataset_simple = GDALOpen(memdsetpath, GA_Update);
    // double adfGeoTransform[6];
    // GDALSetGeoTransform(dem_dataset_simple, adfGeoTransform);
    // GDALSetProjection(dem_dataset_simple, pszProjection); // Can not


    double observer_x = 76.4998611;
    double observer_y = 20.5001389;
    double observer_height;
    band->RasterIO(GF_Read, observer_x, observer_y, 1, 1, &observer_height, 1, 1, GDT_Float64, 0, 0);
    // GDALApplyGeoTransform(geoTransform, observer_x_map, observer_x_map, &observer_x, &observer_y);
    cout<<"observer_x:"<<observer_x<<", observer_y:"<<observer_y<<endl;
    cout<<"observer_height: "<<observer_height<<endl;

    GDALProgressFunc pfnProgress = GDALTermProgress;
    void* pProgressArg = NULL;
    char** papszOptions = NULL;
    GDALRasterBand* viewshed_band = output_ds->GetRasterBand(1);
    // GDALDatasetH viewshed_dataset = GDALCreate("viewshed.tif", dem_dataset->GetRasterXSize(), dem_dataset->GetRasterYSize(), 1, eOutputType, papszOptions);
    /*
    GDALDatasetH GDALViewshedGenerate(GDALRasterBandH hBand, const char *pszDriverName, 
    const char *pszTargetRasterName, CSLConstList papszCreationOptions, 
    double dfObserverX, double dfObserverY, double dfObserverHeight, 
    double dfTargetHeight, double dfVisibleVal, double dfInvisibleVal, double dfOutOfRangeVal, 
    double dfNoDataVal, double dfCurvCoeff, GDALViewshedMode eMode, double dfMaxDistance, 
    GDALProgressFunc pfnProgress, void *pProgressArg, GDALViewshedOutputType heightMode, 
    CSLConstList papszExtraOptions)
    'void* GDALViewshedGenerate(GDALRasterBandH, const char*, 
    const char*, CSLConstList, 
    double, double, double, 
    double, double, double, double, 
    double, double, GDALViewshedMode, double, 
    GDALProgressFunc, void*, GDALViewshedOutputType, 
    CSLConstList)
    */
    // GDALDatasetH hDstDS = GDALViewshedGenerate(dem_band_simple, "GTIFF",
    //     pszDstFilename, NULL, 
    //     observer_x, observer_y, observer_height, 
    //     target_height, min_angle, max_angle, max_distance, 
    //     0.0, 0.0, GVM_Edge, 60, 
    //     pfnProgress, pProgressArg, GVOT_NORMAL, 
    //     papszOptions);

    //double dfObserverHeight = 2.0;
    double dfTargetHeight = 200.0;
    double dfMaxDistance = 10.0;
    //double dfObserverX = 0.0;
    //double dfObserverY = 0.0;
    double dfVisibleVal = 255; // pixel value for visibility 
    double dfInvisibleVal = 100; // pixel value for invisibility
    double dfOutOfRangeVal = -1;
    double dfNoDataVal = 10; // 0
    double dfCurvCoeff = 0; // 0.85714

    GDALDatasetH hDstDS = GDALViewshedGenerate(dem_band_simple, "GTIFF",
        pszDstFilename, NULL, 
        observer_x, observer_y, observer_height, 
        dfTargetHeight, dfVisibleVal, dfInvisibleVal, dfOutOfRangeVal, 
        dfNoDataVal, dfCurvCoeff, GVM_Edge, dfMaxDistance, 
        pfnProgress, pProgressArg, GVOT_NORMAL, 
        papszOptions);

    bool bSuccess = hDstDS != nullptr;
    if (bSuccess == FALSE) {
        cout<<"Viewshed computation mei ERROR ho gaya broo!"<<endl;
        return 0;
    }


    // GDALClose(dem_dataset);
    // GDALClose(viewshed_dataset);
    cout<<"Exiting ...."<<endl;
    // Cleanup
    CPLFree(wkt);
    GDALClose(output_ds);
    return 0;
}