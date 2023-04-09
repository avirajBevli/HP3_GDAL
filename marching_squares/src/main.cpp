#include <stdio.h>
#include <stdlib.h>

const int N = 500, M = 500;
const int CSZ = 3;
const int ISOVALUE = 1200;

void readDem(float data[N][M])
{
    FILE *fptr;
    fptr = fopen("../data/dem_file.txt", "r");
    if(fptr == NULL)
    {
        printf("Error reading file!\n");
        exit(1);
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            fscanf(fptr, "%f", &data[i][j]);

    fclose(fptr);
}

void applyThresh(int binData[N][M], float data[N][M], int isovalue = ISOVALUE)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            binData[i][j] = (data[i][j] > isovalue);
        }
    }
    
}

int makeContourUtil(int a, int b, int c, int d)
{
    return (a * 8) + (b * 4) + (c * 2) + (d * 1);
}

void makeContourGrid(int contourGrid[N - 1][M - 1], int binData[N][M])
{
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < M - 1; j++)
        {
            contourGrid[i][j] = makeContourUtil(binData[i][j], 
                                                binData[i][j + 1],
                                                binData[i + 1][j + 1],
                                                binData[i + 1][j]);
        }
    } 
}


void drawContour(int contourPlot[CSZ][CSZ], int contourNum)
{
    for (int i = 0; i < CSZ; i++)
    {
        for (int j = 0; j < CSZ; j++)
        {
            contourPlot[i][j] = 0;
        }
        
    }
    
    if(contourNum == 0 || contourNum == 15)
    {
        return;
    }
    else if(contourNum == 1 || contourNum == 14)
    {
        for(int i = CSZ/2; i < CSZ; i++)
            contourPlot[i][i - CSZ/2] = 1;
    }
    else if(contourNum == 2 || contourNum == 13)
    {
        for(int i = CSZ/2; i < CSZ; i++)
            contourPlot[i][CSZ - i + 1] = 1;
        
    }
    else if(contourNum == 3 || contourNum == 12)
    {
        for (int i = 0; i < CSZ; i++)
            contourPlot[CSZ / 2][i] = 1;
    }
    else if(contourNum == 4 || contourNum == 11)
    {
        for (int i = 0; i <= CSZ/2; i++)
            contourPlot[i][CSZ/2 + i] = 1;
    }
    else if(contourNum == 6 || contourNum == 9)
    {
        for(int i = 0; i < CSZ; i++)
            contourPlot[i][CSZ/2] = 1;
    }
    else if(contourNum == 7 || contourNum == 8)
    {
        for(int i = 0; i <= CSZ / 2; i++)
            contourPlot[i][CSZ/2 - i] = 1;
    }
    else if(contourNum == 5)
    {
        for(int i = 0; i <= CSZ / 2; i++)
        {
            contourPlot[i][CSZ/2 - i] = 1;
            contourPlot[i + CSZ/2][CSZ - i - 1] = 1;
        }
    }
    else if(contourNum == 10)
    {
        for (int i = 0; i <= CSZ/2; i++)
        {
            contourPlot[i][CSZ/2 + i] = 1;
            contourPlot[i + CSZ/2][i] = 1;
        }
    }
    else
    {
        printf("Invalid contour number\n");
        exit(1);
    }
    
}

void drawContours(short int contourPlot[(N - 1) * CSZ][(M - 1) * CSZ], int contourGrid[N - 1][M - 1])
{
    int contour[CSZ][CSZ];

    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < M - 1; j++)
        {
            drawContour(contour, contourGrid[i][j]);
            for (int row = 0; row < CSZ; row++)
                for (int col = 0; col < CSZ; col++)
                    contourPlot[i * CSZ + row][j * CSZ + col] = contour[row][col];
        }
    }
}

void writeContourFile(short int contourPlot[(N - 1) * CSZ][(M - 1) * CSZ])
{
    FILE *fptr;
    fptr = fopen("contour_plot.txt", "w");
    for (int i = 0; i < (N - 1) * CSZ; i++)
    {
        for (int j = 0; j < (M - 1) * CSZ; j++)
        {
            fprintf(fptr, "%d", contourPlot[i][j]);
            if(j != (M - 1) * CSZ - 1)
                fprintf(fptr, " ");

        }
        fprintf(fptr, "\n");
    }
}

int main()
{
    float data[N][M];
    readDem(data);

    int binData[N][M];
    applyThresh(binData, data);
    
    int contourGrid[N - 1][M - 1];
    makeContourGrid(contourGrid, binData);
    
    // printf("DEM:\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < M; j++)
    //     {
    //         printf("%f (%d) ", data[i][j], binData[i][j]);
    //     }
    //     printf("\n");
    // }

    // printf("\n Contour Grid:\n");
    // for (int i = 0; i < N - 1; i++)
    // {
    //     for (int j = 0; j < M - 1; j++)
    //     {
    //         printf("%d ", contourGrid[i][j]);
    //     }
    //     printf("\n");
    // }

    short int contourPlot[(N - 1) * CSZ][(M - 1) * CSZ];
    drawContours(contourPlot, contourGrid);

    // // for (int i = 0; i < (N - 1) * CSZ; i++)
    // // {
    // //     for (int j = 0; j < (M - 1) * CSZ; j++)
    // //     {
    // //         printf("%d ", contourPlot[i][j]);
    // //     }
    // //     printf("\n");
        
    // // }

    writeContourFile(contourPlot);
    

    return 0;
}