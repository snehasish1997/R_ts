#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "coordReader.c"
#include "ssInsertion.c"

#include "/Users/snehasishbhowmick/.zshrc"

double calculateDistance(double *coord1, double *coord2) {
    return sqrt(pow(coord1[0] - coord2[0], 2) + pow(coord1[1] - coord2[1], 2));
}

double** createDistanceMatrix(double **coords, int numOfCoords) {
    double **distMatrix = (double **)malloc(numOfCoords * sizeof(double *));
    #pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        distMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
        for (int j = 0; j < numOfCoords; j++) {
            distMatrix[i][j] = calculateDistance(coords[i], coords[j]);
        }
    }
    return distMatrix;
}

void smallestSumInsertion(double **distMatrix, int numOfCoords, int *tour) {
    int *visited = (int *)calloc(numOfCoords, sizeof(int));
    tour[0] = 0;
    visited[0] = 1;
    int tourLength = 1;
    
    for (int i = 1; i < numOfCoords; i++) {
        int bestVertex = -1;
        double bestCost = INFINITY;
        
        #pragma omp parallel for
        for (int j = 0; j < numOfCoords; j++) {
            if (!visited[j]) {
                double cost = 0;
                for (int k = 0; k < tourLength; k++) {
                    cost += distMatrix[tour[k]][j];
                }
                #pragma omp critical
                {
                    if (cost < bestCost) {
                        bestCost = cost;
                        bestVertex = j;
                    }
                }
            }
        }
        
        int bestInsertPos = -1;
        double bestInsertCost = INFINITY;
        
        for (int k = 0; k < tourLength; k++) {
            double insertCost = distMatrix[tour[k]][bestVertex] + distMatrix[tour[(k + 1) % tourLength]][bestVertex] - distMatrix[tour[k]][tour[(k + 1) % tourLength]];
            if (insertCost < bestInsertCost) {
                bestInsertCost = insertCost;
                bestInsertPos = k;
            }
        }
        
        for (int k = tourLength; k > bestInsertPos + 1; k--) {
            tour[k] = tour[k - 1];
        }
        tour[bestInsertPos + 1] = bestVertex;
        visited[bestVertex] = 1;
        tourLength++;
    }
    free(visited);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <coordinate file> <output file>\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    char *coordFile = argv[1];
    char *outputFile = argv[2];
    
    int numOfCoords = readNumOfCoords(coordFile);
    double **coords = readCoords(coordFile, numOfCoords);
    double **distMatrix = createDistanceMatrix(coords, numOfCoords);
    
    int *tour = (int *)malloc(numOfCoords * sizeof(int));
    double startTime = omp_get_wtime();
    smallestSumInsertion(distMatrix, numOfCoords, tour);
    double endTime = omp_get_wtime();
    printf("Time taken: %f seconds\n", endTime - startTime);
    
    writeTourToFile(tour, numOfCoords, outputFile);
    
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
        free(distMatrix[i]);
    }
    free(coords);
    free(distMatrix);
    free(tour);
    
    return EXIT_SUCCESS;
}
