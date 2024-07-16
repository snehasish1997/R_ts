#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<limits.h>
#include "coordReader.c"


int readNumOfCoords(char *filename);
double** readCoords(char *filename, int numOfCoords);
double** createDistanceMatrix(double **coords, int numOfCoords);
void smallestSumInsertion(double **distMatrix, int numOfCoords, int *tour);
void *writeTourToFile(int *tour, int tourLength, char *filename);


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
    smallestSumInsertion(distMatrix, numOfCoords, tour);
    
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
 

double calculateDistance(double *coord1, double *coord2) {
    return sqrt(pow(coord1[0] - coord2[0], 2) + pow(coord1[1] - coord2[1], 2));
}

double** createDistanceMatrix(double **coords, int numOfCoords) {
    double **distMatrix = (double **)malloc(numOfCoords * sizeof(double *));
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
    
    for (int i = 1; i < numOfCoords; i++) {
        int bestVertex = -1;
        double bestCost = INFINITY;
        
        for (int j = 0; j < numOfCoords; j++) {
            if (!visited[j]) {
                double cost = 0;
                for (int k = 0; k < i; k++) {
                    cost += distMatrix[tour[k]][j];
                }
                if (cost < bestCost) {
                    bestCost = cost;
                    bestVertex = j;
                }
            }
        }
        
        int bestInsertPos = -1;
        double bestInsertCost = INFINITY;
        
        for (int k = 0; k < i; k++) {
            double insertCost = distMatrix[tour[k]][bestVertex] + distMatrix[tour[k + 1 == i ? 0 : k + 1]][bestVertex] - distMatrix[tour[k]][tour[k + 1 == i ? 0 : k + 1]];
            if (insertCost < bestInsertCost) {
                bestInsertCost = insertCost;
                bestInsertPos = k;
            }
        }
        
        for (int k = i; k > bestInsertPos + 1; k--) {
            tour[k] = tour[k - 1];
        }
        tour[bestInsertPos + 1] = bestVertex;
        visited[bestVertex] = 1;
    }
    free(visited);
}









/*double euclideanDistance(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

void calculateDistanceMatrix() {
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = 0; j < numOfCoords; j++) {
            if (i == j) {
                distanceMatrix[i][j] = 0.0;
            } else {
                distanceMatrix[i][j] = euclideanDistance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            }
        }
    }
}*/