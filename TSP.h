#ifndef TSP_H
#define TSP_H

//  Define Constants: Size of chromosome population, number of cities and size of elite
#define POPSIZE       40000
#define NUM_CITIES    250
#define ELITSIZE      ((int)(POPSIZE * 0.1f))

// Variable used to generate pseudo-random numbers
unsigned int seed;

// Type definition of a Point
typedef struct
{
    int x;
    int y;
} Point;

// Structure base of chromosome. A chromosome is interpreted as a path
typedef struct
{
    float distance;
    int   gen[NUM_CITIES];
} Chrom;

// Array containing the positions of the cities
Point Cities[NUM_CITIES];

// Helping structure to account for cities that has been already visited
int MASK[NUM_CITIES];

Chrom population[POPSIZE];
Chrom tmp_population[POPSIZE];

#endif
