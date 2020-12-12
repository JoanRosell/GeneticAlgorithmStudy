/**
 * description:
 *  Genetic algorithm for finding heuristic solution of Travelling Salesman Problem
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

// Variable used to generate pseudo-random numbers
unsigned int seed;

// Type definition of a Point
typedef float coord_t;
typedef struct Point
{
    coord_t x;
    coord_t y;
} Point;

// Structure base of chromosome. A chromosome is interpreted as a path
typedef uint8_t tag_t;
typedef struct Chromosome
{
    float distance;
    tag_t*  tour;
} Chromosome;

// Function to generate pseudo-random numbers
int myRandom()
{
    seed = (214013 * seed + 2531011);
    return(seed >> 13);
}

// Function definitions
void generateCities(Point* cities, size_t nCities);
void initPopulation(Chromosome* population, size_t popSize, size_t nCities);
void mutate(Chromosome* population, size_t popSize, size_t nCities);
void computeFitness(Chromosome* population, size_t popSize, Point* cities, size_t nCities);
float computePathDistance(Point* cities, size_t nCities, const tag_t *path);
float distance(Point a, Point b);
void mergeSort(Chromosome* in, size_t size);
void merge(Chromosome* a, size_t aSize, Chromosome* b, size_t bSize, Chromosome* c);
void copyPopulation(Chromosome* in, Chromosome* out, size_t popSize, size_t nCities);
void mate(Chromosome *in, size_t inSize, Chromosome *out, size_t outSize, size_t nCities);
int valid(const tag_t *in, size_t nCities);
void printPath(tag_t *path, Point* cities, size_t nCities);


int main(int argc, char** argv)
{
    size_t    epochs = 500;
    size_t    nCities = 250;
    size_t    popSize = 40000;
    double elitism = 0.1;
    seed = 12345;

    // obtain parameters at run time
    switch (argc)
    {
        case 6: elitism = atof(argv[5]);
        case 5: popSize = atoi(argv[4]);
        case 4: nCities = atoi(argv[3]);
        case 3: seed = atoi(argv[2]);
        case 2: epochs = atoi(argv[1]);
        default: break;
    }

    size_t eliteSize = popSize * elitism;
    // Array containing the positions of the cities
    Point* cities = malloc(nCities * sizeof(*cities));

    // Helping structure to account for cities that has been already visited
    int*   mask = malloc(nCities * sizeof(*mask)); // TODO: Change type to char, as a mask this array could be implemented as a simple bitmask

    Chromosome* population = malloc(popSize * sizeof(*population));
    Chromosome* tmpPopulation = malloc(popSize * sizeof(*tmpPopulation));

    printf("Find shortest path for %ld cities. %ld Epochs. population Size: %ld\n", nCities, epochs, popSize);
    printf("Size of Chromosome: %lu bytes\n", sizeof(Chromosome));
    printf("Size of city vector: %lu bytes\n", sizeof(Point) * nCities);
    printf("Size of a permutation: %ld bytes\n", sizeof(tag_t) * nCities);
    printf("Population memory footprint: %ld bytes\n", popSize * (sizeof(Chromosome) + sizeof(tag_t) * nCities));

    generateCities(cities, nCities); // generate random cities and initialize genetic population
    initPopulation(population, popSize, nCities);
    initPopulation(tmpPopulation, popSize, nCities);
    // generate random mutations into initial population
    for (size_t i = 0; i < 10; i++)
    {
        mutate(population, popSize, nCities);
    }
    // compute fitness and sort population by lower fitness, to generate elite
    computeFitness(population, popSize, cities, nCities);
    mergeSort(population, popSize);

    // generate new populations from initial population
    for (size_t i = 0; i < epochs; i++)
    {
        copyPopulation(population, tmpPopulation, eliteSize, nCities);                             // copy elite population to new generation
        mate(population, eliteSize, tmpPopulation + eliteSize, popSize - eliteSize, nCities); // mate from elite
        mutate(tmpPopulation + eliteSize, popSize - eliteSize, nCities);                            // do not affect elite
        copyPopulation(tmpPopulation, population, popSize, nCities);
        computeFitness(population, popSize, cities, nCities);
        mergeSort(population, popSize);                // sort population by lower fitness, to generate new elite

        // display progress
        if (i % 50 == 1)
        {
            // print current best individual
            printf("Fitness: %f\n", population[0].distance);

            // sanity check
            if (!valid(population[0].tour, nCities))
            {
                printf("ERROR: tour is not a valid permutation of cities");
                exit(1);
            }
        }
    }

    // print final result
    printPath(population[0].tour, cities, nCities);
    // sanity check
    if (!valid(population[0].tour, nCities))
    {
        printf("ERROR: tour is not a valid permutation of cities");
        exit(1);
    }

    // free all allocated memory
    for (size_t i = 0; i < popSize; ++i)
    {
        free(population[i].tour);
        free(tmpPopulation[i].tour);
    }

    free(mask);
    free(cities);
    free(population);
    free(tmpPopulation);

    return(0);
}

// Generate random positions for the cities
void generateCities(Point* cities, size_t nCities)
{
    for (size_t i = 0; i < nCities; i++)
    {
        cities[i].x = myRandom() % 4096;
        cities[i].y = myRandom() % 4096;
    }
}

// Initialize a population of popSize Chromosomes
void initPopulation(Chromosome* population, size_t popSize, size_t nCities)
{
    for (size_t i = 0; i < popSize; i++)
    {
        population[i].distance = 0;
        population[i].tour = malloc(nCities * sizeof(*population[i].tour));
        for (size_t j = 0; j < nCities; j++)
        {
            population[i].tour[j] = j;
        }
    }
}

// mutate population: swap cities from two random positions in Chromosome
void mutate(Chromosome* population, size_t popSize, size_t nCities)
{
    for (size_t m = 0; m < popSize; m++)
    {     // generate 2 random positions to swap
        size_t aPos = myRandom() % nCities;
        size_t bPos = myRandom() % nCities;
        int cityA = population[m].tour[aPos];
        population[m].tour[aPos] = population[m].tour[bPos];
        population[m].tour[bPos] = cityA;
    }
}

//  Calculate each individual fitness in population. Fitness is path distance
void computeFitness(Chromosome* population, size_t popSize, Point* cities, size_t nCities)
{
    for (size_t i = 0; i < popSize; i++)
    {
        population[i].distance = computePathDistance(cities, nCities, population[i].tour);
    }
}

// A path is a permutation of cities
// Calculate the total distance of a path
float computePathDistance(Point* cities, size_t nCities, const tag_t *path)
{
    float dist = 0.0f;
    Point citiesInPath[nCities];
    for (size_t i = 0; i < nCities; ++i)
    {
        citiesInPath[i] = cities[path[i]];
    }

    for (size_t i = 1; i < nCities; i++)
    {
        dist = dist + distance(citiesInPath[i - 1], citiesInPath[i]);
    }
    return(dist + distance(citiesInPath[nCities - 1], citiesInPath[0]));
}

// Distance between two Points
float distance(const Point a, const Point b)
{
    return(sqrtf((float)(a.x - b.x) * (float)(a.x - b.x) +
                 (float)(a.y - b.y) * (float)(a.y - b.y)));
}

// Sort array A with n Chromosomes using recursive merge-sort algorithm
void mergeSort(Chromosome* in, size_t size)
{
    if (size < 2)
    {
        return;         // the array is sorted when n=1
    }
    if (size == 2)
    {
        if (in[0].distance > in[1].distance)
        {         // swap values of A[0] and A[1]
            Chromosome temp = in[1];
            in[1] = in[0];
            in[0] = temp;
        }
        return;         // elements sorted
    }

    // divide A into two arrays, A1 and A2
    size_t n1 = size / 2;      // number of elements in A1
    size_t n2 = size - n1;     // number of elements in A2
    Chromosome* firstH = malloc(sizeof(*firstH) * n1);
    Chromosome* secondH = malloc(sizeof(*secondH) * n2);

    // move first n/2 elements to A1
    for (size_t i = 0; i < n1; i++)
    {
        firstH[i] = in[i];         // copy full entry
    }
    // move the rest to A2
    for (size_t i = 0; i < n2; i++)
    {
        secondH[i] = in[i + n1];         // copy full entry
    }

    // recursive calls
    mergeSort(firstH, n1);
    mergeSort(secondH, n2);

    // merge
    merge(firstH, n1, secondH, n2, in);

    // free allocated memory
    free(firstH);
    free(secondH);
}

// Merge two sorted arrays, A with szA Chromosomes and B with szB Chromosomes, into a sorted array C
void merge(Chromosome* a, size_t aSize, Chromosome* b, size_t bSize, Chromosome* c)
{
    size_t i = 0, j = 0;

    while (i + j < aSize + bSize)
    {
        if (j == bSize || ((i < aSize) && (a[i].distance <= b[j].distance)))
        {
            c[i + j] = a[i];             // copy full struct
            i++;
        }
        else
        {
            c[i + j] = b[j];             // copy full struct
            j++;
        }
    }
}

// copy input population to output population
void copyPopulation(Chromosome* in, Chromosome* out, size_t popSize, size_t nCities)
{
    for (size_t i = 0; i < popSize; i++)
    {
        for (size_t j = 0; j < nCities; j++)
        {
            out[i].tour[j] = in[i].tour[j];
        }
    }
}

// mate randomly the elite population in in into P_out
void mate(Chromosome *in, size_t inSize, Chromosome *out, size_t outSize, size_t nCities)
{
    // Declare local mask
    tag_t mask[nCities];

    // mate the elite population to generate new genes
    for (size_t m = 0; m < outSize; m++)
    {
        // Create new gene in Output population by mating to genes from the elite input population
        // select two random genes from elite population and mate them at random position pos
        size_t i1 = myRandom() % inSize;
        size_t i2 = myRandom() % inSize;
        size_t pos = myRandom() % nCities;
        const tag_t* parentA = in[i1].tour;
        const tag_t* parentB = in[i2].tour;
        tag_t* child = out[m].tour;

        // Clear mask of already visited cities
        memset(mask, 0, nCities * sizeof(*mask));

        // Copy first part of parent A to child
        for (size_t i = 0; i < pos; i++)
        {
            tag_t city = parentA[i];
            child[i] = city;
            mask[city] = 1;
        }

        // Copy cities in parent B to last part of child, maintaining the ordering in parent B
        // Copy those cities that are not in the first part of parent A
        size_t j = 0;         // Points to the consecutive positions in parent B
        tag_t city = parentB[j];
        for (size_t i = pos; i < nCities; i++)
        {
            while (mask[city] == 1)               // skip cities in parent B already visited
            {
                j++;
                city = parentB[j];
            }

            mask[city] = 1;                     // mark city as seen
            child[i] = city;            // copy city to child
        }
    }

}

// Checks is a path is valid: does not contain repeats
int valid(const tag_t *in, size_t nCities)
{
    // clear mask
    tag_t mask[nCities];
    // Clear mask of already visited cities
    memset(mask, 0, nCities * sizeof(*mask));

    // check if city has been already visited, otherwise insert city in mask
    for (size_t i = 0; i < nCities; i++)
    {
        if (mask[in[i]] == 0)
        {
            mask[in[i]] = 1;
        }
        else
        {
            return(0);
        }
    }

    return(1);
}

// Display path into screen
void printPath(tag_t *path, Point* cities, size_t nCities)
{
    size_t   i;
    float dist = 0.0f;

    for (i = 1; i < nCities; i++)
    {
        printf("%d,", path[i - 1]);
        dist = dist + distance(cities[path[i - 1]], cities[path[i]]);
    }
    printf("%d\nTotal Distance: %f\n", path[i - 1],
           dist + distance(cities[path[i - 1]], cities[path[0]]));
}
