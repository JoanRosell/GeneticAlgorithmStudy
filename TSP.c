/**
 * description:
 *  Genetic algorithm for finding heuristic solution of Travelling Salesman Problem
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

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

int comp(const void* lhs, const void* rhs)
{
    Chromosome a = *((Chromosome* ) lhs);
    Chromosome b = *((Chromosome* ) rhs);
    return (a.distance > b.distance) - (a.distance < b.distance);
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
int valid(const tag_t *in, size_t nCities);
void printPath(tag_t *path, Point* cities, size_t nCities);

void swapChromosomes(Chromosome* a, Chromosome* b)
{
    float tmpDistance = a->distance;
    a->distance = b->distance;
    b->distance = tmpDistance;

    tag_t* tmpTour = a->tour;
    a->tour = b->tour;
    b->tour = tmpTour;
}

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
    Point* cities = malloc(nCities * sizeof(*cities));

    Chromosome* population = malloc(popSize * sizeof(*population));

    printf("Find shortest path for %ld cities. %ld Epochs. population Size: %ld\n", nCities, epochs, popSize);

    generateCities(cities, nCities); // generate random cities and initialize genetic population
    initPopulation(population, popSize, nCities);

    // generate random mutations into initial population
    for (size_t i = 0; i < 10; i++)
    {
        mutate(population, popSize, nCities);
    }

    // compute fitness and sort population by lower fitness, to generate elite
    computeFitness(population, popSize, cities, nCities);
    mergeSort(population, popSize);
    size_t outSize = popSize - eliteSize;

    uint32_t abMateVector[outSize];
    tag_t cMateVector[outSize];
    uint16_t abMutateVector[outSize];
    Chromosome* candidates = population + eliteSize;
    size_t lastSwap;
    #pragma omp parallel num_threads(8)
    for (size_t e = 0; e < epochs; e++)
    {
        #pragma omp single
        {
            for (size_t m = 0; m < outSize; m++)
            {
                uint16_t a = myRandom() % eliteSize;
                uint16_t b = myRandom() % eliteSize;
                abMateVector[m] = ((uint32_t) a << 16) + b;
                cMateVector[m] = myRandom() % nCities;
            }

            for (size_t m = 0; m < outSize; m++)
            {
                tag_t a = myRandom() % nCities;
                tag_t b = myRandom() % nCities;
                abMutateVector[m] = ((uint16_t) a << 8) + b;
            }
        }

        // mate the elite population to generate new genes
        #pragma omp for schedule(dynamic)
        for (size_t m = 0; m < outSize; m++)
        {
            tag_t mask[nCities];
            memset(mask, 0xFF, nCities * sizeof(*mask));
            uint16_t i1 = (uint16_t) (abMateVector[m] >> 16);
            uint16_t i2 = (uint16_t) abMateVector[m];
            tag_t pos = cMateVector[m];
            const tag_t* parentA = population[i1].tour;
            const tag_t* parentB = population[i2].tour;
            tag_t child[nCities];

            // Copy first part of parent A to child
            memcpy(child, parentA, pos * sizeof(*child));

            for (size_t i = 0; i < pos; i++)
            {
                mask[parentA[i]] = 0;
            }

            size_t k = pos;
            for (size_t i = 0; i < nCities; ++i)
            {
                tag_t tmp = mask[parentB[i]];
                child[k] = (parentB[i] & tmp) | (child[k] & ~tmp);
                k += tmp & 1;
            }

            // Mutate recently generated child
            tag_t aPos = (tag_t) (abMutateVector[m] >> 8);
            tag_t bPos = (tag_t) abMutateVector[m];
            tag_t cityA = child[aPos];
            child[aPos] = child[bPos];
            child[bPos] = cityA;

            memcpy(population[m + eliteSize].tour, child, nCities * sizeof(*child));
            memset(mask, 0xFF, nCities * sizeof(*mask));
        }

        computeFitness(population + eliteSize, popSize - eliteSize, cities, nCities);

        #pragma omp single
        {
            float worstEliteScore = population[eliteSize - 1].distance;
            size_t j = 0;

            for (size_t i = 0; i < outSize; i++)
            {
                while (candidates[j].distance >= worstEliteScore) j++;
                if (j == outSize) break;
                swapChromosomes(&candidates[i], &candidates[j]);
                j++;
                lastSwap = i + 1;
            }

            mergeSort(population, eliteSize + lastSwap);

            if (e % 50 == 1) {
                // print current best individual
                printf("Fitness: %f\n", population[0].distance);

                // sanity check
                if (!valid(population[0].tour, nCities)) {
                    printf("ERROR: tour is not a valid permutation of cities");
                    exit(1);
                }
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
    }


    free(cities);
    free(population);

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
    #pragma omp for schedule(static)
    for (size_t i = 0; i < popSize; i++)
    {
        population[i].distance = computePathDistance(cities, nCities, population[i].tour);
    }
}

// A path is a permutation of cities
// Calculate the total distance of a path
// IDEA:
//  Store distances as a lower triangular matrix instead of computing them in every iteration
//  For matrix n by n you need array (n+1)*n/2 length and transition rule is Matrix[i][j] = Array[i*(i+1)/2+j]
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
    if (size <= 512)
    {
        qsort(in, size, sizeof(*in), comp);
        return;
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
    #pragma omp task
    mergeSort(firstH, n1);
    #pragma omp task
    mergeSort(secondH, n2);
    #pragma omp taskwait
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
