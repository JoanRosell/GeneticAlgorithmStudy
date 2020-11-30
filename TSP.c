/**
 * description:
 *  Genetic algorithm for finding heuristic solution of Travelling Salesman Problem
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Variable used to generate pseudo-random numbers
unsigned int seed;

// Type definition of a Point
typedef struct Point
{
    int x;
    int y;
} Point;

// Structure base of chromosome. A chromosome is interpreted as a path
typedef struct Chromosome
{
    float distance;
    int*  tour;
} Chromosome;

// Function to generate pseudo-random numbers
int myRandom()
{
    seed = (214013 * seed + 2531011);
    return(seed >> 13);
}

// Function definitions
void generateCities(Point* cities, int nCities);
void initPopulation(Chromosome* population, int popSize, int nCities);
void mutate(Chromosome* population, int popSize, int nCities);
void computeFitness(Chromosome* population, int popSize, Point* cities, int nCities);
float computePathDistance(Point* cities, int nCities, const int* path);
float distance(Point a, Point b);
void mergeSort(Chromosome* in, int size);
void merge(Chromosome* a, int aSize, Chromosome* b, int bSize, Chromosome* c);
void copyPopulation(Chromosome* in, Chromosome* out, int popSize, int nCities);
void mate(Chromosome* in, int inSize, Chromosome* out, int outSize, int* mask, int nCities);
int valid(const int* in, int* mask, int nCities);
void printPath(int* path, Point* cities, int nCities);


int main(int argc, char** argv)
{
    int    epochs = 500;
    int    nCities = 250; // TODO: Change type to size_t
    int    popSize = 40000; // TODO: Change type to size_t
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

    int eliteSize = popSize * elitism; // TODO: Change type to size_t
    // Array containing the positions of the cities
    Point* cities = malloc(nCities * sizeof(*cities));

    // Helping structure to account for cities that has been already visited
    int*   mask = malloc(nCities * sizeof(*mask)); // TODO: Change type to char, as a mask this array could be implemented as a simple bitmask

    Chromosome* population = malloc(popSize * sizeof(*population));
    Chromosome* tmpPopulation = malloc(popSize * sizeof(*tmpPopulation));

    printf("Find shortest path for %d cities. %d Epochs. population Size: %d\n", nCities, epochs, popSize);
    printf("Size of Chromosome: %lu bytes\n", sizeof(Chromosome));
    printf("Size of Point: %lu bytes\n", sizeof(Point));

    generateCities(cities, nCities); // generate random cities and initialize genetic population
    initPopulation(population, popSize, nCities);
    initPopulation(tmpPopulation, popSize, nCities);
    // generate random mutations into initial population
    for (int i = 0; i < 10; i++)
    {
        mutate(population, popSize, nCities);
    }
    // compute fitness and sort population by lower fitness, to generate elite
    computeFitness(population, popSize, cities, nCities);
    mergeSort(population, popSize);

    // generate new populations from initial population
    for (int i = 0; i < epochs; i++)
    {
        copyPopulation(population, tmpPopulation, eliteSize, nCities);                             // copy elite population to new generation
        mate(population, eliteSize, tmpPopulation + eliteSize, popSize - eliteSize, mask, nCities); // mate from elite
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
            if (!valid(population[0].tour, mask, nCities))
            {
                printf("ERROR: tour is not a valid permutation of cities");
                exit(1);
            }
        }
    }

    // print final result
    printPath(population[0].tour, cities, nCities);
    // sanity check
    if (!valid(population[0].tour, mask, nCities))
    {
        printf("ERROR: tour is not a valid permutation of cities");
        exit(1);
    }

    // free all allocated memory
    for (int i = 0; i < popSize; ++i)
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
void generateCities(Point* cities, int nCities)
{
    for (int i = 0; i < nCities; i++)
    {
        cities[i].x = myRandom() % 4096;
        cities[i].y = myRandom() % 4096;
    }
}

// Initialize a population of popSize Chromosomes
void initPopulation(Chromosome* population, int popSize, int nCities)
{
    for (int i = 0; i < popSize; i++)
    {
        population[i].distance = 0;
        population[i].tour = malloc(nCities * sizeof(*population[i].tour));
        for (int j = 0; j < nCities; j++)
        {
            population[i].tour[j] = j;
        }
    }
}

// mutate population: swap cities from two random positions in Chromosome
void mutate(Chromosome* population, int popSize, int nCities)
{
    for (int m = 0; m < popSize; m++)
    {     // generate 2 random positions to swap
        int aPos = myRandom() % nCities;
        int bPos = myRandom() % nCities;
        int cityA = population[m].tour[aPos];
        population[m].tour[aPos] = population[m].tour[bPos];
        population[m].tour[bPos] = cityA;
    }
}

//  Calculate each individual fitness in population. Fitness is path distance
void computeFitness(Chromosome* population, int popSize, Point* cities, int nCities)
{
    for (int i = 0; i < popSize; i++)
    {
        population[i].distance = computePathDistance(cities, nCities, population[i].tour);
    }
}

// A path is a permutation of cities
// Calculate the total distance of a path
float computePathDistance(Point* cities, int nCities, const int* path)
{
    float dist = 0.0f;

    for (int i = 1; i < nCities; i++)
    {
        dist = dist + distance(cities[path[i - 1]], cities[path[i]]);
    }
    return(dist + distance(cities[path[nCities - 1]], cities[path[0]]));
}

// Distance between two Points
float distance(const Point a, const Point b)
{
    return(sqrtf((float)(a.x - b.x) * (float)(a.x - b.x) +
                 (float)(a.y - b.y) * (float)(a.y - b.y)));
}

// Sort array A with n Chromosomes using recursive merge-sort algorithm
void mergeSort(Chromosome* in, int size)
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
    int n1 = size / 2;      // number of elements in A1
    int n2 = size - n1;     // number of elements in A2
    Chromosome* firstH = malloc(sizeof(*firstH) * n1);
    Chromosome* secondH = malloc(sizeof(*secondH) * n2);

    // move first n/2 elements to A1
    for (int i = 0; i < n1; i++)
    {
        firstH[i] = in[i];         // copy full entry
    }
    // move the rest to A2
    for (int i = 0; i < n2; i++)
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
void merge(Chromosome* a, int aSize, Chromosome* b, int bSize, Chromosome* c)
{
    int i = 0, j = 0;

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
void copyPopulation(Chromosome* in, Chromosome* out, int popSize, int nCities)
{
    for (int i = 0; i < popSize; i++)
    {
        for (int j = 0; j < nCities; j++)
        {
            out[i].tour[j] = in[i].tour[j];
        }
    }
}

// mate randomly the elite population in in into P_out
void mate(Chromosome* in, int inSize, Chromosome* out, int outSize, int* mask, int nCities)
{
    // mate the elite population to generate new genes
    for (int m = 0; m < outSize; m++)
    {
        // Create new gene in Output population by mating to genes from the elite input population
        // select two random genes from elite population and mate them at random position pos
        int i1 = myRandom() % inSize;
        int i2 = myRandom() % inSize;
        int pos = myRandom() % nCities;

        // Clear mask of already visited cities
        for (int i = 0; i < nCities; i++)
        {
            mask[i] = 0;
        }

        // Copy first part of input gene i1 to output gene
        for (int i = 0; i < pos; i++)
        {
            out[m].tour[i] = in[i1].tour[i];
        }

        // Mark all cities in first part of output gene i1
        for (int i = 0; i < pos; i++)
        {
            int city = out[m].tour[i];
            mask[city] = 1;
        }

        // copy cities in input gene i2 to last part of output gene,
        //    maintaining the ordering in gene i2
        // copy those cities that are not in the first part of gene i1
        int j = 0;         // Points to the consecutive positions in tour i2
        int city = in[i2].tour[j];
        for (int i = pos; i < nCities; i++)
        {
            while (mask[city] == 1)               // skip cities in tour i2 already visited
            {
                j++;
                city = in[i2].tour[j];
            }

            mask[city] = 1;                     // mark city as seen
            out[m].tour[i] = city;            // copy city to output gene
        }
    }
}

// Checks is a path is valid: does not contain repeats
int valid(const int* in, int* mask, int nCities)
{
    // clear mask
    for (int i = 0; i < nCities; i++)
    {
        mask[i] = 0;
    }

    // check if city has been already visited, otherwise insert city in mask
    for (int i = 0; i < nCities; i++)
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
void printPath(int* path, Point* cities, int nCities)
{
    int   i;
    float dist = 0.0f;

    for (i = 1; i < nCities; i++)
    {
        printf("%d,", path[i - 1]);
        dist = dist + distance(cities[path[i - 1]], cities[path[i]]);
    }
    printf("%d\nTotal Distance: %f\n", path[i - 1],
           dist + distance(cities[path[i - 1]], cities[path[0]]));
}
