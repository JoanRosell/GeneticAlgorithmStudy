/**
 * description:
 * 	Genetic algorithm for finding heuristic solution of Traveller Salesman Problem
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//  Define Constants: Size of chromosome population, number of cities and size of elite
#define POPSIZE	    40000
#define NUM_CITIES  250
#define ELITSIZE    ((int)(POPSIZE*0.1f))

// Variable used to generate pseudo-random numbers
unsigned int seed;


// Function to generate pseudo-random numbers
inline int myRandom() {
  seed = (214013*seed+2531011);
  return (seed>>13);
}

// Type definition of a point
typedef struct
{
  int x, y;
} point;


// Distance between two points
float distance ( const point P1, const point P2 )
{
  return sqrtf ( (float)(P1.x-P2.x)*(float)(P1.x-P2.x) +
                 (float)(P1.y-P2.y)*(float)(P1.y-P2.y) );
}

// Array containing the positions of the cities
point Cities[NUM_CITIES];

// Generate random positions for the cities
void Random_Cities()
{
  int i;
  for(i=0; i<NUM_CITIES; i++)
  {
    Cities[i].x = myRandom() % 4096;
    Cities[i].y = myRandom() % 4096;
  }
}

// A path is a permutation of cities
// Calculate the total distance of a path
float compute_path_distance( int *path )
{
  int i;
  float dist=0.0f;
  for (i=1; i<NUM_CITIES; i++)
    dist = dist + distance( Cities[path[i-1]], Cities[path[i]] );
  return dist + distance( Cities[path[NUM_CITIES-1]], Cities[path[0]] );
}

// Display path into screen
void print_path( int *path )
{
  int i;
  float dist=0.0f;
  for(i=1; i<NUM_CITIES; i++)
  {
    printf("%d,", path[i-1] );
    dist = dist + distance( Cities[path[i-1]], Cities[path[i]] );
  }
  printf("%d\nTotal Distance: %f\n", path[i-1],
         dist + distance( Cities[path[i-1]], Cities[path[0]]));
}


// Structure base of chromosome. A chromosome is interpreted as a path
typedef struct
{
  float distance;
  int   gen[NUM_CITIES];
} chrom;


// Initialize a population of N chromosomes
void init_population ( chrom *P, int N )
{
  int i, j;
  for (i=0; i<N; i++)
  {
    P[i].distance = 0;
    for (j=0; j<NUM_CITIES; j++)
      P[i].gen[j] = j;
  }
}


//  Calculate each individual fitness in population. Fitness is path distance
void compute_fitness ( chrom *P, int N )
{
  int i;
  for (i=0; i<N; i++)
    P[i].distance = compute_path_distance( P[i].gen );
}


// Merge two sorted arrays, A with szA chromosomes and B with szB chromosomes, into a sorted array C
void merge (chrom *A, int szA, chrom *B, int szB, chrom *C)
{
  int i= 0, j= 0;
  while (i+j < szA+szB)
  {
    if (j==szB || ((i<szA) && (A[i].distance <= B[j].distance)))
    {
      C[i+j] = A[i];  // copy full struct
      i++;
    }
    else
    {
      C[i+j] = B[j];  // copy full struct
      j++;
    }
  }
}


// Sort array A with n chromosomes using recursive merge-sort algorithm
void merge_sort (chrom *A, int n)
{
  int i, n1, n2;
  chrom *A1, *A2;

  if (n < 2)
    return;   // the array is sorted when n=1
  if (n==2)
  {
    if (A[0].distance > A[1].distance)
    {  // swap values of A[0] and A[1]
      chrom T = A[1];
      A[1] = A[0];
      A[0] = T;
    }
    return; // elements sorted
  }

  // divide A into two arrays, A1 and A2
  n1 = n / 2;   // number of elements in A1
  n2 = n - n1;  // number of elements in A2
  A1 = (chrom *)malloc(sizeof(chrom) * n1);
  A2 = (chrom *)malloc(sizeof(chrom) * n2);

  // move first n/2 elements to A1
  for (i=0; i < n1; i++) {
    A1[i] = A[i]; // copy full entry
  }
  // move the rest to A2
  for (i=0; i < n2; i++) {
    A2[i] = A[i+n1]; // copy full entry
  }

  // recursive calls
  merge_sort(A1, n1);
  merge_sort(A2, n2);

  // merge
  merge(A1, n1, A2, n2, A);

  // free allocated memory
  free(A1);
  free(A2);
}


// copy input population to output population
void copy_population ( chrom *P_in, chrom *P_out, int N )
{
  int i, j;
  for (i=0; i<N; i++)
    for (j=0; j<NUM_CITIES; j++)
      P_out[i].gen[j] = P_in[i].gen[j];
}


// Helping structure to account for cities that has been already visited
int MASK[NUM_CITIES];

// Checks is a path is valid: does not contain repeats
int check_valid( int *IN )
{
  int i;

  // clear mask
  for (i=0; i<NUM_CITIES; i++)
    MASK[i]=0;

  // check if city has been already visited, otherwise insert city in mask
  for (i=0; i<NUM_CITIES; i++)
    if (MASK[IN[i]] == 0)
      MASK[IN[i]] = 1;
    else
      return 0;

  return 1;
}


// mate randomly the elite population in P_in into P_out
void mate ( chrom *P_in, int Nin, chrom *P_out, int Nout )
{
  int m, i, i1, i2, pos, city, j;

  // mate the elite population to generate new genes
  for (m=0; m < Nout; m++)
  {
    // Create new gene in Output population by mating to genes from the elite input population
    // select two random genes from elite population and mate them at random position pos
    i1  = myRandom () % Nin;
    i2  = myRandom () % Nin;
    pos = myRandom () % NUM_CITIES;

    // Clear mask of already visited cities
    for (i=0; i<NUM_CITIES; i++)
      MASK[i]=0;

    // Copy first part of input gene i1 to output gene
    for (i = 0; i < pos; i++)
      P_out[m].gen[i] = P_in[i1].gen[i];

    // Mark all cities in first part of output gene i1
    for (i = 0; i < pos; i++)
    {
      city = P_out[m].gen[i];
      MASK[city] = 1;
    }

    // copy cities in input gene i2 to last part of output gene,
    //    maintaining the ordering in gene i2
    // copy those cities that are not in the first part of gene i1
    j=0; // points to the consecutive positions in gen i2
    for (i = pos; i<NUM_CITIES; i++)
    {
      do { // skip cities in gen i2 already visited
        city = P_in[i2].gen[j];
        j++;
      } while ( MASK[city] ==1 );

      MASK[city] = 1;         // mark city as seen
      P_out[m].gen[i] = city; // copy city to output gene
    }
  }
}


// mutate population: swap cities from two random positions in chromosome
void mutate ( chrom *Pop, int N )
{
  int m;

  for (m=0; m < N; m++)
  { // generate 2 random positions to swap
    int apos = myRandom () % NUM_CITIES;
    int bpos = myRandom () % NUM_CITIES;
    int CityA   = Pop[m].gen[apos];
    Pop[m].gen[apos] = Pop[m].gen[bpos];
    Pop[m].gen[bpos] = CityA;
  }
}

chrom Population[POPSIZE], TempPopulation[POPSIZE];

int main (int argc, char **argv)
{
  int i, R, EPOCHS=500;
  seed = 12345;

  // obtain parameters at run time
  if (argc>1) {  EPOCHS = atoi(argv[1]); }
  if (argc>2) {  seed   = atoi(argv[2]); }

  printf("Find shortest path for %d cities. %d Epochs. Population Size: %d\n",
         NUM_CITIES, EPOCHS, POPSIZE);

  // generate random cities and initialize genetic population
  Random_Cities();
  init_population ( Population, POPSIZE );
  // generate random mutations into initial population
  for (i=0; i<10; i++)
    mutate ( Population, POPSIZE );
  // compute fitness and sort population by lower fitness, to generate elite
  compute_fitness ( Population, POPSIZE );
  merge_sort      ( Population, POPSIZE );

  // generate new populations from initial population
  for (i=0; i<EPOCHS; i++)
  {
    copy_population ( Population, TempPopulation, ELITSIZE );    // copy elite population to new generation
    mate            ( Population, ELITSIZE, TempPopulation+ELITSIZE, POPSIZE-ELITSIZE ); // mate from elite
    mutate          ( TempPopulation+ELITSIZE,    POPSIZE-ELITSIZE );  // do not affect elite
    copy_population ( TempPopulation, Population, POPSIZE );
    compute_fitness ( Population, POPSIZE );
    merge_sort      ( Population, POPSIZE );  // sort population by lower fitness, to generate new elite

    // display progress
    if (i % 50 == 1)
    {
      // print current best individual
      printf ("Fitness: %f\n", Population[0].distance);

      // sanity check
      if ( !check_valid( Population[0].gen ) )
      {
        printf("ERROR: gen is not a valid permutation of Cities");
        exit(1);
      }
    }
  }

  // print final result
  print_path( Population[0].gen );
  // sanity check
  if ( !check_valid( Population[0].gen ) )
  {
    printf("ERROR: gen is not a valid permutation of Cities");
    exit(1);
  }

  return 0;
}
