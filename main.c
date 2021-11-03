////////////////////////////////////////////////////////////////////////////////////////////////////
//
// AED, 2019/2020
//
// 93221 - Pedro Santos
// 93283 - Pedro Amaral
// 95278 - Diogo Cunha
//
// Brute-force solution of the assignment problem (https://en.wikipedia.org/wiki/Assignment_problem)
//
// Compile with "cc -Wall -O2 assignment.c -lm" or equivalent
//
// In the assignment problem we will solve here we have n agents and n tasks; assigning agent
//   a
// to task
//   t
// costs
//   cost[a][t]
// The goal of the problem is to assign one agent to each task such that the total cost is minimized
// The total cost is the sum of the costs
//
// Things to do:
//   0. (mandatory)
//      Place the student numbers and names at the top of this file
//   1. (highly recommended)
//      Read and understand this code
//   2. (mandatory)
//      Modify the function generate_all_permutations to solve the assignment problem
//      Compute the best and worst solutions for all problems with sizes n=2,...,14 and for each
//      student number of the group
//   3. (mandatory)
//      Calculate and display an histogram of the number of occurrences of each cost
//      Does it follow approximately a normal distribution?
//      Note that the maximum possible cost is n * t_range
//   4. (optional)
//      For each problem size, and each student number of the group, generate one million (or more!)
//      random permutations and compute the best and worst solutions found in this way; compare
//      these solutions with the ones found in item 2
//      Compare the histogram computed in item 3 with the histogram computed using the random
//      permutations
//   5. (optional)
//      Try to improve the execution time of the program (use the branch-and-bound technique)
//   6. (optional)
//      Surprise us, by doing something more!
//   7. (mandatory)
//      Write a report explaining what you did and presenting your results
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define NDEBUG  // uncomment to skip disable asserts (makes the code slightly faster)
#include <assert.h>


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// problem data
//
// max_n ........ maximum problem size
// cost[a][t] ... cost of assigning agent a to task t
//

//
// if your compiler complains about srandom() and random(), replace #if 0 by #if 1
//
#if 0
# define srandom srand
# define random  rand
#endif

#define max_n    32           // do not change this (maximum number of agents, and tasks)
#define range    20           // do not change this (for the pseudo-random generation of costs)
#define t_range  (3 * range)  // do not change this (maximum cost of an assignment)

static int cost[max_n][max_n];
static int seed; // place a student number here!

static void init_costs(int n)
{
  if(n == -3)
  { // special case (example for n=3)
    cost[0][0] = 3; cost[0][1] = 8; cost[0][2] = 6;
    cost[1][0] = 4; cost[1][1] = 7; cost[1][2] = 5;
    cost[2][0] = 5; cost[2][1] = 7; cost[2][2] = 5;
    return;
  }
  if(n == -5)
  { // special case (example for n=5)
    cost[0][0] = 27; cost[0][1] = 27; cost[0][2] = 25; cost[0][3] = 41; cost[0][4] = 24;
    cost[1][0] = 28; cost[1][1] = 26; cost[1][2] = 47; cost[1][3] = 38; cost[1][4] = 21;
    cost[2][0] = 22; cost[2][1] = 48; cost[2][2] = 26; cost[2][3] = 14; cost[2][4] = 24;
    cost[3][0] = 32; cost[3][1] = 31; cost[3][2] =  9; cost[3][3] = 41; cost[3][4] = 36;
    cost[4][0] = 24; cost[4][1] = 34; cost[4][2] = 30; cost[4][3] = 35; cost[4][4] = 45;
    return;
  }
  assert(n >= 1 && n <= max_n);
  srandom((unsigned int)seed * (unsigned int)max_n + (unsigned int)n);
  for(int a = 0;a < n;a++)
    for(int t = 0;t < n;t++)
      cost[a][t] = 3 + (random() % range) + (random() % range) + (random() % range); // [3,3*range]
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// code to measure the elapsed time used by a program fragment (an almost copy of elapsed_time.h)
//
// use as follows:
//
//   (void)elapsed_time();
//   // put your code to be time measured here
//   dt = elapsed_time();
//   // put morecode to be time measured here
//   dt = elapsed_time();
//
// elapsed_time() measures the CPU time between consecutive calls
//

#if defined(__linux__) || defined(__APPLE__)

//
// GNU/Linux and MacOS code to measure elapsed time
//

#include <time.h>

static double elapsed_time(void)
{
  static struct timespec last_time,current_time;

  last_time = current_time;
  if(clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&current_time) != 0)
    return -1.0; // clock_gettime() failed!!!
  return            ((double)current_time.tv_sec - (double)last_time.tv_sec)
         + 1.0e-9 * ((double)current_time.tv_nsec - (double)last_time.tv_nsec);
}

#endif

#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)

//
// Microsoft Windows code to measure elapsed time
//

#include <windows.h>

static double elapsed_time(void)
{
  static LARGE_INTEGER frequency,last_time,current_time;
  static int first_time = 1;

  if(first_time != 0)
  {
    QueryPerformanceFrequency(&frequency);
    first_time = 0;
  }
  last_time = current_time;
  QueryPerformanceCounter(&current_time);
  return (double)(current_time.QuadPart - last_time.QuadPart) / (double)frequency.QuadPart;
}

#endif


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// function to generate a pseudo-random permutation
//

void random_permutation(int n,int t[n])
{
  assert(n >= 1 && n <= 1000000);
  for(int i = 0;i < n;i++)
    t[i] = i;
  for(int i = n - 1;i > 0;i--)
  {
    int j = (int)floor((double)(i + 1) * (double)random() / (1.0 + (double)RAND_MAX)); // range 0..i
    assert(j >= 0 && j <= i);
    int k = t[i];
    t[i] = t[j];
    t[j] = k;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// place to store best and worst solutions (also code to print them)
//

static int min_cost,min_cost_assignment[max_n]; // smallest cost information
static int max_cost,max_cost_assignment[max_n]; // largest cost information
static long n_visited; // number of permutations visited (examined)
static int histogram[14*t_range]; // histogram global variable
static double cpu_time;

#define minus_inf  -1000000000  // a very small integer
#define plus_inf   +1000000000  // a very large integer

static void reset_solutions(void)
{
  min_cost = plus_inf;
  max_cost = minus_inf;
  n_visited = 0l;
  // place your histogram initialization code here
  memset(histogram, 0, sizeof(histogram));
  cpu_time = 0.0;
}

#define show_info_1        (1 << 0)
#define show_info_2        (1 << 1)
#define show_costs         (1 << 2)
#define show_min_solution  (1 << 3)
#define show_max_solution  (1 << 4)
#define show_histogram     (1 << 5)
#define show_all           (0xFFFF)

static void show_solutions(int n,char *header,int what_to_show)
{
  printf("%s\n",header);
  if((what_to_show & show_info_1) != 0)
  {
    printf("  seed .......... %d\n",seed);
    printf("  n ............. %d\n",n);
  }
  if((what_to_show & show_info_2) != 0)
  {
    printf("  visited ....... %ld\n",n_visited);
    printf("  cpu time ...... %.3fs\n",cpu_time);
  }
  if((what_to_show & show_costs) != 0)
  {
    printf("  costs .........");
    for(int a = 0;a < n;a++)
    {
      for(int t = 0;t < n;t++)
        printf(" %2d",cost[a][t]);
      printf("\n%s",(a < n - 1) ? "                 " : "");
    }
  }
  if((what_to_show & show_min_solution) != 0)
  {
    printf("  min cost ...... %d\n",min_cost);
    if(min_cost != plus_inf)
    {
      printf("  assignment ...");
      for(int i = 0;i < n;i++)
        printf(" %d",min_cost_assignment[i]);
      printf("\n");
    }
  }
  if((what_to_show & show_max_solution) != 0)
  {
    printf("  max cost ...... %d\n",max_cost);
    if(max_cost != minus_inf)
    {
      printf("  assignment ...");
      for(int i = 0;i < n;i++)
        printf(" %d",max_cost_assignment[i]);
      printf("\n");
    }
  }
  if((what_to_show & show_histogram) != 0)
  {
    // place your code to print the histogram here
   
	printf("\nHistogram n = %d\n", n);
	printf("i\tn[i]\n");
    for (int i = 0; i < n*t_range; i++)
    {
		printf("%d\t%d\n", i, histogram[i]);
	}
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// code used to generate all permutations of n objects (brute force)
//
// n ........ number of objects
// m ........ index where changes occur (a[0], ..., a[m-1] will not be changed)
// a[idx] ... the number of the object placed in position idx
//

static void generate_all_permutations(int n,int m,int a[n])
{
  if(m < n - 1)
  {
    //
    // not yet at the end; try all possibilities for a[m]
    //
    for(int i = m;i < n;i++)
    {
#define swap(i,j)  do { int t = a[i]; a[i] = a[j]; a[j] = t; } while(0)
      swap(i,m);                            // exchange a[i] with a[m]
      generate_all_permutations(n,m + 1,a); // recurse
      swap(i,m);                            // undo the exchange of a[i] with a[m]
#undef swap
    }
  }
  else
  {
    //
    // visit the permutation
    //
    int totalcost=0;
    for (int j=0; j<n; j++) {
        totalcost += cost[j][a[j]];
    }

    n_visited++;

    // place your code to update the best and worst solutions, and to update the histogram here
    if (max_cost < totalcost) {max_cost = totalcost; memcpy(max_cost_assignment, a, sizeof(max_cost_assignment));}
    if (min_cost > totalcost) {min_cost = totalcost; memcpy(min_cost_assignment, a, sizeof(min_cost_assignment));}
	histogram[totalcost]++;

  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// code used to solve the problem with branch and bound
//
// n ........ number of objects
// m ........ index where changes occur (a[0], ..., a[m-1] will not be changed)
// a[idx] ... the number of the object placed in position idx
// custoparcial
//

static void generate_all_permutations_branch_and_bound(int n,int m,int a[n], int custoparcial)
{
    
    if(m < n - 1)
    {
        for(int i = m;i < n;i++)
        {   
            if (custoparcial > min_cost) break;
#define swap(i,j)  do { int t = a[i]; a[i] = a[j]; a[j] = t; } while(0)
            swap(i,m);                            // exchange a[i] with a[m]
            generate_all_permutations_branch_and_bound(n,m + 1,a, custoparcial + cost[m][a[m]]); // recurse
            swap(i,m);                            // undo the exchange of a[i] with a[m]
#undef swap
        }
    }
    else
    {
    //
    // visit the permutation and update the best solution
    //
    if (custoparcial + cost[m][a[m]] < min_cost) {
        min_cost = custoparcial + cost[m][a[m]];
        memcpy(min_cost_assignment, a, sizeof(min_cost_assignment));
    }
    
    n_visited++;
    histogram[custoparcial + cost[m][a[m]]]++;
    
  }
    
}


///////////////////////////////////////////////////////////////////////////////////////////////////////

//code to solve  the problem with random permutations
//
// n ........ number of objects
// t[n] ..... permutation

void calculate_solutions(int n, int t[n])
{
	for (int i = 0; i < 1000000; i++)
	{
		random_permutation(n, t);
		int totalcost=0;
		for (int j=0; j<n; j++) {
			totalcost += cost[j][t[j]];
		}
		
		n_visited++;
		
		if (max_cost < totalcost) {max_cost = totalcost; memcpy(max_cost_assignment, t, sizeof(max_cost_assignment));}
		if (min_cost > totalcost) {min_cost = totalcost; memcpy(min_cost_assignment, t, sizeof(min_cost_assignment));}
		histogram[totalcost]++;
	}
		
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//code to print the data to a file and use in matlab
//
// n - number of objects
// o - 1->brute force, 2->branch and bound, 3->random permutations

void print_to_file(int n, int o)
{
    FILE *fp;
    switch (o)
    {
		case 1:  //brute force
            // histogram
            if (n==13) {
                fp = fopen("histogram13.txt", "w");
                for (int i = 0; i < n*t_range; i++)
                {
                    fprintf(fp, "%d\t%d\n", i, histogram[i]);
                }
                fclose(fp);
            }
            if (n==14) {
                fp = fopen("histogram14.txt", "w");
                for (int i = 0; i < n*t_range; i++)
                {
                    fprintf(fp, "%d\t%d\n", i, histogram[i]);
                }
                fclose(fp);
            }
            // tempos de execução
			fp = fopen("temposbruteforce.txt", "a");
            fprintf(fp, "%d\t%f\n", n, cpu_time);
            fclose(fp);
            break;
		
		case 2:
			fp = fopen("temposbranchandbound.txt", "a");
            fprintf(fp, "%d\t%f\n", n, cpu_time);
            fclose(fp);
            break;
        case 3:
            fp = fopen("tempospermutacoesaleatorias.txt", "a");
            fprintf(fp, "%d\t%f\n", n, cpu_time);
            fclose(fp);
            break;
            
        case 4:
			if (n==13) {
				fp = fopen("randomHistogram13.txt", "w");
				for (int i = 0; i < n*t_range; i++)
				{
					fprintf(fp, "%d\t%d\n", i, histogram[i]);
				}
				fclose(fp);
			}
            if (n==14) {
                fp = fopen("randomHistogram14.txt", "w");
                for (int i = 0; i < n*t_range; i++)
                {
                    fprintf(fp, "%d\t%d\n", i, histogram[i]);
                }
                fclose(fp);
            }
        case 5:
			if (n==13) {
				fp = fopen("BaBHistogram13.txt", "w");
				for (int i = 0; i < n*t_range; i++)
				{
					fprintf(fp, "%d\t%d\n", i, histogram[i]);
				}
				fclose(fp);
			}
            if (n==14) {
                fp = fopen("BaBHistogram14.txt", "w");
                for (int i = 0; i < n*t_range; i++)
                {
                    fprintf(fp, "%d\t%d\n", i, histogram[i]);
                }
                fclose(fp);
            }
		case 6:
			fp = fopen("minimosMaximosNVisitasBF.txt", "a");
			fprintf(fp, "Para n = %d\n", n );
			fprintf(fp, "Min: %d\n", min_cost);
			fprintf(fp, "Max: %d\n", max_cost);
			fprintf(fp, "n_visited: %ld\n", n_visited);
			fclose(fp);
			break;
		case 7: 
			fp = fopen("minimosMaximosNVisitasRP.txt", "a");
			fprintf(fp, "Para n = %d\n", n );
			fprintf(fp, "Min: %d\n", min_cost);
			fprintf(fp, "Max: %d\n", max_cost);
			fprintf(fp, "n_visited: %ld\n", n_visited);
			fclose(fp);
			break;
		case 8:
			fp = fopen("minimosMaximosNVisitasBaB.txt", "a");
			fprintf(fp, "Para n = %d\n", n );
			fprintf(fp, "Min: %d\n", min_cost);
			fprintf(fp, "Max: %d\n", max_cost);
			fprintf(fp, "n_visited: %ld\n", n_visited);
			fclose(fp);
			break;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// main program
//

int main(int argc,char **argv)
{
  if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'e')
  {
    seed = 0;
    {
      int n = 3;
      init_costs(-3); // costs for the example with n = 3
      int a[n];
      for(int i = 0;i < n;i++)
        a[i] = i;
      reset_solutions();
      (void)elapsed_time();
      generate_all_permutations(n,0,a);
      cpu_time = elapsed_time();
      show_solutions(n,"Example for n=3",show_all);
      printf("\n");
    }
    {
      int n = 5;
      init_costs(-5); // costs for the example with n = 5
      int a[n];
      for(int i = 0;i < n;i++)
        a[i] = i;
      reset_solutions();
      (void)elapsed_time();
      generate_all_permutations(n,0,a);
      cpu_time = elapsed_time();
      show_solutions(n,"Example for n=5",show_all);
      return 0;
    }
  }
  if(argc == 2)
  {
    seed = atoi(argv[1]); // seed = student number
    if(seed >= 0 && seed <= 1000000)
    {
      for(int n = 1;n <= max_n;n++)
      {
        init_costs(n);
        show_solutions(n,"Problem statement",show_info_1 | show_costs);
        //
        // 2.
        //
        if(n <= 14) // use a smaller limit here while developing your code(14)
        {
          int a[n];
          for(int i = 0;i < n;i++)
            a[i] = i; // initial permutation
          reset_solutions();
          (void)elapsed_time();
          generate_all_permutations(n,0,a);
          cpu_time = elapsed_time();
          show_solutions(n,"Brute force",show_info_2 | show_min_solution | show_max_solution | show_histogram);
          print_to_file(n, 1);
          print_to_file(n, 6);
        }
        //
        // place here your code that solves the problem with branch-and-bound
        //
        //less_row_cost = +1000000000;

        if(n <= 14) // use a smaller limit here while developing your code(16)
        {
          int a[n];
          for(int i = 0;i < n;i++)
            a[i] = i; // initial permutation
          reset_solutions();
          (void)elapsed_time();
          generate_all_permutations_branch_and_bound(n,0,a,0);
          cpu_time = elapsed_time();
          show_solutions(n,"Brute force with branch-and-bound",show_info_2 | show_min_solution);
          print_to_file(n, 2);
          print_to_file(n, 5);
          print_to_file(n, 8);
        }

        //
        // place here your code that generates the random permutations
        //

		if (n <= 14)
		{
			int t[n];
			reset_solutions();
			(void)elapsed_time();
			calculate_solutions(n, t);
			cpu_time = elapsed_time();
			show_solutions(n,"Random Permutation",show_info_2 | show_min_solution | show_max_solution);
			print_to_file(n, 3);
			print_to_file(n, 4);
			print_to_file(n, 7);
		}
		
        //
        // done
        //
        printf("\n");
      }
      return 0;
    }
  }
  fprintf(stderr,"usage: %s -e              # for the examples\n",argv[0]);
  fprintf(stderr,"usage: %s student_number\n",argv[0]);
  return 1;
}
