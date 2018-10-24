#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <assert.h>
#include <time.h>
#include "mpi.h"
#include "wire_route.h"

int dim_y;
int dim_x;
int global_delta;
int global_numWires;
const int root = 0;
wire_t* wires;
bend_t* bends;
int* costs;
int* trans_costs;
const int INT_MAX = 0x7fffffff;
double parallelTime;
double waitTIme;
double otherTime;


// Initialize problem
static inline void init(int numRows, int numCols, int delta, int numWires)
{
    dim_y = numRows;
    dim_x = numCols;
    global_delta = delta;
    global_numWires = numWires;
    wires = (wire_t*) malloc(sizeof(wire_t) * numWires);
    bends = (bend_t*) malloc(sizeof(bend_t) * numWires);
    costs = (int*) malloc(sizeof(int) * numRows * numCols);
    trans_costs = (int*) malloc(sizeof(int) * numRows * numCols);
}

// Initialize a given wire
static inline void initWire(int wireIndex, int x1, int y1, int x2, int y2)
{
    wires[wireIndex].x1 = x1;
    wires[wireIndex].y1 = y1;
    wires[wireIndex].x2 = x2;
    wires[wireIndex].y2 = y2;
}

// Return number of rows
static inline int getNumRows()
{
    return dim_y;
}

// Return number of cols
static inline int getNumCols(){
    return dim_x;
}

// Return delta
static inline int getDelta()
{
    return global_delta;
}

// Return number of wires
static inline int getNumWires()
{
    return global_numWires;
}

// Get costs array entry
static inline int getcosts(int row, int col)
{
    return costs[row * dim_x + col];
}

// Get a wire placement. Returns number of points (should be 2-4 for 0-2 bends).
static inline int getWire(int wireIndex, int* x1, int* y1, int* x2, int* y2, int* x3, int* y3, int* x4, int* y4)
{
    int x_start = wires[wireIndex].x1;
    int y_start = wires[wireIndex].y1;
    int x_end = wires[wireIndex].x2;
    int y_end = wires[wireIndex].y2;
    int x_bend = bends[wireIndex].x;
    int y_bend = bends[wireIndex].y;

    *x1 = x_start;
    *y1 = y_start;
    *x3 = -1;
    *y3 = -1;
    *x4 = -1;
    *y4 = -1;
    if ((x_start == x_end && x_start == x_bend) || (y_start == y_end && y_start == y_bend)) {
        *x2 = x_end;
        *y2 = y_end;
        return 2;
    }

    *x2 = x_bend;
    *y2 = y_bend;

    if (x_bend == x_end || y_bend == y_end) {    
        *x3 = x_end;
        *y3 = y_end;
        return 3;
    }

    if (x_bend == x_start) {
        *x3 = x_end;
        *y3 = y_bend;     
    } else if (y_bend == y_start) {
        *x3 = x_bend;
        *y3 = y_end;  
    }

    *x4 = x_end;
    *y4 = y_end;   
    //printf("xbend: %d ybend: %d\n", x_bend, y_bend);
    //printf("x1: %d y1: %d x2: %d y2: %d x3: %d y3: %d x4: %d y4: %d \n", *x1, *y1, *x2, *y2, *x3, *y3, *x4, *y4 );
    return 4;
}




void remove_wire(int x1, int x2, int y1, int y2, bend_t bend) {
    double start = MPI_Wtime();

    int bendx = bend.x;
    int bendy = bend.y;
    int unit;
    // horizontal first
    if (bendy == y1) {
        unit = x1 < bendx? 1 : -1;
        
        for (int x = x1; x != bendx; x+=unit) {
            costs[x + y1 * dim_x] -= 1;
            trans_costs[x * dim_y + y1] -= 1;
        }
        unit = y1 < y2? 1 : -1;
        for (int y = y1; y != y2; y+=unit) {
            costs[bendx+y*dim_x] -= 1;
            trans_costs[bendx * dim_y + y] -= 1;
        }
        unit = bendx < x2? 1 : -1;
        for (int x = bendx; x != x2+unit; x+=unit) {
            costs[x+y2*dim_x] -= 1;
            trans_costs[x * dim_y + y2] -= 1;
        }
    }
    // vertical first
    else {
        unit = y1 < bendy? 1 : -1;
        for (int y = y1; y != bendy; y+=unit) {
            costs[x1+y*dim_x] -= 1;
            trans_costs[x1 * dim_y + y] -= 1;
        }
        unit = x1 < x2? 1 : -1;
        for (int x = x1; x != x2; x+=unit) {
            costs[x+bendy*dim_x] -= 1;
            trans_costs[x * dim_y + bendy] -= 1;
        }
        unit = bendy < y2? 1 : -1;
        for (int y = bendy; y != y2+unit; y+=unit) {
            costs[x2+y*dim_x] -= 1;
            trans_costs[x2 * dim_y + y] -= 1;
        }
    }
    otherTime += MPI_Wtime() - start;
}


void add_wire(int x1, int x2, int y1, int y2, bend_t bend) {
    double start = MPI_Wtime();
    int bendx = bend.x;
    int bendy = bend.y;
    int unit;
    // horizontal first
    if (bendy == y1) {
        unit = x1 < bendx? 1 : -1;
        
        for (int x = x1; x != bendx; x+=unit) {
            costs[x + y1 * dim_x] += 1;
            trans_costs[x * dim_y + y1] += 1;
        }
        unit = y1 < y2? 1 : -1;
        for (int y = y1; y != y2; y+=unit) {
            costs[bendx+y*dim_x] += 1;
            trans_costs[bendx * dim_y + y] += 1;
        }
        unit = bendx < x2? 1 : -1;
        for (int x = bendx; x != x2+unit; x+=unit) {
            costs[x+y2*dim_x] += 1;
            trans_costs[x * dim_y + y2] += 1;
        }
    }
    // vertical first
    else {
        unit = y1 < bendy? 1 : -1;
        for (int y = y1; y != bendy; y+=unit) {
            costs[x1+y*dim_x] += 1;
            trans_costs[x1 * dim_y + y] += 1;
        }
        unit = x1 < x2? 1 : -1;
        for (int x = x1; x != x2; x+=unit) {
            costs[x+bendy*dim_x] += 1;
            trans_costs[x * dim_y + bendy] += 1;
        }
        unit = bendy < y2? 1 : -1;
        for (int y = bendy; y != y2+unit; y+=unit) {
            costs[x2+y*dim_x] += 1;
            trans_costs[x2 * dim_y + y] += 1;
        }
    }
    otherTime += MPI_Wtime() - start;
}

void calculate_cost(int x1, int x2, int y1, int y2, int bendx, int bendy, 
    int *max_cost, int *cost_sum) {

  int unit;
  int local_max_cost = 0;
  int local_cost_sum = 0;
  int new_cost;

  // horizontal first
  if (bendy == y1) {
    unit = x1 < bendx? 1 : -1;
    for (int x = x1; x != bendx; x+=unit) {
      new_cost = costs[x+y1*dim_x] + 1;
      local_max_cost = new_cost > local_max_cost? new_cost : local_max_cost;
      local_cost_sum += new_cost*new_cost;
    }

    unit = y1 < y2? 1 : -1;
    for (int y = y1; y != y2; y+=unit) {
      new_cost = trans_costs[bendx*dim_y+y] + 1;
      local_max_cost = new_cost > local_max_cost? new_cost : local_max_cost;
      local_cost_sum += new_cost*new_cost;
    }

    unit = bendx < x2? 1 : -1;
    for (int x = bendx; x != x2+unit; x+=unit) {
      new_cost = costs[x+y2*dim_x] + 1;
      local_max_cost = new_cost > local_max_cost? new_cost : local_max_cost;
      local_cost_sum += new_cost*new_cost;
    }
  }
  // vertical first
  else {
    unit = y1 < bendy? 1 : -1;
    for (int y = y1; y != bendy; y+=unit) {
      new_cost = trans_costs[dim_y*x1+y] + 1;
      local_max_cost = new_cost > local_max_cost? new_cost : local_max_cost;
      local_cost_sum += new_cost*new_cost;
    }

    unit = x1 < x2? 1 : -1;
    for (int x = x1; x != x2; x+=unit) {
      new_cost = costs[x+bendy*dim_x] + 1;
      local_max_cost = new_cost > local_max_cost? new_cost : local_max_cost;
      local_cost_sum += new_cost*new_cost;
    }

    unit = bendy < y2? 1 : -1;
    for (int y = bendy; y != y2+unit; y+=unit) {
      new_cost = trans_costs[dim_y*x2+y] + 1;
      local_max_cost = new_cost > local_max_cost? new_cost : local_max_cost;
      local_cost_sum += new_cost*new_cost;
    }
  }

  *cost_sum = local_cost_sum;
  *max_cost = local_max_cost;
}

void iterations(float SA_prob, int SA_iters, int nproc, int procID) {
   
    int half_delta = global_delta/2;
    int randMax = 1000;

    //Arrays for result combinations
    int* received = (int*) malloc(4 * sizeof(int) * nproc);

    // go over all wires
    srand((unsigned)time(0));
    for (int it = 0; it < SA_iters; it++) {
        for (int i = 0; i < global_numWires; i++) {
            int x1 = wires[i].x1;
            int x2 = wires[i].x2;
            int y1 = wires[i].y1;
            int y2 = wires[i].y2;
            
            if (it > 0) {
                remove_wire(x1, x2, y1, y2, bends[i]);
            }
            
            int randP;
            if (procID == root) {
                randP = rand() % randMax;
            }
            MPI_Bcast(&randP, 1, MPI_INT, root, MPI_COMM_WORLD);
            
            if (randP >= randMax * SA_prob) {
                int low_hor;
                int high_hor;
                int low_ver;
                int high_ver;
                // horizontal
                if (x1 < x2) {
                    low_hor = x1 - half_delta >= 0 ? x1 - half_delta : 0;
                    high_hor = x2 + half_delta + 1 < dim_x ? x2 + half_delta + 1 : dim_x - 1;
                } else {
                    low_hor = x2 - half_delta >= 0 ? x2 - half_delta : 0;
                    high_hor = x1 + half_delta + 1 < dim_x ? x1 + half_delta + 1 : dim_x - 1;
                }
                
                // vertical
                if (y1 < y2) {
                    low_ver = y1 - half_delta >= 0 ? y1 - half_delta : 0;
                    high_ver = y2 + half_delta + 1 < dim_y ? y2 + half_delta + 1 : dim_y - 1;
                } else {
                    low_ver = y2 - half_delta >= 0 ? y2 - half_delta : 0;
                    high_ver = y1 + half_delta + 1 < dim_y ? y1 + half_delta + 1 : dim_y - 1;
                }
                int numRoutes = high_hor - low_hor + high_ver - low_ver;

                int startIndex = procID;
    
                int min_max_cost = INT_MAX;
                int min_cost_sum = INT_MAX;
                int max_cost = -1;
                int cost_sum = -1;
                int best_r = -1;
                int best_c = -1;
                double parallelStart = MPI_Wtime();
                for (int route = startIndex; route < numRoutes; route += nproc) {
                    int r;
                    int c;
                    if (route < high_hor - low_hor) {
                        c = route + low_hor;
                        r = y1;
                        if (c == x1) {
                            continue;
                        }
                        calculate_cost(x1, x2, y1, y2, c, y1,
                                       &max_cost, &cost_sum);
                    } 
                    else {
                        r = route - (high_hor - low_hor) + low_ver;
                        c = x1;
                        if (r == y1) {
                            continue;
                        }
                        calculate_cost(x1, x2, y1, y2, x1, r,
                                       &max_cost, &cost_sum);
                    }
                    //printf("r:%d c:%d max_cost: %d cost_sum: %d\n", r, c, max_cost, cost_sum );
                    if (max_cost < min_max_cost || (max_cost == min_max_cost && cost_sum < min_cost_sum)) {
                        min_max_cost = max_cost;
                        min_cost_sum = cost_sum;
                        best_r = r;
                        best_c = c;
                    }
                }
                double parallelEnd = MPI_Wtime();
                parallelTime += parallelEnd - parallelStart;
                // root combining results
                int sent[4];
                sent[0] = min_max_cost;
                sent[1] = min_cost_sum;
                sent[2] = best_r;
                sent[3] = best_c;

                double waitStart = MPI_Wtime();
                MPI_Gather(sent, 4, MPI_INT, received, 4, MPI_INT, root, MPI_COMM_WORLD);
                double waitEnd = MPI_Wtime();
                waitTIme += waitEnd - waitStart;

                if (procID == root) {
                    for (int j = 0; j < nproc; j++) {
                        max_cost = received[4 * j];
                        cost_sum = received[4 * j + 1];
                        if (max_cost < min_max_cost || (max_cost == min_max_cost && cost_sum < min_cost_sum)) {
                            min_max_cost = max_cost;
                            min_cost_sum = cost_sum;
                            best_r = received[4 * j + 2];
                            best_c = received[4 * j + 3];
                        }
                    }
                    bends[i].x = best_c;
                    bends[i].y = best_r;

                }
            } else {
                if (procID == root) {
                    int lox;
                    int loy;
                    int hix;
                    int hiy;
                    if (x1 < x2) {
                        lox = x1 - half_delta >= 0 ? x1 - half_delta : 0;
                        hix = x2 + half_delta + 1 <= dim_x ? x2 + half_delta + 1: dim_x;
                    } else {
                        lox = x2 - half_delta >= 0 ? x2 - half_delta : 0;
                        hix = x1 + half_delta + 1<= dim_x ? x1 + half_delta + 1: dim_x;
                    }

                    if (y1 < y2) {
                        loy = y1 - half_delta >= 0 ? y1 - half_delta : 0;
                        hiy = y2 + half_delta + 1<= dim_y ? y2 + half_delta + 1: dim_y;
                    } else {
                        loy = y2 - half_delta >= 0 ? y2 - half_delta : 0;
                        hiy = y1 + half_delta + 1 <= dim_y ? y1 + half_delta + 1 : dim_y;
                    }
                    while (1) {
                        //printf("%d %d %d %d %d %d %d %d denom: %d \n", x1, x2, y1, y2, hix, hiy, lox, loy, (hix + hiy - lox - loy - 2));
                        int randBend = rand() % (hix + hiy - lox - loy);
                        
                        // Random among all possible hix + hiy - lox - loy - 2 options
                        // (excluding the bending at the start point).
                        if (randBend < hix - lox) {
                            bends[i].x = lox + randBend;
                            bends[i].y = y1;
                        } else {
                            int index = randBend - (hix - lox);
                            bends[i].x = x1;
                            bends[i].y = loy + index;
                        }

                        if (!(bends[i].x == x1 && bends[i].y == y1)) {
                            break;
                        }
                    }
                }
            }

            MPI_Bcast(&bends[i], 2, MPI_INT, root, MPI_COMM_WORLD);
            add_wire(x1, x2, y1, y2, bends[i]);

        }
    }
}


// Read input file
void readInput(char* inputFilename)
{
    FILE* fp = fopen(inputFilename, "r");
    int numRows;
    int numCols;
    int delta;
    int numWires;
    int wireIndex;
    int x1;
    int y1;
    int x2;
    int y2;
    if (fp == NULL) {
        fprintf(stderr, "Failed to read input file %s\n", inputFilename);
        exit(-1);
    }
    if (fscanf(fp, "%d %d", &numCols, &numRows) != 2) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if ((numRows <= 0) || (numCols <= 0)) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (fscanf(fp, "%d", &delta) != 1) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (fscanf(fp, "%d", &numWires) != 1) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (numWires <= 0) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    init(numRows, numCols, delta, numWires);
    for (wireIndex = 0; wireIndex < numWires; wireIndex++) {
        if (fscanf(fp, "%d %d %d %d", &x1, &y1, &x2, &y2) != 4) {
            fprintf(stderr, "Invalid input file format\n");
            exit(-1);
        }
        initWire(wireIndex, x1, y1, x2, y2);
    }
    fclose(fp);
}

// Write costs array file based on input filename
void writecosts(char* inputFilename, int nproc)
{
    char* dname = strdup(inputFilename);
    char* bname = strdup(inputFilename);
    char* costsFilename = malloc(strlen(inputFilename) + 100);
    FILE* fp;
    int numRows = getNumRows();
    int numCols = getNumCols();
    int row;
    int col;
    assert(costsFilename != NULL);
    sprintf(costsFilename, "%s/costs_%s_%d.txt", dirname(dname), basename(bname), nproc);
    fp = fopen(costsFilename, "w");
    printf("file name: %s\n", costsFilename);
    if (fp == NULL) {
        fprintf(stderr, "Failed to write c file %s\n", costsFilename);
        exit(-1);
    }
    fprintf(fp, "%d %d\n", numCols, numRows);
    for (row = 0; row < numRows; row++) {
        for (col = 0; col < numCols; col++) {
            fprintf(fp, "%d ", getcosts(row, col));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    free(costsFilename);
    free(bname);
    free(dname);
}

// Write wire output file based on input filename
void writeOutput(char* inputFilename, int nproc)
{
    char* dname = strdup(inputFilename);
    char* bname = strdup(inputFilename);
    char* outputFilename = malloc(strlen(inputFilename) + 100);
    FILE* fp;
    int numRows = getNumRows();
    int numCols = getNumCols();
    int delta = getDelta();
    int numWires = getNumWires();
    int wireIndex;
    int numPoints;
    int x1;
    int y1;
    int x2;
    int y2;
    int x3;
    int y3;
    int x4;
    int y4;
    assert(outputFilename != NULL);
    sprintf(outputFilename, "%s/output_%s_%d.txt", dirname(dname), basename(bname), nproc);
    fp = fopen(outputFilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write output file %s\n", outputFilename);
        exit(-1);
    }
    fprintf(fp, "%d %d\n", numCols, numRows);
    fprintf(fp, "%d\n", delta);
    fprintf(fp, "%d\n", numWires);
    for (wireIndex = 0; wireIndex < numWires; wireIndex++) {
        numPoints = getWire(wireIndex, &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4);
        switch (numPoints) {
            case 2:
                fprintf(fp, "%d %d %d %d\n", x1, y1, x2, y2);
                break;

            case 3:
                fprintf(fp, "%d %d %d %d %d %d\n", x1, y1, x2, y2, x3, y3);
                break;

            case 4:
                fprintf(fp, "%d %d %d %d %d %d %d %d\n", x1, y1, x2, y2, x3, y3, x4, y4);
                break;

            default:
                assert(0); // invalid number of points
                break;
        }
    }
    fclose(fp);
    free(outputFilename);
    free(bname);
    free(dname);
}

void print_result(){
    int max_cost = 0;
    int sum_cost = 0;
    for (int row = 0; row < dim_y; row += 1) {
        for (int col = 0; col < dim_x; col += 1) {
            max_cost = costs[row * dim_x + col] > max_cost ? costs[row * dim_x + col] : max_cost;
            sum_cost += costs[row * dim_x + col] * costs[row * dim_x + col];
        }
    }
    printf("Max cost: %d\n", max_cost);
    printf("Sum cost: %d\n", sum_cost);
}

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations)
{
    readInput(inputFilename);
    parallelTime = 0;
    waitTIme = 0;
    otherTime = 0;
    iterations(prob, numIterations, nproc, procID);
    
    if (procID == root) {
        writecosts(inputFilename, nproc);
        writeOutput(inputFilename, nproc);
        print_result();
    }

    printf("procID: %d parallelTime: %f\n", procID, parallelTime);
    printf("procID: %d waitTIme: %f\n", procID, waitTIme);
    printf("procID: %d otherTime: %f\n", procID, otherTime);
}
