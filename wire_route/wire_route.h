#ifndef _WIRE_ROUTE_H
#define _WIRE_ROUTE_H

typedef struct
{
  /* Define the data structure for wire here */
  int x1;
  int x2;
  int y1;
  int y2;
} wire_t;

typedef struct
{
	int x;
	int y;
} bend_t;

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations);

// Read input file
void readInput(char* inputFilename);

// Write cost array file based on input filename
void writeCost(char* inputFilename, int nproc);

// Write wire output file based on input filename
void writeOutput(char* inputFilename, int nproc);

#endif // _WIRE_ROUTE_H
