#include "typedef.h"
#include "Input.h"
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
  MPI_Init(NULL, NULL);
  MPI_Comm comm = MPI_COMM_WORLD;

  UInt rank;
  MPI_Comm_rank(comm, (int*)&rank );
  UInt nprocs;
  MPI_Comm_size(comm, (int*)&nprocs);

  Input input;
  input.openFile( argv[1] );

  for( UInt p=0; p<input.getNumProblems(); p++ )
  {
    Problem prob;
    prob.read( input, p );
    prob.solve();
  }

  MPI_Finalize();
  return 0;
}
