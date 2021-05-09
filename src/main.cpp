#include "mpi.h"
#include "IPGlasma.h"

#include <iostream>

#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
using namespace std;

// main program 1
int main(int argc, char *argv[]) {
    int rank;
    int size;

    int nev = 1;
    if (argc == 3) {
        nev = atoi(argv[2]);
    }

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes

    IPGlasma IPGgenerator(rank, size, nev, "input");

    // event loop starts ...
    for (int iev = 0; iev < nev; iev++) {
        IPGgenerator.generateAnEvent(iev);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (IPGgenerator.get_h5Flag() && rank == 0) {
        int status = 0;
        stringstream collect_command;
        collect_command << "python3 utilities/combine_events_into_hdf5.py ."
                        << " --output_filename RESULTS"
                        << " --combine_hdf5_files_only";
        status = system(collect_command.str().c_str());
        std::cout << "finished system call to python script with status: "
                  << status << std::endl;
    }
    MPI_Finalize();
    return 1;
}


