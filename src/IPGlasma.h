
#ifndef IPGlasma_H_
#define IPGlasma_H_

#include "pretty_ostream.h"
#include "Parameters.h"
#include "Setup.h"
#include "Random.h"

#include <string>

class IPGlasma {
 private:
    int rank_;
    int size_;
    int nev_;

    int h5Flag_;

    Parameters *param;
    Random* random;
    pretty_ostream messager;

 public:
    IPGlasma(int rank, int size, int nev, std::string inputFilename);
    ~IPGlasma();

    void readInput(Setup *setup, std::string inputFilename);

    void generateAnEvent(int iev);

};


#endif  // IPGlasma_H_
