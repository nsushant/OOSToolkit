
#include <armadillo>
#include <iostream>
#include <string>

#include "Nbody.hpp"



int main(int argc, char *argv[])
{
  std::cout << "------------------Running Simulation---------------------" << std::endl;

  force_model fmodel(true,false);

  run_simulation("WalkerDelta.csv", 31536000, fmodel);

  std::cout << "finished running sim" << "\n";

  return 0;
}
