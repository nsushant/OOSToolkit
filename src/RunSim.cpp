
#include <armadillo>
#include <iostream>
#include <string>

#include "Nbody.hpp"



int main(int argc, char *argv[])
{

  force_model fmodel(true,false);

  std::cout << " Running Simulation with params ----- " << "\n";
  std::cout << " - J2 : "<< fmodel.includeJ2 << "(1 or 0)"<<" \n ";
  std::cout << "- Mutual : "<< fmodel.includeMutual << " (1 or 0) "<<"\n ";
  std::cout << "- dt : "<< tstep_size << " seconds "<<"\n ";
  std::cout << "- t_final : "<< t_final << " seconds"<<" \n";
  std::cout << "------------------------------------- " <<"\n";

  run_simulation("WalkerDelta.csv", t_final, fmodel);

  std::cout << "finished running sim" << "\n";

  return 0;
}
