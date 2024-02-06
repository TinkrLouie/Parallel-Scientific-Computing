#include <iomanip>
#include <chrono>
#include "NBodySimulationVectorised.cpp"

/**
 * You can compile this file with
 *   make step-3-g++   // Uses the GNU Compiler Collection.
 *   make step-3-icpx  // Uses the Intel compiler.
 * and run it with
 *   ./step-3-g++
 *   ./step-3-icpx
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

/**
 * Main routine.
 *
 * No major changes are needed in the assignment. You can add initialisation or
 * or remove input checking, if you feel the need to do so. But keep in mind
 * that you may not alter what the program writes to the standard output.
 */
int main (int argc, char** argv) {
  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulationVectorised nbs;
  nbs.setUp(argc,argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  auto t1 = std::chrono::high_resolution_clock::now();
  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }
  auto t2 = std::chrono::high_resolution_clock::now();

  nbs.printSummary();
  nbs.closeParaviewVideoFile();
  
 
  // Calculating total time taken by the program. 
  std::cout << "Time taken by program is : " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " millisec " << std::endl;
  return 0;
}
