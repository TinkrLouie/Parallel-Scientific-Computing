#include <iomanip>

#include "NBodySimulation.h"

/**
 * You can compile this file with
 *   make step-2-g++   // Uses the GNU Compiler Collection.
 *   make step-2-icpx  // Uses the Intel compiler.
 * and run it with
 *   ./step-2-g++
 *   ./step-2-icpx
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */


class NBodySimulationMolecularForces : public NBodySimulation {
  
  public:
  double distance (int i, int j) {  // Euclidean distance between bodies i and j
    const double distance = sqrt(
      (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
      (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
      (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
      );
    return distance;
  }

  double acceleration (int i, int j, double distance, int direction){ // acceleration of body i through body j
    const double distance3 = distance * distance * distance;
    minDx = std::min( minDx,distance );  // can be taken out

    return (x[i][direction]-x[j][direction]) * mass[i] / distance3;
  }

  void updateBody () {
    timeStepCounter++;
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();
    double c = 0.01/NumberOfBodies;
    double dist, nr = 1.0/6;
    double* x2 = new double[3];  // second order x,y,z
    double* x3 = new double[3];  // third order x,y,z
    double* x4 = new double[3];  // fourth order x,y,z
    double* v2 = new double[3];
    double* v3 = new double[3];
    double* v4 = new double[3];
    double* a1 = new double[3];
    double* a2 = new double[3];
    double* a3 = new double[3];
    double* a4 = new double[3];
    double* d = new double[3];

    for (int i=0; i<NumberOfBodies; i++) {
      for (int j = 0; j < NumberOfBodies; j++) {
        if (i == j) continue;
        
        // Step 1
        dist = distance(i,j);
        for (int dim = 0; dim < 3; dim++) a1[dim] = acceleration(j,i,dist,dim);

        // Step 2
        for (int dim = 0; dim < 3; dim++) v2[dim] = v[i][dim] + a1[dim]*timeStepSize*0.5;  // compute 2nd order v
        for (int dim = 0; dim < 3; dim++) x2[dim] = x[i][dim] + v[i][dim]*timeStepSize*0.5; // compute 2nd order x
        for (int dim = 0; dim < 3; dim++) d[dim] = x[j][dim] - x2[dim];  // compute dx,dy,dz of 2nd order x
        dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 2nd order x to j
        for (int dim = 0; dim < 3; dim++) a2[dim] = d[dim]*mass[j]/(dist*dist*dist);  // 2nd order acceleration of i

        // Step 3
        for (int dim = 0; dim < 3; dim++) v3[dim] = v[i][dim] + a2[dim]*timeStepSize*0.5;  // compute 3rd order v
        for (int dim = 0; dim < 3; dim++) x3[dim] = x[i][dim] + v2[dim]*timeStepSize*0.5; // compute 3rd order x
        for (int dim = 0; dim < 3; dim++) d[dim] = x[j][dim] - x3[dim];  // compute dx,dy,dz of 3rd order x
        dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 3rd order x to j
        for (int dim = 0; dim < 3; dim++) a3[dim] = d[dim]*mass[j]/(dist*dist*dist);  // 3rd order acceleration of i

        // Step 4
        for (int dim = 0; dim < 3; dim++) v4[dim] = v[i][dim] + a3[dim]*timeStepSize*0.5;  // compute 4th order v
        for (int dim = 0; dim < 3; dim++) x4[dim] = x[i][dim] + v3[dim]*timeStepSize*0.5; // compute 4th order x
        for (int dim = 0; dim < 3; dim++) d[dim] = x[j][dim] - x4[dim];  // compute dx,dy,dz of 4th order x
        dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 4th order x to j
        for (int dim = 0; dim < 3; dim++) a4[dim] = d[dim]*mass[j]/(dist*dist*dist);  // 4th order acceleration of i


        for (int dim = 0; dim < 3; dim++) x[i][dim] = x[i][dim] + nr*(v[i][dim] + 2*v2[dim] + 2*v3[dim] + v4[dim]) * timeStepSize;
        for (int dim = 0; dim < 3; dim++) v[i][dim] = v[i][dim] + nr*(a1[dim] + 2*a2[dim] + 2*a3[dim] + a4[dim]) * timeStepSize;

        
        maxV = std::max(std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]), maxV);
      }
    }
    t += timeStepSize;
    delete[] x2;
    delete[] x3;
    delete[] x4;
    delete[] v2;
    delete[] v3;
    delete[] v4;
    delete[] a1;
    delete[] a2;
    delete[] a3;
    delete[] a4;
    delete[] d ;
  }
};


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
  NBodySimulationMolecularForces nbs;
  nbs.setUp(argc,argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}
