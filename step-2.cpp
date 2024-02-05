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
  void updateBody () {
    int numBuckets = 10;
    double vBucket = maxV / numBuckets;
    timeStepCounter++;
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();
    double c = 0.01/NumberOfBodies;
    // force0 = force along x direction
    // force1 = force along y direction
    // force2 = force along z direction
    double* force0 = new double[NumberOfBodies];
    double* force1 = new double[NumberOfBodies];
    double* force2 = new double[NumberOfBodies];

    int* particleInBucket = new int[NumberOfBodies]();
    for (int b = 1; b < numBuckets+1; b++) {
      for (int p = 0; p < NumberOfBodies; p++) {
        if (particleInBucket[p]) continue;
        else if (v[p][0]*v[p][0] + v[p][1]*v[p][1] + v[p][2]*v[p][2] <
	         (b+1)*vBucket*vBucket) particleInBucket[p] = b;
      }
    }
    // Initialise forces to 0 
    for (int i = 0; i < NumberOfBodies; i++) {
          force0[i] = 0.0;
          force1[i] = 0.0;
          force2[i] = 0.0;
    }

    for (int bucket = 1; bucket < numBuckets+1; bucket++) {
    int i = 0;
    double oldTimeStepSize = timeStepSize;
    timeStepSize /= std::pow(2, bucket - 1);
    // find all the particles belonging to the bucket
    for (int particle = i; particle < NumberOfBodies; particle++) {
      if (particleInBucket[particle] != bucket) {
	      continue;
      }
      i = particle;
      // split dt into particleInBucket[j] portions, then perform them
      for (int iter = 0; iter < std::pow(2, bucket-1); iter++) {
	      // zero force every partial timestep
	      force0[i] = 0.0;
        force1[i] = 0.0;
        force2[i] = 0.0;
	      for (int j = 0; j < NumberOfBodies; j++) {
	        if (i == j) continue;
          force0[i] += force_calculation(j,i,0);    // force0[0] to [i]
          force1[i] += force_calculation(j,i,1);
          force2[i] += force_calculation(j,i,2);

          double distance = sqrt(
                                 (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                                 (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                                 (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                                 );

	        // check for collisions
	        if (distance <= c*(mass[i] + mass[j])){
            // Momentum calculations
            x[i][0] = (mass[i]*x[i][0] + mass[j]*x[j][0]) / (mass[i]+mass[j]);
            x[i][1] = (mass[i]*x[i][1] + mass[j]*x[j][1]) / (mass[i]+mass[j]);
            x[i][2] = (mass[i]*x[i][2] + mass[j]*x[j][2]) / (mass[i]+mass[j]);

            v[i][0] = (mass[i]*v[i][0] + mass[j]*v[j][0]) / (mass[i]+mass[j]);
            v[i][1] = (mass[i]*v[i][1] + mass[j]*v[j][1]) / (mass[i]+mass[j]);
            v[i][2] = (mass[i]*v[i][2] + mass[j]*v[j][2]) / (mass[i]+mass[j]);

            // Mass of merged object
            mass[i] += mass[j];


            const int l = --NumberOfBodies;
            // if (NumberOfBodies < 2) {     // Print summary and exit if merge is between last 2 bodies
	          // std::cout << "Two remaining bodies merged." << std::endl;
            // printSummary();
            // closeParaviewVideoFile();
	          // std::exit(0);
            // }
            // Remove other merged object from list
            for (int dim = 0; dim < 3; dim++) {
	          x[j][dim] = x[l][dim];
	          v[j][dim] = v[l][dim];
            } 
            force0[j] = force0[l];
            force1[j] = force1[l];
            force2[j] = force2[l];
            mass[j] = mass[l];
            j--;
	        }
	      }
	      for (int i = 0; i < NumberOfBodies; i++) {
          x[i][0] = x[i][0] + timeStepSize * v[i][0];
          x[i][1] = x[i][1] + timeStepSize * v[i][1];
          x[i][2] = x[i][2] + timeStepSize * v[i][2];

          v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
          v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
          v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

          maxV = std::max(std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]), maxV);
        }
    }
    }
    timeStepSize = oldTimeStepSize;
    }
    t += timeStepSize;

    delete[] force0;
    delete[] force1;
    delete[] force2;
    delete[] particleInBucket;
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
