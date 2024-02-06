#include "NBodySimulation.h"
# include <omp.h>

// TODO: Try using structs istead

class NBodySimulationVectorised : public NBodySimulation {
  public:
  #pragma omp declare simd
  double force_calculation (int i, int j, int direction){
    // Euclidean distance
    const double distance = sqrt(
                                 (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                                 (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                                 (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                                 );
    const double distance3 = distance * distance * distance;
    minDx = std::min( minDx,distance );

    return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;
  }

  void updateBody () {
    timeStepCounter++;
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();

    // force0 = force along x direction
    // force1 = force along y direction
    // force2 = force along z direction
    double* force0 = new double[NumberOfBodies];
    double* force1 = new double[NumberOfBodies];
    double* force2 = new double[NumberOfBodies];


    // Initialise forces to 0 
    for (int i = 0; i < NumberOfBodies; i++) {
          force0[i] = 0.0;
          force1[i] = 0.0;
          force2[i] = 0.0;
    }


    #pragma omp parallel for simd reduction(+:force0[:NumberOfBodies], force1[:NumberOfBodies], force2[:NumberOfBodies]) reduction(min:minDx)
    for (int i=0; i<NumberOfBodies; i++) {
      //pragma omp simd 
      for (int j = i+1; j < NumberOfBodies; j++) {
        double f0, f1, f2;
        f0 = force_calculation(j,i,0);    // Calculate force between i and j
        f1 = force_calculation(j,i,1);
        f2 = force_calculation(j,i,2);

        // Force acting i by j
        force0[i] += f0;
        force1[i] += f1;
        force2[i] += f2;
        // Force acting on j by i, opposite magnitude hence negative sign
        force0[j] += -f0;
        force1[j] += -f1;
        force2[j] += -f2;
      }
    }

    // Update velocity and position  
    #pragma omp parallel for simd reduction(max:maxV)
    for (int i = 0; i < NumberOfBodies; i++) {
      //#pragma omp simd
      for (int dim = 0; dim < 3; dim++) {
        x[i][dim] = x[i][dim] + timeStepSize * v[i][dim];
        v[i][dim] = v[i][dim] + timeStepSize * force0[i] / mass[i];
      }
      maxV = std::max(std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]), maxV);
    }


    
    t += timeStepSize;

    delete[] force0;
    delete[] force1;
    delete[] force2;
  }
};
