#include "NBodySimulation.h"
# include <omp.h>

// TODO: Try using structs istead

class NBodySimulationVectorised : public NBodySimulation {
  public:
    // Runge-Kutta 4th Order
    void rk4 (int i, int j) {
      double xTemp[4][3];
      double vTemp[4][3];
      double aTemp[4][3];
      double d[3];
      int dim;
      double dist, nr = 1.0/6;
      double c = 1e-2;
      
      // Step 1-------------------------------------------------------------
      dist = sqrt((x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                  (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                  (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                 );
      // Collision detection
      if (dist <= (c/NumberOfBodies)*(mass[i] + mass[j])){
        // Momentum update
        #pragma omp simd private(dim)
        for (dim = 0; dim < 3; dim++) {
          x[i][dim] = (mass[i]*x[i][dim] + mass[j]*x[j][dim]) / (mass[i]+mass[j]);
          v[i][dim] = (mass[i]*v[i][dim] + mass[j]*v[j][dim]) / (mass[i]+mass[j]);
        }

        // Mass of merged object
        mass[i] += mass[j];

        // Decrement n bodies
        const int l = --NumberOfBodies;

        // Remove other merged object from list
        #pragma omp simd private(dim)
        for (dim = 0; dim < 3; dim++) {
	      x[j][dim] = x[l][dim];
	      v[j][dim] = v[l][dim];
        }
        mass[j] = mass[l];
        j--;
        return;
      }

      #pragma omp simd private(dim)
      for (dim = 0; dim < 3; dim++) aTemp[0][dim] = (x[j][dim]-x[i][dim]) * mass[j] / (dist*dist*dist);
      minDx = std::min(minDx,dist);

        // Step 2-------------------------------------------------------------
      #pragma omp simd private(dim)
      for (dim = 0; dim < 3; dim++) {
        vTemp[1][dim] = v[i][dim] + aTemp[0][dim]*timeStepSize*0.5;  // compute 2nd order v
        xTemp[1][dim] = x[i][dim] + v[i][dim]*timeStepSize*0.5; // compute 2nd order x
        d[dim] = x[j][dim] - xTemp[1][dim];  // compute dx,dy,dz of 2nd order x
      }
      //#pragma omp barrier
      dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 2nd order x to j
      #pragma omp simd private(dim)
      for (dim = 0; dim < 3; dim++) aTemp[1][dim] = d[dim]*mass[j]/(dist*dist*dist);  // 2nd order acceleration of i

      // Step 3-------------------------------------------------------------
      #pragma omp simd private(dim)
      for (dim = 0; dim < 3; dim++) {
        vTemp[2][dim] = v[i][dim] + aTemp[1][dim]*timeStepSize*0.5;  // compute 3rd order v
        xTemp[2][dim] = x[i][dim] + vTemp[1][dim]*timeStepSize*0.5; // compute 3rd order x
        d[dim] = x[j][dim] - xTemp[2][dim];  // compute dx,dy,dz of 3rd order x
      }
      //#pragma omp barrier
      dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 3rd order x to j
      #pragma omp simd private(dim)
      for (dim = 0; dim < 3; dim++) aTemp[2][dim] = d[dim]*mass[j]/(dist*dist*dist);  // 3rd order acceleration of i

      // Step 4-------------------------------------------------------------
      #pragma omp simd private(dim)
      for (dim = 0; dim < 3; dim++) {
        vTemp[3][dim] = v[i][dim] + aTemp[2][dim]*timeStepSize*0.5;  // compute 4th order v
        xTemp[3][dim] = x[i][dim] + vTemp[2][dim]*timeStepSize*0.5; // compute 4th order x
        d[dim] = x[j][dim] -  xTemp[3][dim];  // compute dx,dy,dz of 4th order x
      }
      //#pragma omp barrier
      dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 4th order x to j
      #pragma omp simd private(dim)
      for (dim = 0; dim < 3; dim++) aTemp[3][dim] = d[dim]*mass[j]/(dist*dist*dist);  // 4th order acceleration of i


      // Update x and v-----------------------------------------------------
      //#pragma omp barrier
      #pragma omp simd private(dim)
      for (dim = 0; dim < 3; dim++) {
        x[i][dim] = x[i][dim] + nr*(v[i][dim] + 2*vTemp[1][dim] + 2*vTemp[2][dim] + vTemp[3][dim]) * timeStepSize;
        v[i][dim] = v[i][dim] + nr*(aTemp[0][dim] + 2*aTemp[1][dim] + 2*aTemp[2][dim] + aTemp[3][dim]) * timeStepSize;
      }
    }

    void updateBody () {
      timeStepCounter++;
      maxV   = 0.0;
      minDx  = std::numeric_limits<double>::max();
      for (int i=0; i<NumberOfBodies; i++) {
            // Possible vectorisation and parallism on inner loop
        for (int j = 0; j < NumberOfBodies; j++) {
          // Calculate position and velocity by performing RK4 on i and j
          if (i==j) continue;
          rk4(i,j);
          maxV = std::max(std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]), maxV);
        }
      }
      t += timeStepSize;
    }
};
