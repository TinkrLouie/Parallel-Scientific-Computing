#include "NBodySimulation.h"
# include <omp.h>

// TODO: Try using structs istead

class NBodySimulationVectorised : public NBodySimulation {
  protected:
  public:
    #pragma omp declare simd
    void collision(int i, int j) {
      // Momentum calculations
      #pragma omp simd
      for (int dim = 0; dim < 3; dim++) {
        x[i][dim] = (mass[i]*x[i][dim] + mass[j]*x[j][dim]) / (mass[i]+mass[j]);
        v[i][dim] = (mass[i]*v[i][dim] + mass[j]*v[j][dim]) / (mass[i]+mass[j]);
      }
      
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
      #pragma omp simd
      for (int dim = 0; dim < 3; dim++) {
	    x[j][dim] = x[l][dim];
	    v[j][dim] = v[l][dim];
      }
      mass[j] = mass[l];
      j--;
    }
  
    // Runge-Kutta 4th Order
    #pragma omp declare simd
    bool rk4 (int i, int j) {
        double* x2 = new double[3];  // 2nd order x
        double* x3 = new double[3];  // ...
        double* x4 = new double[3];  // 4th order x
        double* v2 = new double[3];  // 2nd order v
        double* v3 = new double[3];  // ...
        double* v4 = new double[3];  // 4th order v
        double* a1 = new double[3];  // 1st order acceleration
        double* a2 = new double[3];  // ...
        double* a3 = new double[3];  // ...
        double* a4 = new double[3];  // 4th order acceleration
        double* d = new double[3];   // dx, dy, dz
        double dist, nr = 1.0/6;
        double tolerance = 0.06;
        
        // Step 1-------------------------------------------------------------
        dist = distance(i,j);
        // Collision detection
        if (dist <= (0.01/NumberOfBodies)*(mass[i] + mass[j]) + tolerance){
          collision(i,j);
          return false;
        }

        #pragma omp simd
        for (int dim = 0; dim < 3; dim++) a1[dim] = acceleration(j,i,dist,dim);
  
        // Step 2-------------------------------------------------------------
        #pragma omp simd
        for (int dim = 0; dim < 3; dim++) {
          v2[dim] = v[i][dim] + a1[dim]*timeStepSize*0.5;  // compute 2nd order v
          x2[dim] = x[i][dim] + v[i][dim]*timeStepSize*0.5; // compute 2nd order x
          d[dim] = x[j][dim] - x2[dim];  // compute dx,dy,dz of 2nd order x
        }
        dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 2nd order x to j
        #pragma omp simd
        for (int dim = 0; dim < 3; dim++) a2[dim] = d[dim]*mass[j]/(dist*dist*dist);  // 2nd order acceleration of i
  
        // Step 3-------------------------------------------------------------
        #pragma omp simd
        for (int dim = 0; dim < 3; dim++) {
          v3[dim] = v[i][dim] + a2[dim]*timeStepSize*0.5;  // compute 3rd order v
          x3[dim] = x[i][dim] + v2[dim]*timeStepSize*0.5; // compute 3rd order x
          d[dim] = x[j][dim] - x3[dim];  // compute dx,dy,dz of 3rd order x
        }
        dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 3rd order x to j
        #pragma omp simd
        for (int dim = 0; dim < 3; dim++) a3[dim] = d[dim]*mass[j]/(dist*dist*dist);  // 3rd order acceleration of i
  
        // Step 4-------------------------------------------------------------
        for (int dim = 0; dim < 3; dim++) {
          v4[dim] = v[i][dim] + a3[dim]*timeStepSize*0.5;  // compute 4th order v
          x4[dim] = x[i][dim] + v3[dim]*timeStepSize*0.5; // compute 4th order x
          d[dim] = x[j][dim] - x4[dim];  // compute dx,dy,dz of 4th order x
        }
        dist = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);  // distance between 4th order x to j
        #pragma omp simd
        for (int dim = 0; dim < 3; dim++) a4[dim] = d[dim]*mass[j]/(dist*dist*dist);  // 4th order acceleration of i
  
  
        // Update x and v-----------------------------------------------------
        #pragma omp simd
        for (int dim = 0; dim < 3; dim++) {
          x[i][dim] = x[i][dim] + nr*(v[i][dim] + 2*v2[dim] + 2*v3[dim] + v4[dim]) * timeStepSize;
          v[i][dim] = v[i][dim] + nr*(a1[dim] + 2*a2[dim] + 2*a3[dim] + a4[dim]) * timeStepSize;
        }
        
        // Deallocate memory
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
        return true;
    }

    #pragma omp declare simd
    double distance (int i, int j) {  // Euclidean distance between bodies i and j
      const double distance = sqrt(
        (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
        (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
        (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
        );
      return distance;
    }

    #pragma omp declare simd
    double acceleration (int i, int j, double distance, int direction){ // acceleration of body i through body j
      const double distance3 = distance * distance * distance;
      minDx = std::min( minDx,distance );  // can be taken out
  
      return (x[i][direction]-x[j][direction]) * mass[i] / distance3;
    }
  
    void updateBody () {
      timeStepCounter++;
      maxV   = 0.0;
      minDx  = std::numeric_limits<double>::max();
      //double c = 0.01/NumberOfBodies;
      #pragma omp parallel for
      for (int i=0; i<NumberOfBodies; i++) {
        #pragma omp simd reduction(max:maxV) reduction(min:minDx)
        for (int j = i+1; j < NumberOfBodies; j++) {
          // Calculate position and velocity by performing RK4 on i and j, if no collision, do the reverse
          if (rk4(i,j)) rk4(j,i);
          maxV = std::max(std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]), maxV);
        }
      }
      t += timeStepSize;
    }
};
