#include "glCanvas.h"

#include <fstream>
#include "cloth.h"
#include "argparser.h"
#include "utils.h"
#include <cmath>

// ================================================================================
// ================================================================================

Cloth::Cloth(ArgParser *_args) {
  args =_args;

  // open the file
  std::ifstream istr(std::string(args->path+'/'+args->cloth_file).c_str());
  assert (istr.good());
  std::string token;

  // read in the simulation parameters
  istr >> token >> k_structural; assert (token == "k_structural");  // (units == N/m)  (N = kg*m/s^2)
  istr >> token >> k_shear; assert (token == "k_shear");
  istr >> token >> k_bend; assert (token == "k_bend");
  istr >> token >> damping; assert (token == "damping");
  // NOTE: correction factor == .1, means springs shouldn't stretch more than 10%
  //       correction factor == 100, means don't do any correction
  istr >> token >> provot_structural_correction; assert (token == "provot_structural_correction");
  istr >> token >> provot_shear_correction; assert (token == "provot_shear_correction");

  // the cloth dimensions
  istr >> token >> nx >> ny;
  assert (token == "m");
  assert (nx >= 2 && ny >= 2);

  // the corners of the cloth
  // (units == meters)
  glm::vec3 a,b,c,d;
  istr >> token >> a.x >> a.y >> a.z; assert (token == "p");
  istr >> token >> b.x >> b.y >> b.z; assert (token == "p");
  istr >> token >> c.x >> c.y >> c.z; assert (token == "p");
  istr >> token >> d.x >> d.y >> d.z; assert (token == "p");

  // fabric weight  (units == kg/m^2)
  // denim ~300 g/m^2
  // silk ~70 g/m^2
  double fabric_weight;
  istr >> token >> fabric_weight; assert (token == "fabric_weight");
  double area = AreaOfTriangle(a,b,c) + AreaOfTriangle(a,c,d);

  // create the particles
  particles = new ClothParticle[nx*ny];
  double mass = area*fabric_weight / double(nx*ny);
  for (int i = 0; i < nx; i++) {
    double x = i/double(nx-1);
    glm::vec3 ab = float(1-x)*a + float(x)*b;
    glm::vec3 dc = float(1-x)*d + float(x)*c;
    for (int j = 0; j < ny; j++) {
      double y = j/double(ny-1);
      ClothParticle &p = getParticle(i,j);
      glm::vec3 abdc = float(1-y)*ab + float(y)*dc;
      p.setOriginalPosition(abdc);
      p.setPosition(abdc);
      p.setVelocity(glm::vec3(0,0,0));
      p.setMass(mass);
      p.setFixed(false);
    }
  }

  // the fixed particles
  while (istr >> token) {
    assert (token == "f");
    int i,j;
    double x,y,z;
    istr >> i >> j >> x >> y >> z;
    ClothParticle &p = getParticle(i,j);
    p.setPosition(glm::vec3(x,y,z));
    p.setFixed(true);
  }

  computeBoundingBox();

  if (args->method == "midpoint") {
    printf("Using MidPoint Integration\n");
    method = IntegrationMethod::MidPoint;
  } else if (args->method == "rungekutta") {
    printf("Using RungeKutta Integration\n");
    method = IntegrationMethod::RungeKutta;
  }
  else {
    printf("Using Euler Integration\n");
    method = IntegrationMethod::Euler;
  }
}

// ================================================================================

void Cloth::computeBoundingBox() {
  box = BoundingBox(getParticle(0,0).getPosition());
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      box.Extend(getParticle(i,j).getPosition());
      box.Extend(getParticle(i,j).getOriginalPosition());
    }
  }
}

// ================================================================================



void Cloth::Animate() {


  // *********************************************************************
  // ASSIGNMENT:
  //
  // Compute the forces on each particle, and update the state
  // (position & velocity) of each particle.
  //
  // Also, this is where you'll put the Provot correction for super-elasticity
  //
  // *********************************************************************
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      ClothParticle &p = particles[i + j*nx];
      updateState(p, i, j);
    }
  }

  for (int k=0; k<10; k++) {
    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        ClothParticle &p = particles[i + j*nx];
        provotCorrection(p, i, j);
      }
    }
  }
  // commented out because an animated bounding box can be weird
  //computeBoundingBox();
}

// ================================================================================
// Assignment functions
// ================================================================================

void Cloth::provotCorrection(ClothParticle& p, int i, int j) {
  //printf("Basic: %f, Shear: %f\n", provot_structural_correction, provot_shear_correction);
  if (isParticle(i-1, j)) {correct(p, getParticle(i-1, j), provot_structural_correction);}
  if (isParticle(i, j-1)) {correct(p, getParticle(i, j-1), provot_structural_correction);}
  if (isParticle(i+1, j)) {correct(p, getParticle(i+1, j), provot_structural_correction);}
  if (isParticle(i, j+1)) {correct(p, getParticle(i, j+1), provot_structural_correction);}

  if (isParticle(i-1, j-1)) {correct(p, getParticle(i-1, j-1), provot_shear_correction);}
  if (isParticle(i-1, j+1)) {correct(p, getParticle(i-1, j+1), provot_shear_correction);}
  if (isParticle(i+1, j-1)) {correct(p, getParticle(i+1, j-1), provot_shear_correction);}
  if (isParticle(i+1, j+1)) {correct(p, getParticle(i+1, j+1), provot_shear_correction);}


}

void Cloth::correct(ClothParticle& p1, ClothParticle& p2, double p_correct) {
  double restLength = glm::length(p2.getOriginalPosition() - p1.getOriginalPosition());
  double dist = glm::length(p2.getPosition() - p1.getPosition());
  double stretch = dist/restLength;
  if (stretch > (1.0+p_correct) && p2.isFixed() == false) {
    double correction = (stretch - (1.0+p_correct))/2.0; // The amount each particle needs to move

    glm::vec3 p1_pos = p1.getPosition();
    glm::vec3 p2_pos = p2.getPosition();

    if (p1.isFixed() == false) {
      p1.setPosition(p1_pos + normalize(p2_pos - p1_pos)*(float)restLength*(float)(correction));
      p2.setPosition(p2_pos + normalize(p1_pos - p2_pos)*(float)restLength*(float)(correction));
    }
    else {
      p2.setPosition(p2_pos + normalize(p1_pos - p2_pos)*(float)restLength*(float)(correction));
    }
  }

}

std::string Cloth::stringifyVector(const glm::vec3& vec) {
    return "(" + std::to_string(vec.x) + ","+ std::to_string(vec.y) + "," + std::to_string(vec.z) +")";
}

void Cloth::updateState(ClothParticle& p, int i_index, int j_index) {

  switch (method) {
      case (IntegrationMethod::Euler):
          integrationEuler(p);
          break;
      case (IntegrationMethod::MidPoint):
          integrationMidpoint(p, i_index, j_index);
          break;
      case (IntegrationMethod::RungeKutta):
          integrationRungeKutta(p, i_index, j_index);
          break;
      default:
          printf("Cannot Integrate, no known method\n");
          break;
  }

  updateForces(p, i_index, j_index);
}

// Euler takes the current pos, moves it, then adjusts velocity
void Cloth::integrationEuler(ClothParticle& p) {
  if (!p.isFixed()) {
    p.setPosition(p.getPosition()+(p.getVelocity()*(float)args->timestep));
    p.setVelocity(p.getVelocity()+(p.getAcceleration()*(float)args->timestep));
  } else {
    p.setVelocity(glm::vec3(0.0f, 0.0f, 0.0f));
  }
}

void Cloth::integrationMidpoint(ClothParticle& p, int i, int j) {
  if (!p.isFixed()) {
    glm::vec3 mid = p.getPosition()+(p.getVelocity()*(float)args->timestep*0.5f);
    glm::vec3 midAcc = getForces(p, mid, i, j);
    glm::vec3 midVel = p.getVelocity()+(midAcc*(float)args->timestep*0.5f);
    p.setPosition(mid+(midVel*(float)args->timestep*0.5f));
    p.setVelocity(midVel);
  } else {
    p.setVelocity(glm::vec3(0.0f, 0.0f, 0.0f));
  }
}

void Cloth::integrationRungeKutta(ClothParticle& p, int i, int j) {
  if (!p.isFixed()) {
    glm::vec3 mid = p.getPosition()+(p.getVelocity()*(float)args->timestep*0.5f);
    glm::vec3 midAcc = getForces(p, mid, i, j);
    glm::vec3 midVel = p.getVelocity()+(midAcc*(float)args->timestep*0.5f);
    glm::vec3 end = p.getPosition()+(p.getVelocity()*(float)args->timestep);
    glm::vec3 endAcc = getForces(p, end, i, j);
    glm::vec3 endVel = midVel+(endAcc*(float)args->timestep*0.5f);

    glm::vec3 k1 = p.getVelocity()*(float)args->timestep;
    glm::vec3 k2 = (midVel+(0.5f*k1))*(float)args->timestep;
    glm::vec3 k3 = (midVel+(0.5f*k2))*(float)args->timestep;
    glm::vec3 k4 = (endVel+k3)*(float)args->timestep;

    p.setPosition(p.getPosition()+((k1 + 2.0f*k2 + 2.0f*k3 + k4)/6.0f));
    p.setVelocity(p.getVelocity()+(p.getAcceleration()*(float)args->timestep));


  } else {
    p.setVelocity(glm::vec3(0.0f, 0.0f, 0.0f));
  }
}

void Cloth::updateForces(ClothParticle& p, int i_index, int j_index) {
  glm::vec3 totalForce = glm::vec3(0.0f, 0.0f, 0.0f);

  totalForce += args->gravity;
  totalForce += (calcSpring(Spring::Structural, p, i_index, j_index)/(float)p.getMass());
  totalForce += (calcSpring(Spring::Shear, p, i_index, j_index)/(float)p.getMass());
  totalForce += (calcSpring(Spring::Flexion, p, i_index, j_index)/(float)p.getMass());

  totalForce += calcDamping(p.getVelocity())/(float)p.getMass();

  p.setAcceleration(totalForce);
}

// HAHA Im passing by copy and no one can stop me!
glm::vec3 Cloth::getForces(ClothParticle p, glm::vec3 mid, int i_index, int j_index) {
  glm::vec3 totalForce = glm::vec3(0.0f, 0.0f, 0.0f);
  p.setPosition(mid);

  totalForce += args->gravity;
  totalForce += (calcSpring(Spring::Structural, p, i_index, j_index)/(float)p.getMass());
  totalForce += (calcSpring(Spring::Shear, p, i_index, j_index)/(float)p.getMass());
  totalForce += (calcSpring(Spring::Flexion, p, i_index, j_index)/(float)p.getMass());

  totalForce += calcDamping(p.getVelocity())/(float)p.getMass();

  return totalForce;
}

glm::vec3 Cloth::calcSpring(Spring s, ClothParticle& pi, int i, int j) {
  glm::vec3 totalForce = glm::vec3(0.0f,0.0f,0.0f);
  switch (s) {
    case (Spring::Structural):
      if (isParticle(i-1, j)) {totalForce += calcSpringForce(k_structural, pi, getParticle(i-1, j));}
      if (isParticle(i, j-1)) {totalForce += calcSpringForce(k_structural, pi, getParticle(i, j-1));}
      if (isParticle(i+1, j)) {totalForce += calcSpringForce(k_structural, pi, getParticle(i+1, j));}
      if (isParticle(i, j+1)) {totalForce += calcSpringForce(k_structural, pi, getParticle(i, j+1));}
      break;
    case (Spring::Shear):
      if (isParticle(i-1, j-1)) {totalForce += calcSpringForce(k_shear, pi, getParticle(i-1, j-1));}
      if (isParticle(i-1, j+1)) {totalForce += calcSpringForce(k_shear, pi, getParticle(i-1, j+1));}
      if (isParticle(i+1, j-1)) {totalForce += calcSpringForce(k_shear, pi, getParticle(i+1, j-1));}
      if (isParticle(i+1, j+1)) {totalForce += calcSpringForce(k_shear, pi, getParticle(i+1, j+1));}
      break;
    case (Spring::Flexion):
      if (isParticle(i-2, j)) {totalForce += calcSpringForce(k_bend, pi, getParticle(i-2, j));}
      if (isParticle(i, j-2)) {totalForce += calcSpringForce(k_bend, pi, getParticle(i, j-2));}
      if (isParticle(i+2, j)) {totalForce += calcSpringForce(k_bend, pi, getParticle(i+2, j));}
      if (isParticle(i, j+2)) {totalForce += calcSpringForce(k_bend, pi, getParticle(i, j+2));}
      break;
    default:
      return glm::vec3(0.0f, 0.0f, 0.0f);
      break;
  }

  return totalForce;
}

glm::vec3 Cloth::calcSpringForce(double k, ClothParticle pi, ClothParticle pj) {
  if (pi.getPosition() == pj.getPosition()) {return glm::vec3(0.0f, 0.0f, 0.0f);}
  double restLength = std::abs(glm::length(pj.getOriginalPosition() - pi.getOriginalPosition()));
  double distance = glm::length(pj.getPosition() - pi.getPosition());
  glm::vec3 force = (float)(k * (restLength-distance)*-1.0) * ( (pj.getPosition() - pi.getPosition()) / (float)distance);
  //printf("restLength: %f, distance: %f, Spring Force: (%f,%f,%f)\n", restLength, distance, force.x, force.y, force.z);
  return force;
}

glm::vec3 Cloth::calcDamping(glm::vec3 vel) {
  return (-1.0f * vel * (float)damping);
}
