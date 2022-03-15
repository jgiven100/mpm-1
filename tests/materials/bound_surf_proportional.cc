#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include <cmath>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "materials/material.h"
#include "node.h"
#include "particle.h"

//! Check Bounding Surface Plasticity class in 3D
TEST_CASE("Bound Surface Plasticity (proportional) is checked in 3D", "[material][BoundSurfPlasticity][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim>>(pid, coords);

  // Initialise material
  Json jmaterial;

  // Testing
  jmaterial["density"] = 2000.;
  jmaterial["friction"] = 28.;
  jmaterial["G0"] = 100.;
  jmaterial["poisson_ratio"] = 0.33;
  jmaterial["critical_void_ratio"] = 0.72;
  jmaterial["relative_density"] = 0.43;
  jmaterial["hr0"] = 0.25;
  jmaterial["kr0"] = 0.2;
  jmaterial["d0"] = 2.5;
  jmaterial["fp"] = 0.5;
  jmaterial["proportional"] = true;
  jmaterial["cz"] = 20;
  jmaterial["zm"] = 5.5;

  // Check proportional loading
  SECTION("BoundSurfPlasticity stress control anisotropic initial stress") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "BoundSurfPlasticity3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -41000;
    stress(1) = -82001;
    stress(2) = -61500;

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("proportional.txt");

    // Initialise shear strain
    double shear_strain = 0.;

    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
           << shear_strain << '\n';

    double sign = 1.0;
    int count = 0;

    // Loop
    for (unsigned i = 0; i < 60000 + 1; ++i) {

      // Check stress ratio
      if (std::fabs(stress(3)/82001) > 0.2 ) {
        if (count > 10) {
          sign *= -1.;
          count = 0;
        }
      }

      // Set dstrain gamma_{23}
      dstrain(3) = sign * -0.00001;
      shear_strain += dstrain(3);

      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
             << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
             << shear_strain << '\n';

      count += 1;
    }
    myfile.close();
  }
} 
