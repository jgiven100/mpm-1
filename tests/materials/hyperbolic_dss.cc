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

//! Check Hyperbolic
TEST_CASE("Hyperbolic is checked in 3D", "[material][Hyperbolic][3D]") {

  // Dimension
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
  jmaterial["poisson_ratio"] = 0.25;
  jmaterial["youngs_modulus"] = 10.E+6;
  jmaterial["beta"] = 1.;
  jmaterial["s"] = 1.;
  jmaterial["reference_strain"] = 0.01;
  jmaterial["p1"] = 1.;
  jmaterial["p2"] = 1.;
  jmaterial["p3"] = 1.;

  // Check DSS loading
  SECTION("Hyperbolic DSS") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Hyperbolic3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -1.0E+6;
    stress(1) = -1.0E+6;
    stress(2) = -1.0E+6 - 1;

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(4) = -0.00001;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("hyperbolic-dss.txt");
    double eng_shear_strain = 0.;

    // Write initial states
    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
           << eng_shear_strain << '\n';

    // Loop
    for (unsigned i = 0; i < 4500; ++i) {
      if (std::fabs(eng_shear_strain) > 0.00499) dstrain(4) *= -1.;
      eng_shear_strain += dstrain(4);
      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
             << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
             << eng_shear_strain << '\n';
    }
    myfile.close();
  }
}