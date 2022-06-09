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

//! Simulate Bounding Surface Plasticity TXC in 3D
TEST_CASE("Bounding Surface Plasticity TXC", "[TXC][BoundSurfPlasticity]") {

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim>>(pid, coords);

  // Initialise material
  Json jmaterial;

  // Paramters
  jmaterial["density"] = 2000.;
  jmaterial["poisson_ratio"] = 0.2;
  jmaterial["friction"] = 34.0;

  jmaterial["G0"] = 450.0;
  jmaterial["hr0"] = 0.2;
  jmaterial["kr0"] = 0.5;
  jmaterial["d0"] = 10;
  jmaterial["gm0"] = 30.0;
  jmaterial["fp"] = 0.85;

  jmaterial["eta0"] = 5.0;

  // Optional [e_int = 0.9; lambda = 0.03]
  jmaterial["void_ratio_intercept"] = 0.95;
  jmaterial["lambda"] = 0.03;

  // Maximum and minimum void ratio
  jmaterial["maximum_void_ratio"] = 1.000;
  jmaterial["minimum_void_ratio"] = 0.500;

  // Loose parameters
  jmaterial["initial_void_ratio"] = 0.9;

  // Check Bounding Surface Plasticity TXC in loose condition
  SECTION("Bounding Surface Plasticity (loose)") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "BoundSurfPlasticity3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -20000.0;
    stress(1) = -20000.0;
    stress(2) = -20001.0;

    // Check dmatrix
    const double nu = jmaterial["poisson_ratio"];
    const double e_init = jmaterial["initial_void_ratio"];
    const double G0 = jmaterial["G0"];
    const double mean_p = -1. * (stress(0) + stress(1) + stress(2)) / 3.;
    const double p_atm = 100.e+3;
    const double G = p_atm * G0 *
                     (std::pow((2.973 - e_init), 2) / (1. + e_init)) *
                     std::sqrt(mean_p / p_atm);
    const double K = (2. * G * (1. + nu)) / (3. * (1. - 2. * nu));
    const double a1 = K + (4. / 3.) * G;
    const double a2 = K - (2. / 3.) * G;

    // Initialise dstrain (DRAINED)
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.00001;
    dstrain(0) = -1. * dstrain(2) * a2 / (a2 + a1);
    dstrain(1) = -1. * dstrain(2) * a2 / (a2 + a1);

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("bound-surf-loose.txt");
    double axial_strain = 0.;

    // Write initial states
    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
           << axial_strain << '\t' << (state_vars).at("void_ratio") << '\n';

    // Loop
    for (unsigned i = 0; i < 20001; ++i) {
      axial_strain += dstrain(2);
      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
             << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
             << axial_strain << '\t' << (state_vars).at("void_ratio") << '\n';
    }
    myfile.close();
  }

  // Dense parameters
  jmaterial["initial_void_ratio"] = 0.7;

  // Check Bounding Surface Plasticity TXC in dense condition
  SECTION("Bounding Surface Plasticity (dense)") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "BoundSurfPlasticity3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -20000.0;
    stress(1) = -20000.0;
    stress(2) = -20001.0;

    // // Check dmatrix
    const double nu = jmaterial["poisson_ratio"];
    const double e_init = jmaterial["initial_void_ratio"];
    const double G0 = jmaterial["G0"];
    const double mean_p = -1. * (stress(0) + stress(1) + stress(2)) / 3.;
    const double p_atm = 100.e+3;
    const double G = p_atm * G0 *
                     (std::pow((2.973 - e_init), 2) / (1. + e_init)) *
                     std::sqrt(mean_p / p_atm);
    const double K = (2. * G * (1. + nu)) / (3. * (1. - 2. * nu));
    const double a1 = K + (4. / 3.) * G;
    const double a2 = K - (2. / 3.) * G;

    // Initialise dstrain (DRAINED)
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.00001;
    dstrain(0) = -1. * dstrain(2) * a2 / (a2 + a1);
    dstrain(1) = -1. * dstrain(2) * a2 / (a2 + a1);

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("bound-surf-dense.txt");
    double axial_strain = 0.;

    // Write initial states
    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
           << axial_strain << '\t' << (state_vars).at("void_ratio") << '\n';

    // Loop
    for (unsigned i = 0; i < 20001; ++i) {
      axial_strain += dstrain(2);
      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
             << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
             << axial_strain << '\t' << (state_vars).at("void_ratio") << '\n';
    }
    myfile.close();
  }
}