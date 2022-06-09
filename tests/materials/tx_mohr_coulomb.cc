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

//! Simulate Mohr-Coulomb TXC in 3D
TEST_CASE("Mohr-Coulomb TXC", "[TXC][MohrCoulomb]") {

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim>>(pid, coords);

  // Initialise material
  Json jmaterial;

  // Base parameters
  jmaterial["density"] = 2000.0;
  jmaterial["youngs_modulus"] = 10.0E+6;
  jmaterial["poisson_ratio"] = 0.2;
  jmaterial["friction"] = 33.0;
  jmaterial["dilation"] = 0.0;
  jmaterial["cohesion"] = 0.0;
  jmaterial["softening"] = false;
  jmaterial["residual_friction"] = 0.0;
  jmaterial["residual_dilation"] = 0.0;
  jmaterial["residual_cohesion"] = 0.0;
  jmaterial["peak_pdstrain"] = 0.0;
  jmaterial["residual_pdstrain"] = 0.0;
  jmaterial["tension_cutoff"] = 10.0;

  // Check Mohr-Coulomb TXC
  SECTION("Mohr-Coulomb TXC") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -20000.0;
    stress(1) = -20000.0;
    stress(2) = -20001.0;

    // Check dmatrix
    const double poisson_ratio = jmaterial["poisson_ratio"];
    const double youngs_modulus = jmaterial["youngs_modulus"];
    const double K = youngs_modulus / (3. * (1. - (2. * poisson_ratio)));
    const double G = youngs_modulus / (2. * (1. + poisson_ratio));
    const double a1 = K + (4. / 3.) * G;
    const double a2 = K - (2. / 3.) * G;

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.00001;
    dstrain(0) = -1. * dstrain(2) * a2 / (a2 + a1);
    dstrain(1) = -1. * dstrain(2) * a2 / (a2 + a1);

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("mohr-coulomb.txt");
    double axial_strain = 0.;
    double vol_strain = 0.;

    // Write initial states
    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
           << axial_strain << '\t' << vol_strain << '\n';

    // Loop
    for (unsigned i = 0; i < 20001; ++i) {
      axial_strain += dstrain(2);
      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      if ((state_vars).at("yield_state") == 1) {
        dstrain(2) = -0.00001;
        dstrain(0) = -0.5 * dstrain(2);
        dstrain(1) = -0.5 * dstrain(2);
      }

      myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
             << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
             << axial_strain << '\t' << vol_strain << '\n';
    }
    myfile.close();
  }

  // Loose parameters
  jmaterial["friction"] = 35.0;
  jmaterial["dilation"] = 3.0;
  jmaterial["softening"] = true;
  jmaterial["residual_friction"] = 33.0;
  jmaterial["residual_dilation"] = 0.0;
  jmaterial["peak_pdstrain"] = 0.05;
  jmaterial["residual_pdstrain"] = 0.28;

  // Check Mohr-Coulomb Strain Softening TXC in loose condition
  SECTION("Mohr-Coulomb Strain Softening TXC (LOOSE)") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -20000.0;
    stress(1) = -20000.0;
    stress(2) = -20001.0;

    // Check dmatrix
    const double poisson_ratio = jmaterial["poisson_ratio"];
    const double youngs_modulus = jmaterial["youngs_modulus"];
    const double K = youngs_modulus / (3. * (1. - (2. * poisson_ratio)));
    const double G = youngs_modulus / (2. * (1. + poisson_ratio));
    const double a1 = K + (4. / 3.) * G;
    const double a2 = K - (2. / 3.) * G;

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.00001;
    dstrain(0) = -1. * dstrain(2) * a2 / (a2 + a1);
    dstrain(1) = -1. * dstrain(2) * a2 / (a2 + a1);

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("mohr-coulomb-loose.txt");
    double axial_strain = 0.;
    double vol_dstrain = 0.;

    // Write initial states
    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
           << axial_strain << '\t' << vol_dstrain << '\n';

    // Loop
    for (unsigned i = 0; i < 20001; ++i) {
      axial_strain += dstrain(2);

      Eigen::Matrix<double, 6, 6> Dep =
          material->compute_consistent_tangent_matrix(
              stress, stress, dstrain, particle.get(), &state_vars);

      dstrain(0) = -1. * Dep(0, 2) / (Dep(0, 0) + Dep(0, 1)) * dstrain(2);
      dstrain(1) = -1. * Dep(1, 2) / (Dep(1, 0) + Dep(1, 1)) * dstrain(2);

      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      vol_dstrain = dstrain(0) + dstrain(1) + dstrain(2);

      myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
             << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
             << axial_strain << '\t' << vol_dstrain << '\n';
    }
    myfile.close();
  }

  // Dense parameters
  jmaterial["friction"] = 43.0;
  jmaterial["dilation"] = 17.0;
  jmaterial["softening"] = true;
  jmaterial["residual_friction"] = 33.0;
  jmaterial["residual_dilation"] = 0.0;
  jmaterial["peak_pdstrain"] = 0.03;
  jmaterial["residual_pdstrain"] = 0.36;

  // Check Mohr-Coulomb Strain Softening TXC in dense condition
  SECTION("Mohr-Coulomb Strain Softening TXC (DENSE)") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -20000.0;
    stress(1) = -20000.0;
    stress(2) = -20001.0;

    // Check dmatrix
    const double poisson_ratio = jmaterial["poisson_ratio"];
    const double youngs_modulus = jmaterial["youngs_modulus"];
    const double K = youngs_modulus / (3. * (1. - (2. * poisson_ratio)));
    const double G = youngs_modulus / (2. * (1. + poisson_ratio));
    const double a1 = K + (4. / 3.) * G;
    const double a2 = K - (2. / 3.) * G;

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.00001;
    dstrain(0) = -1. * dstrain(2) * a2 / (a2 + a1);
    dstrain(1) = -1. * dstrain(2) * a2 / (a2 + a1);

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("mohr-coulomb-dense.txt");
    double axial_strain = 0.;
    double vol_dstrain = 0.;

    // Write initial states
    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
           << axial_strain << '\t' << vol_dstrain << '\n';

    // Loop
    for (unsigned i = 0; i < 20001; ++i) {
      axial_strain += dstrain(2);

      Eigen::Matrix<double, 6, 6> Dep =
          material->compute_consistent_tangent_matrix(
              stress, stress, dstrain, particle.get(), &state_vars);

      dstrain(0) = -1. * Dep(0, 2) / (Dep(0, 0) + Dep(0, 1)) * dstrain(2);
      dstrain(1) = -1. * Dep(1, 2) / (Dep(1, 0) + Dep(1, 1)) * dstrain(2);

      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      vol_dstrain = dstrain(0) + dstrain(1) + dstrain(2);

      myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
             << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
             << axial_strain << '\t' << vol_dstrain << '\n';
    }
    myfile.close();
  }
}