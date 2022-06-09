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

//! Simulate NorSand TXC in 3D
TEST_CASE("NorSand TXC", "[TXC][NorSand]") {

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim>>(pid, coords);

  // Initialise material
  Json jmaterial;

  // Paramters
  jmaterial["density"] = 2000.0;
  jmaterial["poisson_ratio"] = 0.2;
  jmaterial["friction_cs"] = 33.0;
  jmaterial["p_image_initial"] = 7360;

  // jmaterial["use_bolton_csl"] = true;
  // jmaterial["e_min"] = 0.500;
  // jmaterial["e_max"] = 0.946;
  // jmaterial["crushing_pressure"] = 10E+6;
  jmaterial["gamma"] = 0.91;
  jmaterial["reference_pressure"] = 100.0E+3;

  jmaterial["N"] = 0.2;
  jmaterial["chi"] = 3.0;
  jmaterial["hardening_modulus"] = 450.0;
  jmaterial["lambda"] = 0.03;
  jmaterial["kappa"] = 0.002;

  // jmaterial["bond_model"] = false;
  // jmaterial["p_cohesion_initial"] = 0.0E+0;
  // jmaterial["p_dilation_initial"] = 0.0E+0;
  // jmaterial["m_cohesion"] = 0.0;
  // jmaterial["m_dilation"] = 0.0;
  // jmaterial["m_modulus"] = 0.0;

  // Loose parameter
  jmaterial["void_ratio_initial"] = 0.9;

  // Check NorSand TXC in loose condition
  SECTION("NorSand (loose)") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -20000.0;
    stress(1) = -20000.0;
    stress(2) = -20001.0;

    // Check dmatrix
    const double nu = jmaterial["poisson_ratio"];
    const double kappa = jmaterial["kappa"];
    const double e_init = jmaterial["void_ratio_initial"];
    const double p = (stress(0) + stress(1) + stress(2)) / 3.0;
    const double K = (1. + e_init) / (kappa * p);
    const double G = 3. * K * (1. - 2. * nu) / (2. * (1. + nu));
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
    myfile.open("norsand-loose.txt");
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

  // Dense parameter
  jmaterial["void_ratio_initial"] = 0.7;

  // Check NorSand TXC in dense condition
  SECTION("NorSand (dense)") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -20000.0;
    stress(1) = -20000.0;
    stress(2) = -20001.0;

    // Check dmatrix
    const double nu = jmaterial["poisson_ratio"];
    const double kappa = jmaterial["kappa"];
    const double e_init = jmaterial["void_ratio_initial"];
    const double p = (stress(0) + stress(1) + stress(2)) / 3.0;
    const double K = (1. + e_init) / (kappa * p);
    const double G = 3. * K * (1. - 2. * nu) / (2. * (1. + nu));
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
    myfile.open("norsand-dense.txt");
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