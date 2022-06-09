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

// MUST ADD TO CMAKE TO GENERATE OUTPUTS!!

//! Simulate NorSand TXC in 3D
TEST_CASE("NorSand TXC Calibration", "[TXC][NorSand][Calibration]") {

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

  // jmaterial["use_bolton_csl"] = true;
  // jmaterial["e_min"] = 0.500;
  // jmaterial["e_max"] = 0.946;
  // jmaterial["crushing_pressure"] = 10E+6;
  jmaterial["gamma"] = 0.91;
  jmaterial["reference_pressure"] = 100.0E+3;

  jmaterial["N"] = 0.16;
  jmaterial["chi"] = 5.3;
  jmaterial["hardening_modulus"] = 150.0;
  jmaterial["lambda"] = 0.03;
  jmaterial["kappa"] = 0.002;

  // jmaterial["bond_model"] = false;
  // jmaterial["p_cohesion_initial"] = 0.0E+0;
  // jmaterial["p_dilation_initial"] = 0.0E+0;
  // jmaterial["m_cohesion"] = 0.0;
  // jmaterial["m_dilation"] = 0.0;
  // jmaterial["m_modulus"] = 0.0;

  // Loose parameter
  jmaterial["void_ratio_initial"] = 0.839;
  jmaterial["p_image_initial"] = 7360;
  // Check NorSand TXC in loose condition
  SECTION("NorSand, loose, 20 kPa") {
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
    myfile.open("tokyotest10-1.txt");
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

  // Check NorSand TXC in loose condition
  jmaterial["void_ratio_initial"] = 0.836;
  jmaterial["p_image_initial"] = 14720;
  SECTION("NorSand, loose, 40 kPa") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -40000.0;
    stress(1) = -40000.0;
    stress(2) = -40001.0;

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
    myfile.open("tokyotest10-2.txt");
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

  // Check NorSand TXC in loose condition
  jmaterial["void_ratio_initial"] = 0.828;
  jmaterial["p_image_initial"] = 29440;
  SECTION("NorSand, loose, 80 kPa") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -80000.0;
    stress(1) = -80000.0;
    stress(2) = -80001.0;

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
    myfile.open("tokyotest10-3.txt");
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

  // Medium parameter
  jmaterial["void_ratio_initial"] = 0.739;
  jmaterial["p_image_initial"] = 7360;
  // Check NorSand TXC in medium condition
  SECTION("NorSand, medium, 20 kPa") {
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
    myfile.open("tokyotest11-1.txt");
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

  // Check NorSand TXC in medium condition
  jmaterial["void_ratio_initial"] = 0.735;
  jmaterial["p_image_initial"] = 14720;
  SECTION("NorSand, medium, 40 kPa") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -40000.0;
    stress(1) = -40000.0;
    stress(2) = -40001.0;

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
    myfile.open("tokyotest11-2.txt");
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

  // Check NorSand TXC in medium condition
  jmaterial["void_ratio_initial"] = 0.724;
  jmaterial["p_image_initial"] = 29440;
  SECTION("NorSand, medium, 80 kPa") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -80000.0;
    stress(1) = -80000.0;
    stress(2) = -80001.0;

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
    myfile.open("tokyotest11-3.txt");
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
  jmaterial["void_ratio_initial"] = 0.663;
  jmaterial["p_image_initial"] = 7360;
  // Check NorSand TXC in dense condition
  SECTION("NorSand, dense, 20 kPa") {
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
    myfile.open("tokyotest12-1.txt");
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

  // Check NorSand TXC in dense condition
  jmaterial["void_ratio_initial"] = 0.647;
  jmaterial["p_image_initial"] = 14720;
  SECTION("NorSand, dense, 40 kPa") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -40000.0;
    stress(1) = -40000.0;
    stress(2) = -40001.0;

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
    myfile.open("tokyotest12-2.txt");
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

  // Check NorSand TXC in dense condition
  jmaterial["void_ratio_initial"] = 0.641;
  jmaterial["p_image_initial"] = 29440;
  SECTION("NorSand, dense, 80 kPa") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -80000.0;
    stress(1) = -80000.0;
    stress(2) = -80001.0;

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
    myfile.open("tokyotest12-3.txt");
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