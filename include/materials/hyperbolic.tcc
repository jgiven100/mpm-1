//! Read material properties
template <unsigned Tdim>
mpm::Hyperbolic<Tdim>::Hyperbolic(unsigned id, const Json& material_properties)
    : InfinitesimalElastoPlastic<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();

    // Material Properties for viscoelastic damping
    beta_ = material_properties.at("beta").template get<double>();
    s_ = material_properties.at("s").template get<double>();
    reference_strain_ =
        material_properties.at("reference_strain").template get<double>();
    p1_ = material_properties.at("p1").template get<double>();
    p2_ = material_properties.at("p2").template get<double>();
    p3_ = material_properties.at("p3").template get<double>();

  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::Hyperbolic<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {// Yield state: 0: elastic, 1: yield
                               {"yield_state", 0}};

  return state_vars;
}

//! Initialise state variables
template <unsigned Tdim>
std::vector<std::string> mpm::Hyperbolic<Tdim>::state_variables() const {
  const std::vector<std::string> state_vars = {"yield_state"};
  return state_vars;
}

//! Return elastic tensor
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6> mpm::Hyperbolic<Tdim>::compute_tensor(
    const Vector6d& strain) {
  // Shear modulus
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));
  // Calculate bulk modulus
  bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

  Vector6d G_mod_ = Vector6d::Zero();
  if (init_loading_) {
    G_mod_ =
        (G * (1 + beta_ * pow(strain.array().abs() / reference_strain_, s_))
                 .inverse())
            .matrix();
  } else {
    // Strain Vectors and subtracting reversal strain
    Vector6d strain_dif = (strain - reversal_strain_).cwiseAbs();

    // Calculating Reduction Factor
    double reduction_factor = (p1_ - p2_ * pow(1. - secant_modulus_ / G, p3_));

    // Calculating Second Term
    Vector6d term2 =
        (G *
         (1. + beta_ * pow(strain_dif.array() / (2 * reference_strain_), s_))
             .inverse())
            .matrix();

    // Calculating Third Term
    Vector6d term3 =
        (G *
         (1. +
          beta_ * pow((maximum_strain_.array() - reference_strain_).abs(), s_))
             .inverse())
            .matrix();

    G_mod_ = reduction_factor * (term2 - term3) + term3;
  }

  std::cout << G_mod_(0) << '\t' << G_mod_(1) << '\t' << G_mod_(2) << '\t'
            << G_mod_(3) << '\t' << G_mod_(4) << '\t' << G_mod_(5) << '\n';
  Matrix6x6 de = Matrix6x6::Zero();
  // clang-format off
  // compute elasticityTensor
  de(0,0)=a1;    de(0,1)=a2;    de(0,2)=a2;
  de(1,0)=a2;    de(1,1)=a1;    de(1,2)=a2;
  de(2,0)=a2;    de(2,1)=a2;    de(2,2)=a1;
  de(3,3)=G_mod_(3);
  de(4,4)=G_mod_(4);
  de(5,5)=G_mod_(5);
  // clang-format on
  return de;
}

//! Compute plastic tensor
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6>
    mpm::Hyperbolic<Tdim>::compute_elasto_plastic_tensor(
        const Vector6d& stress, const Vector6d& dstrain,
        const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars,
        bool hardening) {

  const Matrix6x6 dep = Matrix6x6::Zero();
  return dep;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Hyperbolic<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  auto strain = ptr->strain();

  if (((strain.cwiseEqual(0)).all()) or (init_dir_.cwiseEqual(0).all())) {
    // Determine directionality of initial strain
    init_dir_ = dstrain.cwiseSign();
    init_loading_ = true;
  } else {
    new_dir_ = dstrain.cwiseSign();
    if (new_dir_(3) != init_dir_(3)) {
      init_loading_ = false;
      reversal_strain_ = strain;
      reversal_stress_ = stress;
      if ((maximum_strain_.cwiseEqual(0)).all()) {
        maximum_strain_ = strain;
        secant_modulus_ = stress(3) / strain(3);
      } else {
        if ((abs(strain(3)) > abs(maximum_strain_(3)))) {
          maximum_strain_ = strain;
          secant_modulus_ = stress(3) / strain(3);
        }
      }
    }
  }
  Matrix6x6 de_ = compute_tensor(strain);

  const Vector6d dstress = de_ * dstrain;
  Vector6d b_mod = Vector6d::Zero();
  if (!init_loading_) {
    double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

    // Strain Vectors and subtracting reversal strain
    Vector6d strain_dif = (strain - reversal_strain_);

    // Calculating Reduction Factor
    double reduction_factor = (p1_ - p2_ * pow(1. - secant_modulus_ / G, p3_));

    // Calculating Second Term
    Vector6d term2 =
        (G * strain_dif.array() *
         (1. + beta_ * pow((strain_dif.array() + dstrain.array()).abs() /
                               (2 * reference_strain_),
                           s_))
             .inverse())
            .matrix();

    // Calculating Third Term
    Vector6d term3 =
        (G * strain_dif.array() *
         (1. +
          beta_ * pow((maximum_strain_.array() - reference_strain_).abs(), s_))
             .inverse())
            .matrix();

    // Calculate intercept value
    Vector6d b_mod = reduction_factor * (term2 - term3);
  }
  return (stress + dstress + b_mod);
}

//! Compute consistent tangent matrix
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6>
    mpm::Hyperbolic<Tdim>::compute_consistent_tangent_matrix(
        const Vector6d& stress, const Vector6d& prev_stress,
        const Vector6d& dstrain, const ParticleBase<Tdim>* ptr,
        mpm::dense_map* state_vars) {
  const Matrix6x6 de = this->compute_tensor(dstrain);
  return de;
}