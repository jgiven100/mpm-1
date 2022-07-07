//! Read material properties
template <unsigned Tdim>
mpm::Hyperbolic<Tdim>::Hyperbolic(unsigned id, const Json& material_properties)
    : InfinitesimalElastoPlastic<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    // Input elastic moduli
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();
    // Compute other elastic moduli
    shear_modulus_ = youngs_modulus_ / (2. * (1. + poisson_ratio_));
    bulk_modulus_ = youngs_modulus_ / (3. * (1. - 2. * poisson_ratio_));

    // Material properties for viscoelastic damping
    beta_ = material_properties.at("beta").template get<double>();
    s_ = material_properties.at("s").template get<double>();
    strain_ref_ =
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
Eigen::Matrix<double, 6, 6> mpm::Hyperbolic<Tdim>::compute_elastic_tensor() {
  const double G = this->shear_modulus_;
  const double K = this->bulk_modulus_;
  const double a1 = K + (4.0 / 3.0) * G;
  const double a2 = K - (2.0 / 3.0) * G;

  Matrix6x6 de = Matrix6x6::Zero();
  // clang-format off
  // compute elasticityTensor
  de(0,0)=a1;    de(0,1)=a2;    de(0,2)=a2;
  de(1,0)=a2;    de(1,1)=a1;    de(1,2)=a2;
  de(2,0)=a2;    de(2,1)=a2;    de(2,2)=a1;
  de(3,3)=G;     de(4,4)=G;     de(5,5)=G;
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

  // Update current total strain
  this->strain_total_ += dstrain;

  // WARNING : hard coded for strain to be in 1-2 [Y-Z] plane
  const unsigned dir = 4;
  const double gamma = this->strain_total_(dir);

  // Check for new maximum strain
  if ((std::fabs(this->strain_max_(dir)) - std::fabs(gamma)) < 1.E-12) {
    this->strain_max_ = this->strain_total_;
  } else {
    this->init_loading_ = false;
  }

  // Compute incremental tangent shear factor 
  double factor = 0.;
  if (this->init_loading_) {
    factor = 1. / std::pow(1 + (std::fabs(gamma) / strain_ref_), 2);
  } else {
    const double gamma_rev = this->strain_max_(dir);
    factor = 4. / std::pow((std::fabs(gamma - gamma_rev) / strain_ref_) + 2, 2);
  }

  // Elastic stiffness tensor
  Matrix6x6 de = this->compute_elastic_tensor();

  // Updated shear component of stiffness matrix
  de(dir, dir) *= factor;

  // Compute updated stress
  const Vector6d updated_stress = stress + (de * dstrain);

  return updated_stress;
}
