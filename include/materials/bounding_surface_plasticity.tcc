//! Constructor with id and material properties
template <unsigned Tdim>
mpm::BoundSurfPlasticity<Tdim>::BoundSurfPlasticity(
    unsigned id, const Json& material_properties)
    : InfinitesimalElastoPlastic<Tdim>(id, material_properties) {
  try {
    // Density
    density_ = material_properties.at("density").template get<double>();
    // Friction angle (input in degrees, saved in radians)
    friction_ = material_properties.at("friction").template get<double>() * M_PI / 180.;
    // Poisson ratio
    poisson_ = material_properties.at("poisson_ratio").template get<double>();
    // Relative density (decimal)
    relative_density_ = material_properties.at("relative_density").template get<double>();
    if (relative_density_ > 1. or relative_density_ < 0.) {
      throw std::runtime_error("Relative density out of the bounds [0,1]!");
    }
    // Initial plastic shear modulus parameter
    hr0_ = material_properties.at("hr0").template get<double>();
    // Initial plastic bulk modulus parameter
    kr0_ = material_properties.at("kr0").template get<double>();
    // Initial plastic bulk modulus power
    d0_ = material_properties.at("d0").template get<double>();
    // Initial UNKNOWN gamma parameter
    gamma0_ = material_properties.at("unknown_gamma_parameter").template get<double>(); // UNKNOWN
    // Initial UNKNOWN eta (maybe ita) parameter
    eta0_ = material_properties.at("unknown_eta_parameter").template get<double>(); // UNKNOWN
    // Declared p min
    pmin_ = material_properties.at("pmin").template get<double>();
    // Initial maximum stress ratio
    Rm0_ = material_properties.at("Rm0").template get<double>();
    // Initial stress ratio
    R0_ = material_properties.at("R0").template get<double>();
    // Critical mean pressure
    pc_ = material_properties.at("pc").template get<double>();
    // Reference pressure
    p_atm_ = material_properties.at("reference_pressure").template get<double>();
    // Initial mean pressure
    p0_ = material_properties.at("p0").template get<double>();
    // fp ratio
    fp_ = material_properties.at("fp").template get<double>();




    // Initial failure surface
    const double sin_friction = sin(friction_);
    Rf0_ = 2. * std::sqrt(3.) * sin_friction / (3. - sin_friction);
    // Failure surface
    Rf_ = scale(1.1, relative_density_) * Rf0_;
    // Plastic shear modulus parameter
    hr_ = scale(1.5, relative_density_) * hr0_;
    // Plastic bulk modulus parameter
    kr_ = scale(2.0, relative_density_) * kr0_;
    // Plastic bulk modulus power
    d_ = scale(2.0, relative_density_) * d0_;
    // UNKNOWN gamma parameter
    gamma_ = scale(0.1, relative_density_) * gamma0_;
    // UNKNOWN eta (maybe ita) parameter
    // eta_ = scale(0.05, relative_density_) * eta0_;                           // commented out in code
    eta_ = eta0_;
    // Plastic shear modulus power #1
    ke_ = scale(1.5, relative_density_); 
    // Plastic shear modulus power #2
    ka_ = std::pow(scale(2.0, relative_density_), 2);
    // Plastic bulk modulus power
    ks_ = 2. * std::max((relative_density_ - 0.5), 0.);

    // Rmax
    Rmax_ = 0.99 * Rf_;

    // Default a parameter
    a_ = 1.;
    // Default b parameter
    b_ = 2.;

    // Rp0
    Rp0_ = Rf0_ * fp_;
    // Rp90
    Rp90_ = 0.9 * Rf0_;


    // TODO : remove me ////////////////////////////////////////////////////////
    std::cout << "\n\ndensity_ : " << density_ << std::endl;
    std::cout << "friction_ : " << friction_ << std::endl;
    std::cout << "poisson_ : " << poisson_ << std::endl;
    std::cout << "relative_density_ : " << relative_density_ << std::endl;
    std::cout << "hr0_ : " << hr0_ << std::endl;
    std::cout << "kr0_ : " << kr0_ << std::endl;
    std::cout << "d0_ :  " << d0_ << std::endl;
    std::cout << "gamma0_ :  " << gamma0_ << std::endl;
    std::cout << "eta0_ :  " << eta0_ << std::endl;
    std::cout << "pmin_ :  " << pmin_ << std::endl;
    std::cout << "Rm0_ :  " << Rm0_ << std::endl;
    std::cout << "R0_ :  " << R0_ << std::endl;
    std::cout << "pc_ :  " << pc_ << std::endl;
    std::cout << "p_atm_ :  " << p_atm_ << std::endl;
    std::cout << "p0_ :  " << p0_ << std::endl;
    std::cout << "fp_ :  " << fp_ << std::endl;







    std::cout << "\n\nRf0_ : " << Rf0_ << std::endl;
    std::cout << "Rf_ : " << Rf_ << std::endl;
    std::cout << "hr_ : " << hr_ << std::endl;
    std::cout << "kr_ : " << kr_ << std::endl;
    std::cout << "d_ :  " << d_ << std::endl;
    std::cout << "gamma_ :  " << gamma_ << std::endl;
    std::cout << "eta_ :  " << eta_ << std::endl;
    std::cout << "ke_ :  " << ke_ << std::endl;
    std::cout << "ka_ :  " << ka_ << std::endl;
    std::cout << "ks_ :  " << ks_ << std::endl;


    std::cout << "\n\nRmax_ : " << Rmax_ << std::endl;
    std::cout << "a_ : " << a_ << std::endl;
    std::cout << "b_ : " << b_ << std::endl;
    std::cout << "Rp0_ : " << Rp0_ << std::endl;
    std::cout << "Rp90_ : " << Rp90_ << std::endl;








    std::cout << "\n\n" << std::endl;
    ////////////////////////////////////////////////////////////////////////////

    // Properties
    properties_ = material_properties;

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::BoundSurfPlasticity<Tdim>::initialise_state_variables() {

  mpm::dense_map state_vars = {// Yield state: 0: elastic, 1: yield
                               {"yield_state", 0},
                               // Rho bar
                               {"rho_bar", Rm0_},
                               // Rho
                               {"rho", R0_},
                               // State pressure index
                               {"Ip", (p0_ / pc_)},
                               // Equivalent plastic deviatoric strain
                               {"pdstrain", 0.},
                               // Plastic strain components
                               {"plastic_strain0", 0.},
                               {"plastic_strain1", 0.},
                               {"plastic_strain2", 0.},
                               {"plastic_strain3", 0.},
                               {"plastic_strain4", 0.},
                               {"plastic_strain5", 0.}};

  return state_vars;
}

//! Initialise state variables
template <unsigned Tdim>
std::vector<std::string> mpm::BoundSurfPlasticity<Tdim>::state_variables()
    const {
  const std::vector<std::string> state_vars = {
      "yield_state",     "rho_bar",         "rho",
      "Ip",              "pdstrain",
      "plastic_strain0", "plastic_strain1", "plastic_strain2",
      "plastic_strain3", "plastic_strain4", "plastic_strain5"};
  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6>
    mpm::BoundSurfPlasticity<Tdim>::compute_elastic_tensor(
        const Vector6d& stress, mpm::dense_map* state_vars) {

  Matrix6x6 de = Matrix6x6::Zero();

  return de;
}

//! Compute plastic tensor
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6>
    mpm::BoundSurfPlasticity<Tdim>::compute_elasto_plastic_tensor(
        const Vector6d& stress, const Vector6d& dstrain,
        const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars,
        bool hardening) {

  mpm::boundsurfplasticity::FailureState yield_type =
      yield_type_.at(int((*state_vars).at("yield_state")));
  // Return the updated stress in elastic state
  const Matrix6x6 de = this->compute_elastic_tensor(stress, state_vars);
  if (yield_type == mpm::boundsurfplasticity::FailureState::Elastic) {
    return de;
  }

  // Construct dep matrix
  Matrix6x6 dep = de;
  return dep;
}



//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::BoundSurfPlasticity<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Compute alpha
  const double alpha = 0.;
  const double cos_alpha = cos(alpha);


  // Rp parameter
  const double Rp = Rp90_ + (Rp0_ - Rp90_) * std::fabs(cos_alpha);

  // Dilation surface
  const double Rd = Rp + (Rf_ - Rp) * (*state_vars).at("Ip");



  // Note that stress (tension positive) should be converted to stress_neg
  // (compression positive) for the Wang bounding surface plasticity model
  Vector6d stress_neg = -1 * stress;
  const double temp_p = mpm::materials::p(stress_neg);
  if (temp_p < pmin_) { const double mean_p = pmin_; }
  else { const double mean_p = temp_p; }











  Matrix6x6 de = Matrix6x6::Zero();

  
  // Return updated stress
  return (de * dstrain);
}
