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
    // Initial shear modulus
    G0_ = material_properties.at("G0").template get<double>();
    // Poisson ratio
    poisson_ = material_properties.at("poisson_ratio").template get<double>();
    // Initial void ratio
    e0_ = material_properties.at("initial_void_ratio").template get<double>();
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
    // fp ratio
    fp_ = material_properties.at("fp").template get<double>();
    
    // Critical state slope
    lambda_ = material_properties.at("lambda").template get<double>();
    // Critical void ratio
    ec_ = material_properties.at("critical_void_ratio").template get<double>();


    // Default reference pressure (atmospheric pressure, 100 kPa)
    if (material_properties.find("reference_pressure") != material_properties.end()) 
      p_atm_ = material_properties.at("reference_pressure").template get<double>();
    // Default tolerance
    if (material_properties.find("tolerance") != material_properties.end())
      tolerance_ = material_properties.at("tolerance").template get<double>();
    // Default minimum allowable mean pressure
    if (material_properties.find("pmin") != material_properties.end())
      p_min_ = material_properties.at("pmin").template get<double>();

    
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
    eta_ = eta0_;
    // Plastic shear modulus power #1
    ke_ = scale(1.5, relative_density_) * 2.; 
    // Plastic shear modulus power #2
    ka_ = std::pow(scale(2.0, relative_density_), 2);
    // Plastic bulk modulus power
    ks_ = 2. * std::max((relative_density_ - 0.5), 0.);
    // Maximum allowable R
    R_max_ = 0.99 * Rf_;
    // Default input a parameter
    ia_ = 1.;
    // Default b parameter
    b_ = 2.;


    // Critical state mean pressure
    pc_ = p_atm_ * std::pow(((ec_ - e0_) / lambda_), (10./7.));


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
                               {"yield_state", 0}};

  return state_vars;
}

//! Initialise state variables
template <unsigned Tdim>
std::vector<std::string> mpm::BoundSurfPlasticity<Tdim>::state_variables()
    const {
  const std::vector<std::string> state_vars = {"yield_state"};
  return state_vars;
}

//! Compute stress ratio invariant R
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_R(const Vector6d& dev_r) {

  const double R = std::sqrt(0.5 * (std::pow(dev_r(0), 2) + 
                                    std::pow(dev_r(1), 2) +
                                    std::pow(dev_r(2), 2) +
                                    2 * std::pow(dev_r(3), 2) +
                                    2 * std::pow(dev_r(4), 2) +
                                    2 * std::pow(dev_r(5), 2)));

  return R;
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

  Matrix6x6 dep = Matrix6x6::Zero();

  return dep;
}



//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::BoundSurfPlasticity<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Note that stress (tension positive) should be converted to stress_comp
  // (compression positive) for the Wang bounding surface plasticity model
  Vector6d stress_comp = -1. * stress;

  // Stress-based initial values
  if (first_loop_){
    // Mean pressure
    mean_p_ = check_low(mpm::materials::p(stress_comp));

    // Deviatoric stress ratio vector and invariant
    Vector6d dev_r = stress_comp / mean_p_;
    for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
    double R = compute_R(dev_r);

    // Save for later
    p_max_ = mean_p_;
    dev_a_vector_ = dev_r;    
    dev_a_scalar_ = R;
    Rm_ = R;

    // Check minimum allowalbe mean pressure
    if (mean_p_ < p_min_) { mean_p_ = p_min_; }

    // State index parameter
    double Ip = mean_p_ / pc_;
    
    // Dilation surface
    double Rp = fp_ * Rf_;
    Rd_ = Rp + (Rf_ - Rp) * Ip;

    // Rho
    rho_ = R;
    rho_bar_ = Rm_;

    // Set bool
    first_loop_ = false;    
  }

  // Mean pressure
  mean_p_ = check_low(mpm::materials::p(stress_comp));

  // Check minimum allowalbe mean pressure
  if (mean_p_ < p_min_) { mean_p_ = p_min_; }

  // Deviatoric stress ratio vector and invariant
  Vector6d dev_r = stress_comp / mean_p_;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  double R = compute_R(dev_r);

  if (R > R_max_) {
    // Scaling factor
    const double scale = R_max_ / R;

    // Update stress
    stress_comp *= scale;
    for (unsigned i = 0; i < 3; ++i) stress_comp(i) -= (mean_p_ * (1. - scale));

    // Update deviatoric stress ratio and invariant 
    dev_r = stress_comp / mean_p_;
    for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1;
    R = compute_R(dev_r);
  }


  Vector6d nl;
  Vector6d rb;
  double rb_scalar;

  // Flow direction
  if (R > 0.999999 * Rm_) {
    // New maximum pre-stress surface
    Rm_ = R;

    // Normal direction
    nl = dev_r / R;

    // 
    rb = dev_r;

    //
    rb_scalar = Rm_;

  } else {
    // Compute rho
    Vector6d r_less_a = dev_r - dev_a_vector_;
    rho_ = compute_R(r_less_a);
    
    // Compute rho_bar
    Vector6d a_less_r_norm = (dev_a_vector_ - dev_r).cwiseProduct(dev_a_vector_);
    double temp = 0.;
    for (unsigned i; i < 3; ++i) temp += (a_less_r_norm(i) + 2 * a_less_r_norm(i + 3));
    temp *= (1. / rho_);
    rho_bar_ = temp + std::sqrt(std::fabs(std::pow(temp, 2) + std::pow(Rm_, 2) - std::pow(dev_a_scalar_, 2)));
    rho_p_ = temp + std::sqrt(std::fabs(std::pow(temp, 2) + std::pow(Rd_ , 2) - std::pow(dev_a_scalar_, 2)));

    // Project
    rb = dev_a_vector_ + (rho_bar_ / rho_) * r_less_a;
    rb_scalar = compute_R(rb);

    // Normal
    nl = rb / rb_scalar;

  }


  // Vector6d temp0;
  // for (unsigned i; i < 6; ++i) temp0(i) = i*2;

  // Vector6d temp1;
  // for (unsigned i; i < 6; ++i) temp1(i) = i+3;

  // Vector6d temp2 = temp0.cwiseProduct(temp0) - temp1;
  // std::cout << "temp2(0) : " << temp2(0) << std::endl;
  // std::cout << "temp2(1) : " << temp2(1) << std::endl;
  // std::cout << "temp2(2) : " << temp2(2) << std::endl;
  // std::cout << "temp2(3) : " << temp2(3) << std::endl;
  // std::cout << "temp2(4) : " << temp2(4) << std::endl;
  // std::cout << "temp2(5) : " << temp2(5) << std::endl;














  Matrix6x6 de = Matrix6x6::Zero();

  
  // Return updated stress
  return (de * dstrain);
}
