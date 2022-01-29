//! Constructor with id and material properties
template <unsigned Tdim>
mpm::BoundSurfPlasticity<Tdim>::BoundSurfPlasticity(
    unsigned id, const Json& material_properties)
    : InfinitesimalElastoPlastic<Tdim>(id, material_properties) {
  try {
    // Input parameter: wr function d power
    d0_ = material_properties.at("d0").template get<double>();
    // Density
    density_ = material_properties.at("density").template get<double>();
    // Ratio fp = Rp / Rf
    fp_ = material_properties.at("fp").template get<double>();
    // Friction angle (input in degrees, saved in radians)
    friction_ =
        material_properties.at("friction").template get<double>() * M_PI / 180.;
    // Initial void ratio
    e0_ = material_properties.at("void_ratio_initial").template get<double>();
    // Shear modulus model parameter
    G0_ = material_properties.at("G0").template get<double>();
    // Input parameter: Ci and Cg functions gm coefficient
    gm0_ = material_properties.at("gm0").template get<double>();
    // Input parameter: Hr function hr coefficient
    hr0_ = material_properties.at("hr0").template get<double>();
    // Input parameter: wm function kr coefficient
    kr0_ = material_properties.at("kr0").template get<double>();
    // Poisson ratio
    poisson_ = material_properties.at("poisson_ratio").template get<double>();
    // Relative density (decimal)
    relative_density_ =
        material_properties.at("relative_density").template get<double>();
    if (relative_density_ > 1. or relative_density_ < 0.) {
      throw std::runtime_error("Relative density out of the bounds [0,1]!");
    }

    // OPTIONAL:: Critical state line void ratio intercept
    if (material_properties.find("void_ratio_intercept") !=
        material_properties.end()) {
      e_int_ =
          material_properties.at("void_ratio_intercept").template get<double>();
    } else {
      e_int_ = 0.9;
    }
    // OPTIONAL:: Input parameter: Ci function eta coefficient
    if (material_properties.find("eta0") != material_properties.end()) {
      eta0_ = material_properties.at("eta0").template get<double>();
    } else {
      eta0_ = 10.;
    }
    // OPTIONAL:: Critical state line slope
    if (material_properties.find("lambda") != material_properties.end()) {
      lambda_ = material_properties.at("lambda").template get<double>();
    } else {
      lambda_ = 0.03;
    }
    // OPTIONAL:: Reference (atomspheric) pressure
    if (material_properties.find("reference_pressure") !=
        material_properties.end()) {
      p_atm_ =
          material_properties.at("reference_pressure").template get<double>();
    } else {
      p_atm_ = 100.e+3;
    }
    // OPTIONAL:: Minimum allowable mean pressure
    if (material_properties.find("pmin") != material_properties.end()) {
      p_min_ = material_properties.at("pmin").template get<double>();
    } else {
      p_min_ = 10.;
    }
    // OPTIONAL:: Tolerance
    if (material_properties.find("tolerance") != material_properties.end()) {
      tolerance_ = material_properties.at("tolerance").template get<double>();
    }

    // Model parameter: wm function a power
    a_pow_ = 1.;
    // Model parameter: wm function b powers
    b_pow_ = 2.;
    // Model parameter: wr function d power
    d_ = scale(2.0, relative_density_) * d0_;
    // Cumulative plastic deviatoric strain set zero at start
    dep_ = 0.;
    // Model parameter: Ci function eta coefficient
    eta_ = eta0_;
    // Model parameter: Ci and Cg functions gm coefficient
    gm_ = scale(0.1, relative_density_) * gm0_;
    // Model parameter: Hr function hr coefficient
    hr_ = scale(1.5, relative_density_) * hr0_;
    // Model parameter: m function ka
    ka_ = std::pow(scale(2.0, relative_density_), 2);
    // Model parameter: Ci and Cg function ke power
    ke_ = scale(1.5, relative_density_) * 2.;
    // Model parameter: wm function kr coefficient
    kr_ = scale(2.0, relative_density_) * kr0_;
    // Model parameter: wr function ks power
    ks_ = 2. * std::max((relative_density_ - 0.5), 0.);
    // Critical state mean pressure
    pc_ = std::pow((e0_ - e_int_) / (-lambda_), (10. / 7.)) * p_atm_;
    // Input parameter: failure surface
    const double s_friction = sin(friction_);
    Rf0_ = std::sqrt(2.) * 2. * std::sqrt(3.) * s_friction / (3. - s_friction);
    // Failure surface
    Rf_ = scale(1.1, relative_density_) * Rf0_;
    // Maximum allowable R
    R_max_ = 0.99 * Rf_;
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

//! Compute Frobenius inner product
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::frobenius_prod(const Vector6d& vec_A,
                                                      const Vector6d& vec_B) {
  // Solves:
  // (A_ij B_ij) = (A_00*B_00 + A_11*B11 + A_22*B_22)
  //                + 2*(A_01*B_01 + A_12*B_12 + A_02*B_02)
  // Note that matrices are assumed symmtric and provided in [6x1] form
  double product = 0.;
  for (unsigned i = 0; i < 3; ++i) {
    product += (vec_A(i) * vec_B(i)) + (2. * (vec_A(i + 3) * vec_B(i + 3)));
  }
  return product;
}

//! Compute elastic tensor
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6>
    mpm::BoundSurfPlasticity<Tdim>::compute_elastic_tensor(
        const Vector6d& stress, mpm::dense_map* state_vars) {

  // Constants for elastic moduli
  const double mean_p = check_low(-1. * mpm::materials::p(stress));
  const double p_atm = this->p_atm_;
  const double G0 = this->G0_;
  const double e0 = this->e0_;
  const double nu = this->poisson_;

  // Elastic moduli
  const double G = p_atm * G0 * (std::pow((2.973 - e0), 2) / (1. + e0)) *
                   std::sqrt(mean_p / p_atm);
  const double K = (2. * G * (1. + nu)) / (3 * (1 - 2 * nu));

  // Matrix values
  const double a1 = K + (4.0 / 3.0) * G;
  const double a2 = K - (2.0 / 3.0) * G;
  const double S = 2. * G;

  // compute elastic stiffness matrix
  Matrix6x6 de = Matrix6x6::Zero();
  // clang-format off
  de(0,0)=a1;    de(0,1)=a2;    de(0,2)=a2;    de(0,3)=0;    de(0,4)=0;    de(0,5)=0;
  de(1,0)=a2;    de(1,1)=a1;    de(1,2)=a2;    de(1,3)=0;    de(1,4)=0;    de(1,5)=0;
  de(2,0)=a2;    de(2,1)=a2;    de(2,2)=a1;    de(2,3)=0;    de(2,4)=0;    de(2,5)=0;
  de(3,0)= 0;    de(3,1)= 0;    de(3,2)= 0;    de(3,3)=S;    de(3,4)=0;    de(3,5)=0;
  de(4,0)= 0;    de(4,1)= 0;    de(4,2)= 0;    de(4,3)=0;    de(4,4)=S;    de(4,5)=0;
  de(5,0)= 0;    de(5,1)= 0;    de(5,2)= 0;    de(5,3)=0;    de(5,4)=0;    de(5,5)=S;
  // clang-format on

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

  // Save stress as muteable vector
  Vector6d sigma = stress;

  // Note that dstrain (engineering) should be converted to dstrain (true)
  // for the Wang bounding surface plasticity model
  Vector6d dstrain_t = dstrain;
  for (unsigned i = 0; i < 3; ++i) dstrain_t(i + 3) *= 0.5;

  // Stress-based initial values
  if (first_loop_) {
    // Mean pressure
    double mean_p = check_low(-1. * mpm::materials::p(sigma));

    // Deviatoric stress ratio vector and invariant
    Vector6d dev_r = -1. * sigma / mean_p;
    for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
    const double R = std::sqrt((frobenius_prod(dev_r, dev_r)));

    // Save for later
    p_max_ = mean_p;
    alpha_v_ = dev_r;
    alpha_s_ = R;

    // Maximum pre-stress surface
    Rm_ = R;

    // Check minimum allowalbe mean pressure
    if (mean_p < p_min_) {
      mean_p = p_min_;
    }

    // State index parameter
    const double Ip = mean_p / pc_;

    // Dilation surface
    const double Rp = fp_ * Rf_;
    Rd_ = Rp + (Rf_ - Rp) * Ip;

    // Distance from reversal point to current stress state
    rho_ = R;

    // Distance from reversal point to yield surface
    rho_bar_ = Rm_;

    // Distance from reversal point to dilation surface
    rho_d_ = Rd_;

    // Set bool
    first_loop_ = false;
  }

  // Mean pressure
  double mean_p = check_low(-1. * mpm::materials::p(sigma));

  // Check minimum allowalbe mean pressure
  if (mean_p < p_min_) {
    mean_p = p_min_;
  }

  // Deviatoric stress ratio vector and invariant
  Vector6d dev_r = -1. * sigma / mean_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  double R = std::sqrt((frobenius_prod(dev_r, dev_r)));

  // Check if past maximum allowable R
  if (R > R_max_) {
    // Update stress
    sigma *= (R_max_ / R);
    for (unsigned i = 0; i < 3; ++i) sigma(i) -= (mean_p * (1. - (R_max_ / R)));

    // Update deviatoric stress ratio and invariant
    dev_r = -1. * sigma / mean_p;
    for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1;
    R = std::sqrt((frobenius_prod(dev_r, dev_r)));
  }

  // Normal to yield (pre-stress) surface
  Vector6d nl = Vector6d::Zero();

  // Check if past maximum pre-stress surface
  if (R > 0.999999 * Rm_) {
    // New maximum pre-stress surface
    Rm_ = R;

    // Normal to yield (pre-stress) surface
    nl = dev_r / R;

  } else {
    // NOT CHECKED /////////////////////////////////////////////////////////////
    // Compute rho
    Vector6d r_less_a = dev_r - alpha_v_;
    rho_ = std::sqrt((frobenius_prod(r_less_a, r_less_a)));

    // Compute rho_bar
    Vector6d a_less_r_norm = (alpha_v_ - dev_r).cwiseProduct(alpha_v_);
    double temp = 0.;
    for (unsigned i; i < 3; ++i)
      temp += (a_less_r_norm(i) + 2 * a_less_r_norm(i + 3));
    temp *= (1. / rho_);
    rho_bar_ = temp + std::sqrt(std::fabs(std::pow(temp, 2) + std::pow(Rm_, 2) -
                                          std::pow(alpha_s_, 2)));
    rho_d_ = temp + std::sqrt(std::fabs(std::pow(temp, 2) + std::pow(Rd_, 2) -
                                        std::pow(alpha_s_, 2)));

    // Project to yield (pre-stress) surface
    const Vector6d r_bar = alpha_v_ + (rho_bar_ / rho_) * r_less_a;
    const double r_bar_scalar = std::sqrt((frobenius_prod(r_bar, r_bar)));

    // Normal to yield (pre-stress) surface
    nl = r_bar / r_bar_scalar;
  }

  // Shear and bulk elastic moduli
  const double G = p_atm_ * G0_ * (std::pow((2.973 - e0_), 2) / (1. + e0_)) *
                   std::sqrt(mean_p / p_atm_);
  const double K = (2. * G * (1. + poisson_)) / (3 * (1 - 2 * poisson_));
  const double G_over_K = G / K;

  // Plastic moduli
  const double C_g = std::max(0.05, 1. / (1. + gm_ * std::pow(dep_, ke_)));
  const double C_i = (1. + gm_ * std::pow(dep_, ke_)) /
                     (1. + eta_ * gm_ * std::pow(dep_, ke_));

  double temp = std::max((rho_ / rho_bar_), 1.E-5);
  // double temp_a = (1. - std::exp(1. - mean_p / p_min_));
  double temp_b = std::min(temp, 1.);
  double temp_c = Rm_ / Rf_;
  double temp_d = mean_p / p_max_;

  // w_m parameter
  double w_m = 0.;
  if (kr_ < 100.) {
    if (fp_ < 1.) {
      w_m = std::pow(temp_d, a_pow_) * std::pow(temp_c, b_pow_) *
            ((Rd_ - Rm_) / (Rf_ - Rm_)) * (1. / kr_);
    } else {
      // NOT CHECKED ///////////////////////////////////////////////////////////
      w_m = std::pow(temp_d, a_pow_) * std::pow(temp_c, b_pow_) * (1. / kr_);
    }
  } else {
    // NOT CHECKED /////////////////////////////////////////////////////////////
    w_m = 0.;
  }

  // power m
  double m = std::pow(std::min(2. * Rm_ / rho_bar_, 50.), ka_);
  double tmp = std::pow(temp_b, m);
  double temp_e = 1. / std::max(tmp, 1.0E-10);

  // Hardening (?)
  double Hr = C_g * G * hr_ * (temp_e / temp_c - 1);
  double Hr_over_G = Hr / G;

  // w_r parameter
  double w_r = 0.;
  if (d_ < 100.) {
    double is = std::pow(std::max((p_max_ / p_atm_), 1.), ks_);
    w_r = is * std::pow(temp_c, d_) * (rho_d_ - rho_) * std::pow(rho_, m);
  } else {
    // NOT CHECKED /////////////////////////////////////////////////////////////
    w_r = 0.;
  }

  // pick w_m or w_r
  double factor = 0.;
  if (temp_b < 1.) {
    factor = std::pow(temp_b, 30.);
  } else {
    factor = 1.0;
  }
  const double w = w_r * (1. - factor) + w_m * factor;

  // std::cout << "\nw : " << w << std::endl;

  const double kri = w / K / C_i;  // probably delete
  const double kokr = w / C_i;

  // Elastic stress increment
  const Matrix6x6 de = this->compute_elastic_tensor(sigma, state_vars);
  const Vector6d dsigma_e = -1. * de * dstrain_t;
  const double d_mean_p_e = check_low(mpm::materials::p(dsigma_e));

  // Deviatoric strain increment (compression positive)
  const double vol_strain = -1. * (dstrain_t(0) + dstrain_t(1) + dstrain_t(2));
  Vector6d dev_dstrain = -1. * dstrain_t;
  for (unsigned i = 0; i < 3; ++i) dev_dstrain(i) -= (vol_strain / 3.);

  // dot
  double dot = frobenius_prod(dev_dstrain, dev_r);

  // drm
  double drm = std::sqrt(frobenius_prod(dev_dstrain, dev_dstrain));

  // sum
  double sum = frobenius_prod(dev_dstrain, nl);

  dot *= (1. / drm);

  // yn parameter (yield?)
  double y = std::max(0., (dot * dot) + (Rf_ * Rf_) - (R * R));
  y = -dot + std::sqrt(y);
  const double yn = y / drm;

  // Normal to FAILURE surface
  Vector6d normal = dev_r + (dev_dstrain * yn);
  double normal_temp = std::sqrt(frobenius_prod(normal, normal));
  normal *= (1. / normal_temp);

  //
  const double nd = Rm_ / Rf_;
  normal = ((1. - nd) * normal) + (nd * dev_dstrain / drm);
  normal_temp = std::sqrt(frobenius_prod(normal, normal));
  normal *= (1. / normal_temp);

  // dr_ij
  Vector6d dr_ij = Vector6d::Zero();

  // Either loading or unloading
  if (sum >= 0.) {
    // Loading continues
    // Updated maximum mean pressure
    if (mean_p > p_max_) p_max_ = mean_p;

    // Elasto-plastic coefficients
    const double Ar = 0.5 * Hr_over_G + 1.;
    const double Ap = Hr_over_G * G_over_K * kokr;
    const double Bp = 2. * G_over_K;
    const double Br = frobenius_prod(dev_r, normal);
    const double c2 = (Ar * Bp) - (Ap * Br);

    // Qp
    Vector6d Qp = normal * Bp;
    for (unsigned i; i < 3; ++i) Qp(i) -= Br;

    // Plastic increment
    // Qp
    double dQp = (-1. / c2) * frobenius_prod(dstrain_t, Qp);
    // Pr
    const double temp_pr = 0.5 * kokr * Hr_over_G;
    Vector6d Pr = normal;
    for (unsigned i; i < 3; ++i) Pr(i) += temp_pr;

    // Dp
    Vector6d Dp = 2. * G * Pr * dQp;

    // update stress
    sigma -= dsigma_e + Dp;

    //
    double dp = ((dsigma_e(0) + dsigma_e(1) + dsigma_e(2)) / 3.) -
                (2. * G * dQp * (Pr(0) + Pr(1) + Pr(2)) * (1. / 3.));
    mean_p = check_low(-1. * mpm::materials::p(sigma));

    //
    dr_ij = (dsigma_e - Dp + sigma * dp / mean_p) / mean_p;

  } else {
    // NOT CHECKED /////////////////////////////////////////////////////////////
    // Unloading begins
    // Reset projection center
    for (unsigned i; i < 6; ++i) alpha_v_(i) = -1. * (sigma(i) / mean_p) - 1.;
    alpha_s_ = std::sqrt(frobenius_prod(alpha_v_, alpha_v_));

    // Update stress
    sigma -= dsigma_e;
  }

  //
  double new_p = check_low(-1. * mpm::materials::p(sigma));
  if (new_p < p_min_) {
    for (unsigned i; i < 3; ++i) sigma(i) *= (p_min_ / new_p);
    new_p = p_min_;
  }

  //
  dev_r = -1. * sigma / new_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  R = std::sqrt((frobenius_prod(dev_r, dev_r)));

  //
  if (R > R_max_) {
    sigma *= (R_max_ / R);
    for (unsigned i = 0; i < 3; ++i) sigma(i) -= new_p * (1. - (R_max_ / R));
  }

  //
  const double dif = frobenius_prod(dr_ij, normal);
  const double devp = mean_p * dif / Hr;
  const double ddevp = std::sqrt(2.) * std::fabs(devp);

  // check
  double sum_dep = ddevp;
  if (sum_dep > 1.) sum_dep = 0.;

  //
  if (Rm_ - Rd_ > 0.) dep_ += sum_dep;

  // State index parameter
  double Ip = mean_p / pc_;

  // Dilation surface
  double Rp = fp_ * Rf_;
  Rd_ = Rp + (Rf_ - Rp) * Ip;

  // Return updated stress
  return (sigma);
}
