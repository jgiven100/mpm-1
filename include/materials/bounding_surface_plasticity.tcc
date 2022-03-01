//! Constructor with id and material properties
template <unsigned Tdim>
mpm::BoundSurfPlasticity<Tdim>::BoundSurfPlasticity(
    unsigned id, const Json& material_properties)
    : InfinitesimalElastoPlastic<Tdim>(id, material_properties) {
  try {
    // Input parameter: dz function cz coefficient
    cz_ = material_properties.at("cz").template get<double>();
    // Input parameter: wr function d power
    d0_ = material_properties.at("d0").template get<double>();
    // Density
    density_ = material_properties.at("density").template get<double>();
    // Critical void ratio
    ec_ = material_properties.at("critical_void_ratio").template get<double>();
    // Ratio fp = Rp / Rf
    fp_ = material_properties.at("fp").template get<double>();
    // Friction angle (input in degrees, saved in radians)
    friction_ =
        material_properties.at("friction").template get<double>() * M_PI / 180.;
    // Shear modulus model parameter
    G0_ = material_properties.at("G0").template get<double>();
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
    // Input parameter: zm function cz coefficient
    zm_ = material_properties.at("zm").template get<double>();

    // OPTIONAL:: Proportional loading flag
    if (material_properties.find("proportional") != material_properties.end()) {
      proportional_ =
          material_properties.at("proportional").template get<bool>();
    } else {
      proportional_ = true;
    }
    // OPTIONAL:: Critical state line void ratio intercept
    if (material_properties.find("void_ratio_intercept") !=
        material_properties.end()) {
      e0_ =
          material_properties.at("void_ratio_intercept").template get<double>();
    } else {
      e0_ = 0.9;
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
    // Model parameter: Hr function hr coefficient
    hr_ = scale(1.5, relative_density_) * hr0_;
    // Model parameter: m function ka power
    ka_ = std::pow(scale(2.0, relative_density_), 2);
    // Model parameter: wm function kr coefficient
    kr_ = scale(2.0, relative_density_) * kr0_;
    // Model parameter: wr function ks power
    ks_ = 2. * std::max((relative_density_ - 0.5), 0.);
    // Critical state mean pressure
    pc_ = std::pow((e0_ - ec_) / (lambda_), (10. / 7.)) * p_atm_;
    // Input parameter: failure surface
    const double s_friction = sin(friction_);
    Rf0_ = 2. * std::sqrt(3.) * s_friction / (3. - s_friction);
    // Failure surface
    Rf_ = scale(1.1, relative_density_) * Rf0_;
    // Maximum allowable R
    R_max_ = 0.99 * Rf_;
    // Dilatancy parameter zero at start
    z_ = 0.;
    // Set cos_a_ if proportional loading
    if (proportional_) cos_a_ = 1.;
    // Properties
    properties_ = material_properties;

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variabless
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

//! Compute Frobenius norm
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::frobenius_norm(const Vector6d& vec_A) {
  // Solves:
  // sqrt(A_ij A_ij) = sqrt[(A_00*A_00 + A_11*A_11 + A_22*A_22)
  //                        + 2*(A_01*A_01 + A_12*A_12 + A_02*A_02)]
  // Note that matrices are assumed symmtric and provided in [6x1] form
  double product = 0.;
  for (unsigned i = 0; i < 3; ++i) {
    product += (vec_A(i) * vec_A(i)) + (2. * (vec_A(i + 3) * vec_A(i + 3)));
  }
  const double norm = std::sqrt(product);
  return norm;
}

//! Compute Frobenius inner product
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::frobenius_prod(const Vector6d& vec_A,
                                                      const Vector6d& vec_B) {
  // Solves:
  // (A_ij B_ij) = (A_00*B_00 + A_11*B_11 + A_22*B_22)
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
  const double e0 = this->ec_;
  const double nu = this->poisson_;

  // Elastic moduli
  const double G = p_atm * G0 * (std::pow((2.973 - e0), 2) / (1. + e0)) *
                   std::sqrt(mean_p / p_atm);
  const double K = (2. * G * (1. + nu)) / (3 * (1 - (2 * nu)));

  // Matrix values
  const double a1 = K + (4.0 / 3.0) * G;
  const double a2 = K - (2.0 / 3.0) * G;
  const double S = 2. * G;

  // Compute elastic stiffness matrix
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

  // Note: dstrain (engineering shear strain) is converted to dstrain_t
  //       (true shear strain) for easier tensor multiplication in the Wang
  //       bounding surface plasticity model
  Vector6d dstrain_t = dstrain;
  for (unsigned i = 0; i < 3; ++i) dstrain_t(i + 3) *= 0.5;

  // Stress-based initial values
  if (first_loop_) {
    // Mean pressure
    double mean_p = check_low(-1. * mpm::materials::p(sigma));

    // Deviatoric stress ratio vector and invariant
    Vector6d dev_r = -1. * sigma / mean_p;
    for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
    const double R = std::sqrt(0.5) * frobenius_norm(dev_r);

    // Save for later
    p_max_ = mean_p;
    alpha_ = dev_r;

    // Maximum pre-stress surface
    Rm_ = R;

    // Check minimum allowalbe mean pressure
    if (mean_p < p_min_) {
      mean_p = p_min_;
    }

    // State index parameter
    const double Ip = mean_p / pc_;

    // Dilation surface
    const double Rp90 = 0.9 * Rf_;
    const double Rp0 = fp_ * Rf_;  // TODO : determine if correct
    const double Rp = Rp90 + (Rp0 - Rp90) * std::fabs(cos_a_);
    Rd_ = Rp + (Rf_ - Rp) * Ip;

    // Distance from reversal point to current stress state
    rho_ = std::sqrt(2.) * R;

    // Distance from reversal point to yield surface
    rho_bar_ = std::sqrt(2.) * Rm_;

    // Distance from reversal point to dilation surface
    rho_d_ = std::sqrt(2.) * Rd_;

    // Save stress state for incremental calculation in next loop
    last_sigma_ = stress;

    // Set bool
    first_loop_ = false;
  }

  // Mean pressure
  double mean_p = check_low(-1. * mpm::materials::p(sigma));

  // Check minimum allowalbe mean pressure
  if (mean_p < p_min_) {
    mean_p = p_min_;
  }

  // State index parameter
  const double Ip = mean_p / pc_;

  // Dilation surface
  const double Rp90 = 0.9 * Rf_;
  const double Rp0 = fp_ * Rf_;  // TODO : determine if correct
  const double Rp = Rp90 + (Rp0 - Rp90) * std::fabs(cos_a_);
  Rd_ = Rp + (Rf_ - Rp) * Ip;

  // Deviatoric stress ratio vector and invariant
  Vector6d dev_r = -1. * sigma / mean_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  double R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Check if past maximum allowable R
  if (R > R_max_) {
    // Update stress
    sigma *= (R_max_ / R);
    for (unsigned i = 0; i < 3; ++i) sigma(i) -= (mean_p * (1. - (R_max_ / R)));

    // Update deviatoric stress ratio and invariant
    dev_r = -1. * sigma / mean_p;
    for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1;
    R = std::sqrt(0.5) * frobenius_norm(dev_r);
  }

  // Normal to yield (pre-stress) surface
  Vector6d n_bar = Vector6d::Zero();

  // Check if past maximum pre-stress surface
  if (R > 0.999999 * Rm_) {
    // New maximum pre-stress surface
    Rm_ = R;

    // Normal to yield (pre-stress) surface
    n_bar = dev_r / frobenius_norm(dev_r);

  } else {
    // Compute rho
    const Vector6d w = dev_r - alpha_;
    rho_ = frobenius_norm(w);

    // Compute Frobenius norm of alpha
    const double A = frobenius_norm(alpha_);

    // Compute rho_bar
    const double v = frobenius_prod(-w, alpha_) / rho_;
    rho_bar_ = v + std::sqrt(std::fabs((v * v) + (2 * Rm_ * Rm_) - (A * A)));
    rho_d_ = v + std::sqrt(std::fabs((v * v) + (2 * Rd_ * Rd_) - (A * A)));

    // Project to yield (pre-stress) surface
    const Vector6d r_bar = alpha_ + (rho_bar_ / rho_) * w;

    // Normal to yield (pre-stress) surface
    n_bar = r_bar / frobenius_norm(r_bar);
  }

  // Elastic shear and bulk modulus
  const double G = p_atm_ * G0_ * (std::pow((2.973 - ec_), 2) / (1. + ec_)) *
                   std::sqrt(mean_p / p_atm_);
  const double K = (2. * G * (1. + poisson_)) / (3 * (1 - (2 * poisson_)));

  // Incremental stress
  Vector6d dsigma = sigma - this->last_sigma_;

  // Incremental mean stress
  const double mean_dp = check_low(-1. * mpm::materials::p(dsigma));

  // Incremental deviatoric stress
  Vector6d ds = -1. * dsigma;
  for (unsigned i = 0; i < 3; ++i) ds(i) -= mean_dp;

  // Incremental elastic strain
  Vector6d dstrain_e = (ds / (2. * G));
  for (unsigned i = 0; i < 3; ++i) dstrain_e(i) += (mean_dp / (3. * K));

  // Incremental plastic strain
  Vector6d dstrain_p = (-1. * dstrain_t) - dstrain_e;
  const double dstrain_vp =
      (1. / 3.) * (dstrain_p(0) + dstrain_p(1) + dstrain_p(2));

  // Fabric dilatancy
  const double sign = cos_a_ > 0 ? 1. : -1.;
  const double dsvp = dstrain_vp < 0 ? (-1. * dstrain_vp) : 0.;
  const double dz = cz_ * (sign * zm_ - z_) * dsvp;
  z_ += dz;
  const double Cz = 1 / (1 + (-1. * sign * z_));

  // Save stress state for incremental calculation in next loop
  last_sigma_ = stress;

  // Plastic shear modulus
  // TODO : think about me!!
  const double m = std::pow(std::min(2 * Rm_ / rho_bar_, 50.), ka_);

  const double rho_ratio = std::min(std::max((rho_ / rho_bar_), 1.E-5), 1.);
  const double rho_ratio_pow = 1. / std::max(std::pow(rho_ratio, m), 1.0E-10);

  double Hr = 0.;
  if (proportional_) {
    Hr = Cz * G * hr_ * (rho_ratio_pow * (Rf_ / Rm_) - 1);
  } else {
    Hr = Cz * G * hr_ / R;
  }

  // Plastic bulk modulus
  double w = 0;
  if (proportional_) {
    const double wm = (1. / kr_) * std::pow((mean_p / p_max_), a_pow_) *
                      std::pow((Rm_ / Rf_), b_pow_) *
                      ((Rd_ - Rm_) / (Rf_ - Rm_));
    const double wr = std::pow(std::max((p_max_ / p_atm_), 1.), ks_) *
                      std::pow((Rm_ / Rf_), d_) * (rho_d_ - rho_) *
                      std::pow(rho_, m);

    const double w_factor = (rho_ratio < 1.) ? std::pow(rho_ratio, 30.) : 1.0;
    w = wr * (1. - w_factor) + wm * w_factor;
  } else {
    w = (1. / kr_) * ((Rd_ - R) / (Rf_ - Rm_));
  }

  const double Kr = Cz * K / w;

  // Elastic stress increment
  const Matrix6x6 De = this->compute_elastic_tensor(sigma, state_vars);
  const Vector6d dsigma_e = -1. * De * dstrain_t;

  // Deviatoric strain increment
  const double vol_dstrain = -1. * (dstrain_t(0) + dstrain_t(1) + dstrain_t(2));
  Vector6d dev_dstrain = -1. * dstrain_t;
  for (unsigned i = 0; i < 3; ++i) dev_dstrain(i) -= (vol_dstrain / 3.);

  const double dev_dstrain_norm = frobenius_norm(dev_dstrain);
  const Vector6d dev_dstrain_unit = dev_dstrain / dev_dstrain_norm;

  // Compute scaling factor
  // TODO : better naming
  const double v = frobenius_prod(dev_dstrain_unit, dev_r);
  const double y =
      -v + std::sqrt(std::max(0., (v * v) + (2 * Rf_ * Rf_) - (2 * R * R)));
  const double y_norm = y / dev_dstrain_norm;

  // Project to failure surface
  const Vector6d r_caret = dev_r + (dev_dstrain * y_norm);

  // Normals to failure surface
  const Vector6d n_caret = r_caret / frobenius_norm(r_caret);
  const Vector6d n_tilde = dev_dstrain_unit;

  // Combine n_caret and n_tilde
  const Vector6d Ne = ((1. - (Rm_ / Rf_)) * n_caret) + ((Rm_ / Rf_) * n_tilde);
  const Vector6d n = Ne / frobenius_norm(Ne);

  // dev_dr
  Vector6d dev_dr = Vector6d::Zero();

  // Check loading status
  const bool loading = (frobenius_prod(dev_dstrain, n_bar) > 0.) ? true : false;

  // Either loading or unloading
  if (loading) {
    // Loading continues
    // Updated maximum mean pressure
    if (mean_p > p_max_) p_max_ = mean_p;

    // Plastic coefficients
    const double Ar = 2. * G / Hr + 1.;
    const double Ap = 2. * G / Kr;
    const double Bp = 2. * G / K;
    const double Br = frobenius_prod(dev_r, n);
    const double c2 = (Ar * Bp) - (Ap * Br);

    // Qp plastic tensor
    Vector6d Qp = n * Bp;
    for (unsigned i; i < 3; ++i) Qp(i) -= Br;

    // Qp plastic increment
    const double dQp = (-1. / c2) * frobenius_prod(Qp, dstrain_t);

    // Pr plastic tensor
    Vector6d Pr = (2. * G / Hr) * n;
    for (unsigned i; i < 3; ++i) Pr(i) += (K / Kr);

    // Plastic stiffness matrix
    Vector6d Dp = 2. * G * Pr * dQp;

    // Update stress
    sigma += (-dsigma_e + Dp);

    // Change in mean pressure
    double dp = ((dsigma_e(0) + dsigma_e(1) + dsigma_e(2)) / 3.) -
                (2. * G * dQp * (Pr(0) + Pr(1) + Pr(2)) * (1. / 3.));
    mean_p = check_low(-1. * mpm::materials::p(sigma));

    // Incremental deviatoric stress ratio vector
    dev_dr = (dsigma_e - Dp + sigma * dp / mean_p) / mean_p;

  } else {
    // Unloading begins
    // Reset projection center
    alpha_ = -1. * sigma / mean_p;
    for (unsigned i; i < 3; ++i) alpha_(i) -= 1.;

    // Update stress
    sigma -= dsigma_e;
  }

  // Check minimum allowalbe mean pressure
  double new_p = check_low(-1. * mpm::materials::p(sigma));
  if (new_p < p_min_) {
    for (unsigned i; i < 3; ++i) sigma(i) *= (p_min_ / new_p);
    new_p = p_min_;
  }

  // Deviatoric stress ratio vector and invariant
  dev_r = -1. * sigma / new_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Compute loading direction cos(alpha) for next step if bi-directional load
  if (proportional_) {
    if (loading) {
      cos_a_ = 1.;
    } else {
      cos_a_ = -1.;
    }
  } else {
    cos_a_ = frobenius_prod(dev_dr, dev_r) / frobenius_norm(dev_dr) /
             frobenius_norm(dev_r);
  }

  // Check if past maximum allowable R
  if (R > R_max_) {
    // Update stress
    sigma *= (R_max_ / R);
    for (unsigned i = 0; i < 3; ++i) sigma(i) -= (new_p * (1. - (R_max_ / R)));
  }

  // Loading index
  const double L = mean_p * frobenius_prod(dev_dr, n) / Hr;
  const double ddevp = std::sqrt(2.) * std::fabs(L);

  // If Hr too small, ddevp is too big; use 1. as arbitrary cutoff
  double sum_dep = ddevp;
  if (sum_dep > 1.) sum_dep = 0.;

  // Only consider strain if pre-stress surface is past dilation surface
  if (Rm_ - Rd_ > 0.) dep_ += sum_dep;

  // Return updated stress
  return (sigma);
}
