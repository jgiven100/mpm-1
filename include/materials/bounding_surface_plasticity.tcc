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
    // Maximum void ratio
    e_max_ =
        material_properties.at("maximum_void_ratio").template get<double>();
    // Minimum void ratio
    e_min_ =
        material_properties.at("minimum_void_ratio").template get<double>();
    // Ratio fp = Rp / Rf
    fp_ = material_properties.at("fp").template get<double>();
    // Friction angle (input in degrees, saved in radians)
    friction_ =
        material_properties.at("friction").template get<double>() * M_PI / 180.;
    // Shear modulus parameter
    G0_ = material_properties.at("G0").template get<double>();
    // Input parameter: Ci and Cg functions gm coefficient
    gm0_ = material_properties.at("gm0").template get<double>();
    // Input parameter: Hr function hr coefficient
    hr0_ = material_properties.at("hr0").template get<double>();
    // Input parameter: wm function kr coefficient
    kr0_ = material_properties.at("kr0").template get<double>();
    // Poisson ratio
    poisson_ = material_properties.at("poisson_ratio").template get<double>();
    // Initial void ratio (critical void ratio if undrained)
    void_ratio_initial_ =
        material_properties.at("initial_void_ratio").template get<double>();

    // OPTIONAL:: Critical state line void ratio intercept
    if (material_properties.find("void_ratio_intercept") !=
        material_properties.end()) {
      e0_ =
          material_properties.at("void_ratio_intercept").template get<double>();
    } else {
      e0_ = 0.9;
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
    // OPTIONAL:: Maximum allowable mean pressure
    if (material_properties.find("pmax") != material_properties.end()) {
      p_max_ = material_properties.at("pmin").template get<double>();
    } else {
      p_max_ = 1.e+7;
    }
    // OPTIONAL:: Minimum allowable mean pressure
    if (material_properties.find("pmin") != material_properties.end()) {
      p_min_ = material_properties.at("pmin").template get<double>();
    } else {
      p_min_ = 1.e-5;
    }
    // OPTIONAL:: Tolerance
    if (material_properties.find("tolerance") != material_properties.end()) {
      tolerance_ = material_properties.at("tolerance").template get<double>();
    }

    // Model parameter: wm function a power
    a_pow_ = 1.;
    // Model parameter: wm function b powers
    b_pow_ = 2.;
    // Input parameter: failure surface
    Rf0_ = 2. * std::sqrt(3.) * sin(friction_) / (3. - sin(friction_));
    // Properties
    properties_ = material_properties;

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::BoundSurfPlasticity<Tdim>::initialise_state_variables() {

  mpm::dense_map state_vars = {
      // Yield state: 0: elastic, 1: yield
      {"yield_state", 0},
      // Current void ratio
      {"void_ratio", this->void_ratio_initial_},
      // Cumulative plastic deviatoric strain
      {"dep", 0.},
      // Stress reversal point
      {"alpha0", 0.},
      {"alpha1", 0.},
      {"alpha2", 0.},
      {"alpha3", 0.},
      {"alpha4", 0.},
      {"alpha5", 0.},
      // First loop: 1: first loop, 0: else
      {"first_loop", 1},
      // Maximum mean pressure
      {"p_max", 0.},
      // Maximum pre-stress surface
      {"Rm", 0.},
      // Distance from reversal point to current state
      {"rho", 0.},
      // Distance from reversal point to yield surface
      {"rho_bar", 0.},
      // Distance from reversal point to dilation surface
      {"rho_d", 0.},
      // Equivalent plastic deviatoric strain
      {"pdstrain", 0.}};
  return state_vars;
}

//! Initialise state variables
template <unsigned Tdim>
std::vector<std::string> mpm::BoundSurfPlasticity<Tdim>::state_variables()
    const {
  const std::vector<std::string> state_vars = {
      "yield_state", "void_ratio", "dep",    "alpha0",     "alpha1", "alpha2",
      "alpha3",      "alpha4",     "alpha5", "first_loop", "p_max",  "Rm",
      "rho",         "rho_bar",    "rho_d",  "pdstrain"};
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

  // Elastic moduli
  const double G = compute_G(mean_p, state_vars);
  const double K = (2. * G * (1. + poisson_)) / (3. * (1. - 2. * poisson_));

  // Matrix values
  const double a1 = K + (4.0 / 3.0) * G;
  const double a2 = K - (2.0 / 3.0) * G;

  // Compute elastic stiffness matrix
  Matrix6x6 de = Matrix6x6::Zero();
  // clang-format off
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
    mpm::BoundSurfPlasticity<Tdim>::compute_elasto_plastic_tensor(
        const Vector6d& stress, const Vector6d& dstrain,
        const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars,
        bool hardening) {

  // Get power terms
  const double ka = density_vars_(4);

  // Get distance from stress reversal to current and max pre-stress surfaces
  const double rho = (*state_vars).at("rho");
  const double rho_bar = (*state_vars).at("rho_bar");

  // Get current surfaces
  const double Rm = (*state_vars).at("Rm");

  // Mean stress
  const double mean_p = check_low(-1. * mpm::materials::p(stress));

  // Elastic moduli
  const double G = compute_G(mean_p, state_vars);
  const double K = (2. * G * (1. + poisson_)) / (3. * (1. - (2. * poisson_)));

  // Compute dilation surface
  const double Rd = compute_Rd(mean_p, state_vars);

  // Deviatoric stress ratio vector and invariant
  Vector6d dev_r = -1. * stress / mean_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  double R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Normal to yield surface (incremental strain dependent)
  const Vector6d n = compute_n(dev_r, dstrain, state_vars);

  // Plastic shear modulus
  const double m = compute_m(state_vars);
  const double Hr = compute_Hr(mean_p, state_vars);

  // Plastic bulk modulus
  const double Ci = compute_Ci(state_vars);
  const double w = compute_w(mean_p, m, Rd, state_vars);
  const double Kr = Ci * K / w;

  // Plastic coefficients
  const double Ar = 2. * G / Hr + 1.;
  const double Ap = 2. * G / Kr;
  const double Bp = 2. * G / K;
  const double Br = frobenius_prod(dev_r, n);
  const double c = (Ar * Bp) - (Ap * Br);

  // Qp plastic tensor
  Matrix6x6 Qp = Matrix6x6::Zero();
  for (unsigned i; i < 6; ++i) {
    Qp(i, 0) = n(0) * Bp - Br;
    Qp(i, 1) = n(1) * Bp - Br;
    Qp(i, 2) = n(2) * Bp - Br;
    Qp(i, 3) = n(3) * Bp;
    Qp(i, 4) = n(4) * Bp;
    Qp(i, 5) = n(5) * Bp;
  }

  // Pr plastic tensor
  Matrix6x6 Pr = Matrix6x6::Zero();
  for (unsigned i; i < 6; ++i) {
    Pr(i, i) = (2. * G / Hr) * n(i);
  }
  for (unsigned i; i < 3; ++i) {
    Pr(i, i) += (K / Kr);
  }

  // Plastic stiffness tensor
  const Matrix6x6 Dp = 2. * G * Pr * Qp * (1. / c);

  // Elastic stress tensor
  const Matrix6x6 De = this->compute_elastic_tensor(stress, state_vars);

  // Combine tensors
  const Matrix6x6 Dep = De - Dp;

  return Dep;
}

//! Compute necessary paramters for first loop
template <unsigned Tdim>
void mpm::BoundSurfPlasticity<Tdim>::compute_first_loop(
    const Vector6d& stress, mpm::dense_map* state_vars) {

  // Get params
  const double Rf = density_vars_(8);

  // Mean pressure
  double mean_p = check_low(-1. * mpm::materials::p(stress));

  // Deviatoric stress ratio vector and invariant
  Vector6d dev_r = -1. * stress / mean_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  const double R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Save for later
  (*state_vars).at("p_max") = mean_p;

  // Save stress reversal point
  (*state_vars).at("alpha0") = dev_r(0);
  (*state_vars).at("alpha1") = dev_r(1);
  (*state_vars).at("alpha2") = dev_r(2);
  (*state_vars).at("alpha3") = dev_r(3);
  (*state_vars).at("alpha4") = dev_r(4);
  (*state_vars).at("alpha5") = dev_r(5);

  // Maximum pre-stress surface
  (*state_vars).at("Rm") = R;

  // Check minimum allowalbe mean pressure
  if (mean_p < p_min_) mean_p = p_min_;

  // Set distance from stress reversal point to current stress state
  (*state_vars).at("rho") = std::sqrt(2.) * R;

  // Set distance from stress reversal point to yield surface
  (*state_vars).at("rho_bar") = std::sqrt(2.) * R;

  // set distance from stress reversal point to dilation surface
  const double Rd = compute_Rd(mean_p, state_vars);
  (*state_vars).at("rho_d") = std::sqrt(2.) * Rd;

  // Set bool
  (*state_vars).at("first_loop") = 0;
}

//! Compute Cg coefficient (plastic shear modulus)
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_Cg(mpm::dense_map* state_vars) {

  // Get density dependent params
  const double gm = density_vars_(2);
  const double ke = density_vars_(5);

  // Get history dependent params
  const double dep = (*state_vars).at("dep");

  // Plastic shear modulus coefficient
  const double Cg = std::max(0.05, 1. / (1. + gm * std::pow(dep, ke)));

  return Cg;
}

//! Compute Ci coefficient (plastic bulk modulus)
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_Ci(mpm::dense_map* state_vars) {

  // Get density dependent params
  const double eta = density_vars_(1);
  const double gm = density_vars_(2);
  const double ke = density_vars_(5);

  // Get history dependent params
  const double dep = (*state_vars).at("dep");

  // Plastic bulk modulus coefficient
  const double Ci =
      (1. + gm * std::pow(dep, ke)) / (1. + eta * gm * std::pow(dep, ke));

  return Ci;
}

//! Compute plastic shear modulus G
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_G(const double& mean_p,
                                                 mpm::dense_map* state_vars) {

  // Get inputs
  const double G0 = this->G0_;
  const double void_ratio = (*state_vars).at("void_ratio");

  // Plastic shear modulus
  const double G = p_atm_ * G0 *
                   (std::pow((2.973 - void_ratio), 2) / (1. + void_ratio)) *
                   std::sqrt(mean_p / p_atm_);

  return G;
}

//! Compute plastic shear modulus coefficient Hr
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_Hr(const double& mean_p,
                                                  mpm::dense_map* state_vars) {

  // Get density dependent params
  const double hr = density_vars_(3);

  // Get distance from stress reversal to current and max pre-stress surfaces
  const double rho = (*state_vars).at("rho");
  const double rho_bar = (*state_vars).at("rho_bar");

  // Get current surfaces
  const double Rf = density_vars_(8);
  const double Rm = (*state_vars).at("Rm");

  // Elastic moduli
  const double G = compute_G(mean_p, state_vars);

  // m power
  const double m = compute_m(state_vars);

  // rho ratios
  const double rho_ratio = std::min(std::max((rho / rho_bar), 1.E-5), 1.);
  const double rho_ratio_pow = 1. / std::max(std::pow(rho_ratio, m), 1.0E-10);

  // Plastic shear modulus coefficient
  const double Cg = compute_Cg(state_vars);

  // Plastic shear modulus
  const double Hr = Cg * G * hr * (rho_ratio_pow * (Rf / Rm) - 1);

  return Hr;
}

//! Compute m power
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_m(mpm::dense_map* state_vars) {

  // Get density dependent params
  const double ka = density_vars_(4);

  // Get params
  const double Rm = (*state_vars).at("Rm");
  const double rho_bar = (*state_vars).at("rho_bar");

  // TODO : think about the sqrt(2) here!!
  const double m =
      std::pow(std::min(2 * std::sqrt(2.) * Rm / rho_bar, 50.), ka);

  return m;
}

//! Compute normal to yield surface (dstrain)
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::BoundSurfPlasticity<Tdim>::compute_n(
    const Vector6d& dev_r, const Vector6d& dstrain,
    mpm::dense_map* state_vars) {

  // Get density dependent params
  const double Rf = density_vars_(8);

  // Get params
  const double Rm = (*state_vars).at("Rm");

  // Deviatoric stress invariant
  const double R = std::sqrt(0.5) * frobenius_norm(dev_r);

  //! Convert strain from engineering to true shear strain
  const Vector6d dstrain_t = convert_strain(dstrain);

  // Deviatoric strain (true shear) increment
  const double vol_dstrain = -1. * trace(dstrain_t);
  Vector6d dev_dstrain = -1. * dstrain_t;
  for (unsigned i = 0; i < 3; ++i) dev_dstrain(i) -= (vol_dstrain / 3.);

  // Useful norm and unit normal
  const double dev_dstrain_norm = frobenius_norm(dev_dstrain);
  const Vector6d dev_dstrain_unit = dev_dstrain / dev_dstrain_norm;

  // Compute scaling factor
  // TODO : better naming
  const double v = frobenius_prod(dev_dstrain_unit, dev_r);
  const double y =
      -v + std::sqrt(std::max(0., (v * v) + (2 * Rf * Rf) - (2 * R * R)));
  const double y_norm = y / dev_dstrain_norm;

  // Project to failure surface
  const Vector6d r_caret = dev_r + (dev_dstrain * y_norm);

  // Normals to failure surface
  const Vector6d n_caret = r_caret / frobenius_norm(r_caret);
  const Vector6d n_tilde = dev_dstrain_unit;

  // Combine n_caret & n_tilde and normalize
  const Vector6d Ne = ((1. - (Rm / Rf)) * n_caret) + ((Rm / Rf) * n_tilde);
  const Vector6d n = Ne / frobenius_norm(Ne);

  return n;
}

//! Compute normal to yield (pre-stress) surface
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::BoundSurfPlasticity<Tdim>::compute_n_bar(
    const Vector6d& dev_r, const double& Rd, mpm::dense_map* state_vars) {

  // Get history dependent params
  Vector6d alpha = Vector6d::Zero();
  alpha(0) = (*state_vars).at("alpha0");
  alpha(1) = (*state_vars).at("alpha1");
  alpha(2) = (*state_vars).at("alpha2");
  alpha(3) = (*state_vars).at("alpha3");
  alpha(4) = (*state_vars).at("alpha4");
  alpha(5) = (*state_vars).at("alpha5");

  const double rho_bar = (*state_vars).at("rho_bar");
  const double Rm = (*state_vars).at("Rm");

  // Deviatoric stress invariant
  const double R = sqrt(0.5) * frobenius_norm(dev_r);

  // Normal to yield (pre-stress) surface
  Vector6d n_bar = Vector6d::Zero();

  // Check if past maximum pre-stress surface
  if (R > 0.999999 * Rm) {
    // Set updated max (pre-stress) surface
    (*state_vars).at("Rm") = R;

    // Normal to yield (pre-stress) surface
    n_bar = dev_r / frobenius_norm(dev_r);

  } else {
    // Compute rho
    const Vector6d renameme = dev_r - alpha;  // TODO : better name here
    const double rho = frobenius_norm(renameme);

    // Compute Frobenius norm of alpha
    const double A = frobenius_norm(alpha);

    // Compute rho_bar
    const double v = frobenius_prod(-renameme, alpha) / rho;
    const double rho_bar =
        v + std::sqrt(std::fabs((v * v) + (2 * Rm * Rm) - (A * A)));
    const double rho_d =
        v + std::sqrt(std::fabs((v * v) + (2 * Rd * Rd) - (A * A)));

    // Project to yield (pre-stress) surface
    const Vector6d r_bar = alpha + (rho_bar / rho) * renameme;

    // Set updated distances
    (*state_vars).at("rho") = rho;
    (*state_vars).at("rho_bar") = rho_bar;
    (*state_vars).at("rho_d") = rho_d;

    // Normal to yield (pre-stress) surface
    n_bar = r_bar / frobenius_norm(r_bar);
  }

  return n_bar;
}

//! Compute dilation surface
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_Rd(const double& mean_p,
                                                  mpm::dense_map* state_vars) {

  // Get dilation params
  const double fp = this->fp_;
  const double Rf = density_vars_(8);

  // Get critical pressure parameters
  const double void_ratio = (*state_vars).at("void_ratio");
  const double void_ratio_intercept = this->e0_;
  const double lambda = this->lambda_;

  // Critical mean pressure
  double pc = 0.;
  if (void_ratio > void_ratio_intercept) {
    pc = mean_p;
  } else {
    pc = std::pow((void_ratio_intercept - void_ratio) / (lambda), (10. / 7.)) *
         p_atm_;
  }

  // State index parameter
  const double Ip = mean_p / pc;

  // Dilation surface
  const double Rp = fp * Rf;
  const double Rd = Rp + (Rf - Rp) * Ip;

  return Rd;
}

//! Compute relative density
template <unsigned Tdim>
void mpm::BoundSurfPlasticity<Tdim>::compute_relative_density_terms(
    mpm::dense_map* state_vars) {

  // Get void ratio
  double void_ratio = (*state_vars).at("void_ratio");

  // Relative density
  double relative_density = (e_max_ - void_ratio) / (e_max_ - e_min_);

  if (relative_density > 0.95) {
    relative_density = 0.95;
  } else if (relative_density < 0.10) {
    relative_density = 0.10;
  } else {
    relative_density = relative_density;
  }

  // Model parameter: wr function d power
  density_vars_(0) = scale(2.0, relative_density) * d0_;
  // Model parameter: Ci function eta coefficient
  density_vars_(1) = eta0_;
  // Model parameter: Ci and Cg functions gm coefficient
  density_vars_(2) = scale(0.1, relative_density) * gm0_;
  // Model parameter: Hr function hr coefficient
  density_vars_(3) = scale(1.5, relative_density) * hr0_;
  // Model parameter: m function ka power
  density_vars_(4) = std::pow(scale(2.0, relative_density), 2);
  // Model parameter: Ci and Cg function ke power
  density_vars_(5) = scale(1.5, relative_density) * 2.;
  // Model parameter: wm function kr coefficient
  density_vars_(6) = scale(2.0, relative_density) * kr0_;
  // Model parameter: wr function ks power
  density_vars_(7) = 2. * std::max((relative_density - 0.5), 0.);
  // Failure surface
  // const double Rf = scale(1.1, relative_density) * Rf0_;
  const double Rf = scale(1.0, relative_density) * Rf0_;
  density_vars_(8) = Rf;
  // Maximum allowable R
  density_vars_(9) = 0.99 * Rf;
}

//! Compute plastic bulk modulus coefficient w
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_w(const double& mean_p,
                                                 const double& m,
                                                 const double& Rd,
                                                 mpm::dense_map* state_vars) {

  // Get power terms
  const double a = this->a_pow_;
  const double b = this->b_pow_;

  // Get density dependent params
  const double d = density_vars_(0);
  const double ks = density_vars_(7);
  const double kr = density_vars_(6);

  // Get coefficients
  const double fp = this->fp_;

  // Get distance from stress reversal to current, max pre-stress, and dilation
  // surfaces
  const double rho = (*state_vars).at("rho");
  const double rho_bar = (*state_vars).at("rho_bar");
  const double rho_d = (*state_vars).at("rho_d");

  // Get current surfaces
  const double Rf = density_vars_(8);
  const double Rm = (*state_vars).at("Rm");

  //
  const double p_max = (*state_vars).at("p_max");

  // Useful rho ratio
  const double rho_ratio = std::min(std::max((rho / rho_bar), 1.E-5), 1.);

  double wm = 0.;
  if (kr < 100.) {
    wm = std::pow((mean_p / p_max), a) * std::pow((Rm / Rf), b) * (1. / kr);
    if (fp < 1.) {
      wm *= ((Rd - Rm) / (Rf - Rm));
    }
  }

  double wr = 0.;
  if (d < 100.) {
    wr = std::pow(std::max((p_max / p_atm_), 1.), ks) * std::pow((Rm / Rf), d) *
         (rho_d - rho) * std::pow(rho, m);
  }

  // Determine proportions
  const double w_factor = (rho_ratio < 1.) ? std::pow(rho_ratio, 30.) : 1.0;

  // Set w
  const double w = wr * (1. - w_factor) + wm * w_factor;

  return w;
}

//! Convert strain from engineering to true shear strain
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::BoundSurfPlasticity<Tdim>::convert_strain(
    const Vector6d& dstrain) {

  // Note: dstrain (engineering shear strain) is converted to dstrain_t
  //       (true shear strain) for easier tensor multiplication in the Wang
  //       bounding surface plasticity model
  Vector6d dstrain_t = dstrain;
  for (unsigned i = 0; i < 3; ++i) dstrain_t(i + 3) *= 0.5;

  return dstrain_t;
}

//! Set loading bool
template <unsigned Tdim>
bool mpm::BoundSurfPlasticity<Tdim>::set_loading(const Vector6d& dev_dstrain,
                                                 const Vector6d& n_bar) {
  // Solves:
  // (A_ij B_ij) = (A_00*B_00 + A_11*B_11 + A_22*B_22)
  //                + 2*(A_01*B_01 + A_12*B_12 + A_02*B_02)
  // Note that dev_dstrain is engineering shear strain; therefore there is no
  // need to multiply by 2 for shear terms
  double product = 0.;
  for (unsigned i = 0; i < 3; ++i) {
    product +=
        (dev_dstrain(i) * n_bar(i)) + (dev_dstrain(i + 3) * n_bar(i + 3));
  }

  const bool loading = (product > 0.) ? true : false;

  return loading;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::BoundSurfPlasticity<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Relative density dependent terms
  this->compute_relative_density_terms(state_vars);

  // Stress-based initial values
  if ((*state_vars).at("first_loop")) {
    compute_first_loop(stress, state_vars);
    dstrain_last_ = dstrain;
  }

  // REMOVE ME IF TXC //////////////////////////////////////////////////////////
  dstrain_last_ = dstrain;

  // Save stress as muteable vector
  Vector6d sigma = stress;

  // Mean pressure
  double mean_p = check_low(-1. * mpm::materials::p(sigma));

  // Check minimum allowalbe mean pressure
  if (mean_p < p_min_) {
    // Return mapping to minimum allowalbe mean pressure
    for (unsigned i; i < 3; ++i) sigma(i) *= (p_min_ / mean_p);
    mean_p = p_min_;
  }

  // Check maximum allowalbe mean pressure
  if (mean_p > p_max_) {
    // Return mapping to maximum allowalbe mean pressure
    for (unsigned i; i < 3; ++i) sigma(i) *= (p_max_ / mean_p);
    mean_p = p_max_;
  }

  // Deviatoric stress ratio vector and invariant
  Vector6d dev_r = -1. * sigma / mean_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  double R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Check if past maximum allowable R
  const double R_max = density_vars_(9);
  if (R > R_max) {
    // Return mapping parallel to Pi-plane
    sigma *= (R_max / R);
    for (unsigned i = 0; i < 3; ++i) sigma(i) -= (mean_p * (1. - (R_max / R)));

    // Update deviatoric stress ratio and invariant
    dev_r = -1. * sigma / mean_p;
    for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1;
    R = std::sqrt(0.5) * frobenius_norm(dev_r);
  }

  // Dilation surface
  const double Rd = compute_Rd(mean_p, state_vars);

  // Normal to yield surface (incremental stress dependent)
  Vector6d n_bar = compute_n_bar(dev_r, Rd, state_vars);

  // Normal to yield surface (incremental strain dependent)
  const Vector6d n = compute_n(dev_r, dstrain_last_, state_vars);

  // Plastic shear modulus
  const double Hr = compute_Hr(mean_p, state_vars);

  // Deviatoric strain (engineering shear) increment
  const double vol_dstrain = -1. * trace(dstrain_last_);
  Vector6d dev_dstrain = -1. * dstrain_last_;
  for (unsigned i = 0; i < 3; ++i) dev_dstrain(i) -= (vol_dstrain / 3.);

  // dev_dr
  Vector6d dev_dr = Vector6d::Zero();

  // Check loading status
  const bool loading = set_loading(dev_dstrain, n_bar);

  // Set elastic stiffness tensor
  const Matrix6x6 De = this->compute_elastic_tensor(sigma, state_vars);

  // Either loading or unloading
  if (loading) {
    // Loading continues
    // Updated maximum mean pressure
    if (mean_p > (*state_vars).at("p_max")) (*state_vars).at("p_max") = mean_p;

    // Set elasto-plastic stiffness tensor
    const Matrix6x6 Dep = this->compute_elasto_plastic_tensor(
        sigma, dstrain_last_, ptr, state_vars);

    // ADD ME IF TXC ///////////////////////////////////////////////////////////
    // dstrain_last_(0) =
    //     -1. * Dep(0, 2) / (Dep(0, 0) + Dep(0, 1)) * dstrain_last_(2);
    // dstrain_last_(1) =
    //     -1. * Dep(1, 2) / (Dep(1, 0) + Dep(1, 1)) * dstrain_last_(2);

    // Elasto-plastic stress increment
    Vector6d dsigma_ep = Dep * dstrain_last_;

    // Update stress
    sigma += dsigma_ep;

    // Change in mean pressure
    const double mean_dp = (-1. / 3.) * trace(dsigma_ep);
    mean_p = check_low(-1. * mpm::materials::p(sigma));

    // Incremental deviatoric stress ratio vector
    dev_dr = (-dsigma_ep + sigma * mean_dp / mean_p) / mean_p;

  } else {
    // Unloading begins
    // Reset stress reversal point
    Vector6d alpha = -1. * sigma / mean_p;
    for (unsigned i; i < 3; ++i) alpha(i) -= 1.;

    // Save stress reversal point
    (*state_vars).at("alpha0") = alpha(0);
    (*state_vars).at("alpha1") = alpha(1);
    (*state_vars).at("alpha2") = alpha(2);
    (*state_vars).at("alpha3") = alpha(3);
    (*state_vars).at("alpha4") = alpha(4);
    (*state_vars).at("alpha5") = alpha(5);

    // ADD ME IF TXC ///////////////////////////////////////////////////////////
    // dstrain_last_(0) =
    //     -1. * De(0, 2) / (De(0, 0) + De(0, 1)) * dstrain_last_(2);
    // dstrain_last_(1) =
    //     -1. * De(1, 2) / (De(1, 0) + De(1, 1)) * dstrain_last_(2);

    // Elastic stress increment
    const Vector6d dsigma_e = De * dstrain_last_;

    // Update stress
    sigma += dsigma_e;
  }

  // Check minimum allowalbe mean pressure
  double new_p = check_low(-1. * mpm::materials::p(sigma));
  if (new_p < p_min_) {
    // Return mapping to minimum allowalbe mean pressure
    for (unsigned i; i < 3; ++i) sigma(i) *= (p_min_ / new_p);
    new_p = p_min_;
  }

  // Check maximum allowalbe mean pressure
  if (mean_p > p_max_) {
    // Return mapping to maximum allowalbe mean pressure
    for (unsigned i; i < 3; ++i) sigma(i) *= (p_max_ / mean_p);
    mean_p = p_max_;
  }

  // Deviatoric stress ratio vector and invariant
  dev_r = -1. * sigma / new_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Check if past maximum allowable R
  if (R > R_max) {
    // Return mapping parallel to Pi-plane
    sigma *= (R_max / R);
    for (unsigned i = 0; i < 3; ++i) sigma(i) -= (new_p * (1. - (R_max / R)));
  }

  // Compute incremental plastic strain
  Vector6d dsigma = sigma - stress;
  Vector6d dstrain_p = dstrain - (De.inverse() * dsigma);
  if (Tdim == 2) dstrain_p(4) = dstrain_p(5) = 0.;

  // Update equivalent plastic deviatoric strain
  (*state_vars).at("pdstrain") += mpm::materials::pdstrain(dstrain_p);

  // Loading index
  const double L = mean_p * frobenius_prod(dev_dr, n) / Hr;
  const double ddevp = std::sqrt(2.) * std::fabs(L);

  // If Hr too small, ddevp is too big; use 1. as arbitrary cutoff
  double sum_dep = ddevp;
  if (sum_dep > 1.) sum_dep = 0.;

  // Only consider strain if pre-stress surface is past dilation surface
  if ((*state_vars).at("Rm") - Rd > 0.) (*state_vars).at("dep") += sum_dep;

  // Update void ratio
  double dvolumetric_strain = trace(dstrain_last_);
  double void_ratio =
      check_low((*state_vars).at("void_ratio") +
                (1. + void_ratio_initial_) * dvolumetric_strain);
  (*state_vars).at("void_ratio") = void_ratio;

  // Return updated stress
  return (sigma);
}
