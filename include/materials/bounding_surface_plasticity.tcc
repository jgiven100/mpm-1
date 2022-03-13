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

    // OPTIONAL: Proportional loading flag
    if (material_properties.find("proportional") != material_properties.end()) {
      proportional_ =
          material_properties.at("proportional").template get<bool>();
    } else {
      proportional_ = true;
    }
    // OPTIONAL: Critical state line void ratio intercept
    if (material_properties.find("void_ratio_intercept") !=
        material_properties.end()) {
      e0_ =
          material_properties.at("void_ratio_intercept").template get<double>();
    } else {
      e0_ = 0.9;
    }
    // OPTIONAL: Critical state line slope
    if (material_properties.find("lambda") != material_properties.end()) {
      lambda_ = material_properties.at("lambda").template get<double>();
    } else {
      lambda_ = 0.03;
    }
    // OPTIONAL: Reference (atomspheric) pressure
    if (material_properties.find("reference_pressure") !=
        material_properties.end()) {
      p_atm_ =
          material_properties.at("reference_pressure").template get<double>();
    } else {
      p_atm_ = 100.e+3;
    }
    // OPTIONAL: Minimum allowable mean pressure
    if (material_properties.find("pmin") != material_properties.end()) {
      p_min_ = material_properties.at("pmin").template get<double>();
    } else {
      p_min_ = 10.;
    }
    // OPTIONAL: Tolerance
    if (material_properties.find("tolerance") != material_properties.end()) {
      tolerance_ = material_properties.at("tolerance").template get<double>();
    }

    // Model parameter: wm function a power
    a_pow_ = 1.;
    // Model parameter: wm function b powers
    b_pow_ = 2.;
    // Model parameter: wr function d power
    d_ = scale(2.0, relative_density_) * d0_;
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

  // Get params
  const double p_atm = this->p_atm_;
  const double G0 = this->G0_;
  const double e0 = this->ec_;
  const double nu = this->poisson_;

  // Mean stress
  const double mean_p = check_low(-1. * mpm::materials::p(stress));

  // Elastic moduli
  const double G = p_atm * G0 * (std::pow((2.973 - e0), 2) / (1. + e0)) *
                   std::sqrt(mean_p / p_atm);
  const double K = (2. * G * (1. + nu)) / (3 * (1 - (2 * nu)));

  // Matrix values
  const double a1 = K + (4.0 / 3.0) * G;
  const double a2 = K - (2.0 / 3.0) * G;

  // Compute elastic stiffness matrix
  Matrix6x6 de = Matrix6x6::Zero();
  // clang-format off
  de(0,0)=a1;    de(0,1)=a2;    de(0,2)=a2;    de(0,3)=0;    de(0,4)=0;    de(0,5)=0;
  de(1,0)=a2;    de(1,1)=a1;    de(1,2)=a2;    de(1,3)=0;    de(1,4)=0;    de(1,5)=0;
  de(2,0)=a2;    de(2,1)=a2;    de(2,2)=a1;    de(2,3)=0;    de(2,4)=0;    de(2,5)=0;
  de(3,0)= 0;    de(3,1)= 0;    de(3,2)= 0;    de(3,3)=G;    de(3,4)=0;    de(3,5)=0;
  de(4,0)= 0;    de(4,1)= 0;    de(4,2)= 0;    de(4,3)=0;    de(4,4)=G;    de(4,5)=0;
  de(5,0)= 0;    de(5,1)= 0;    de(5,2)= 0;    de(5,3)=0;    de(5,4)=0;    de(5,5)=G;
  // clang-format on

  return de;
}

//! Compute elasto-plastic tensor
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6>
    mpm::BoundSurfPlasticity<Tdim>::compute_elasto_plastic_tensor(
        const Vector6d& stress, const Vector6d& dstrain,
        const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars,
        bool hardening) {

  // Get constants
  const double p_atm = this->p_atm_;

  // Get inputs
  const double G0 = this->G0_;
  const double e0 = this->ec_;
  const double nu = this->poisson_;

  // Get power terms
  const double ka = this->ka_;

  // Get distance from stress reversal to current and max pre-stress surfaces
  const double rho = this->rho_;
  const double rho_bar = this->rho_bar_;

  // Get current surfaces
  const double Rm = this->Rm_;

  // Mean stress
  const double mean_p = check_low(-1. * mpm::materials::p(stress));

  // Elastic moduli
  const double G = p_atm * G0 * (std::pow((2.973 - e0), 2) / (1. + e0)) *
                   std::sqrt(mean_p / p_atm);
  const double K = (2. * G * (1. + nu)) / (3 * (1 - (2 * nu)));

  // Compute dilation surface
  const double Rd = compute_Rd(mean_p);

  // Deviatoric stress ratio vector and invariant
  Vector6d dev_r = -1. * stress / mean_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  double R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Normal to yield surface (incremental strain dependent)
  const Vector6d n = compute_n(dev_r, dstrain);

  // Compute fabric dilatancy
  const double Cz = compute_fabric_dilatancy(stress, dstrain, G, K);

  // Plastic shear modulus
  // TODO : think about me!!
  const double m = std::pow(std::min(2 * Rm / rho_bar, 50.), ka);
  const double Hr = compute_Hr(Cz, G, R, m);

  // Plastic bulk modulus
  const double w = compute_w(mean_p, R, m, Rd);
  const double Kr = Cz * K / w;

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
  const Matrix6x6 dep = De - Dp;

  return dep;
}

//! Compute fabric dilatancy term Cz
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_fabric_dilatancy(
    const Vector6d& stress, const Vector6d& dstrain, const double& G,
    const double& K) {

  // Get params
  const double cos_a = this->cos_a_;

  // Compute volumetric plastic strain
  const double dstrain_vp =
      compute_volumetric_plastic_dstrain(stress, dstrain, G, K);

  // Save stress state for incremental calculation in next loop
  this->last_sigma_ = stress;

  // Increment z
  const double sign = cos_a > 0 ? 1. : -1.;
  const double dsvp = dstrain_vp < 0 ? (-1. * dstrain_vp) : 0.;
  const double dz = cz_ * (sign * zm_ - z_) * dsvp;
  this->z_ += dz;

  // Fabric dilatancy
  const double Cz = 1 / (1 + (-1. * sign * z_));

  return Cz;
}

//! Compute necessary paramters for first loop
template <unsigned Tdim>
void mpm::BoundSurfPlasticity<Tdim>::compute_first_loop(
    const Vector6d& stress) {

  // Get params
  const double cos_a = this->cos_a_;
  const double pc = this->pc_;
  const double Rf = this->Rf_;

  // Mean pressure
  double mean_p = check_low(-1. * mpm::materials::p(stress));

  // Deviatoric stress ratio vector and invariant
  Vector6d dev_r = -1. * stress / mean_p;
  for (unsigned i = 0; i < 3; ++i) dev_r(i) -= 1.;
  const double R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Set for later
  this->p_max_ = mean_p;
  this->alpha_ = dev_r;

  // Set maximum pre-stress surface
  this->Rm_ = R;

  // Check minimum allowalbe mean pressure
  if (mean_p < p_min_) mean_p = p_min_;

  // Set distance from reversal point to current stress state
  this->rho_ = std::sqrt(2.) * R;

  // Set distance from reversal point to yield surface
  this->rho_bar_ = std::sqrt(2.) * Rm_;

  // Set distance from reversal point to dilation surface
  const double Rd = compute_Rd(mean_p);
  this->rho_d_ = std::sqrt(2.) * Rd;

  // Set stress state for incremental calculation in next loop
  this->last_sigma_ = stress;

  // Set bool
  this->first_loop_ = false;
}

//! Compute plastic shear modulus coefficient Hr
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_Hr(const double& Cz,
                                                  const double& G,
                                                  const double& R,
                                                  const double& m) {

  // Get coefficient
  const double hr = this->hr_;

  // Get loading type
  const bool proportional = this->proportional_;

  // Get distance from stress reversal to current and max pre-stress surfaces
  const double rho = this->rho_;
  const double rho_bar = this->rho_bar_;

  // Get current surfaces
  const double Rf = this->Rf_;
  const double Rm = this->Rm_;

  // Useful rho terms
  const double rho_ratio = std::min(std::max((rho / rho_bar), 1.E-5), 1.);
  const double rho_ratio_pow = 1. / std::max(std::pow(rho_ratio, m), 1.0E-10);

  // Compute Hr
  double Hr = 0.;
  if (proportional) {
    Hr = Cz * G * hr * (rho_ratio_pow * (Rf / Rm) - 1);
  } else {
    Hr = Cz * G * hr / R;
  }

  return Hr;
}

//! Compute normal to yield surface (dstrain)
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::BoundSurfPlasticity<Tdim>::compute_n(
    const Vector6d& dev_r, const Vector6d& dstrain) {

  // Get params
  const double Rf = this->Rf_;
  const double Rm = this->Rm_;

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
    const Vector6d& dev_r, const double& Rd) {

  // Get params
  const Vector6d alpha = this->alpha_;
  const double rho_bar = this->rho_bar_;
  const double Rm = this->Rm_;

  // Deviatoric stress invariant
  const double R = std::sqrt(0.5) * frobenius_norm(dev_r);

  // Normal to yield (pre-stress) surface
  Vector6d n_bar = Vector6d::Zero();

  // Check if past maximum pre-stress surface
  if (R > 0.999999 * Rm) {
    // Set updated max (pre-stress) surface
    this->Rm_ = R;

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
    this->rho_ = rho;
    this->rho_bar_ = rho_bar;
    this->rho_d_ = rho_d;

    // Normal to yield (pre-stress) surface
    n_bar = r_bar / frobenius_norm(r_bar);
  }

  return n_bar;
}

//! Compute dilation surface
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_Rd(const double& mean_p) {

  // Get params
  const double cos_a = this->cos_a_;
  const double fp = this->fp_;
  const double pc = this->pc_;
  const double Rf = this->Rf_;

  // State index parameter
  const double Ip = mean_p / pc;

  // Dilation surface
  const double Rp90 = 0.9 * Rf;
  const double Rp0 = fp * Rf;  // TODO : determine if correct
  const double Rp = Rp90 + (Rp0 - Rp90) * std::fabs(cos_a);
  const double Rd = Rp + (Rf - Rp) * Ip;

  return Rd;
}

//! Compute volumetric plastic incremental strain
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_volumetric_plastic_dstrain(
    const Vector6d& stress, const Vector6d& dstrain, const double& G,
    const double& K) {
  // SIGN CONVENTION
  // ---------------------------------------------------------------------------
  // NAME       | SIGN          | NOTES
  // ---------------------------------------------------------------------------
  // dstrain    | tension +     | true shear strain
  // dstrain_e  | compression + |
  // dstrain_p  | compression + |
  // dstrain_vp | compression + | (>0) -> contraction; (<0) -> dilation
  // stress     | tension +     | mutable stress vector
  // dsigma     | tension +     |
  // ds         | compression + | incremental deviatoric stress
  // mean_dp    | compression + |
  // ---------------------------------------------------------------------------

  // Incremental stress
  Vector6d dsigma = stress - this->last_sigma_;

  // Incremental mean stress
  const double mean_dp = -1. * mpm::materials::p(dsigma);

  // Incremental deviatoric stress
  Vector6d ds = -1. * dsigma;
  for (unsigned i = 0; i < 3; ++i) ds(i) -= mean_dp;

  // Incremental elastic strain
  Vector6d dstrain_e = (ds / (2. * G));
  for (unsigned i = 0; i < 3; ++i) dstrain_e(i) += (mean_dp / (3. * K));

  // Incremental plastic strain
  Vector6d dstrain_p = (-1. * dstrain) - dstrain_e;

  // Incremental volumetric plastic strain
  const double dstrain_vp = trace(dstrain_p);

  return dstrain_vp;
}

//! Compute plastic bulk modulus coefficient w
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::compute_w(const double& mean_p,
                                                 const double& R,
                                                 const double& m,
                                                 const double& Rd) {

  // Get power terms
  const double a = this->a_pow_;
  const double b = this->b_pow_;
  const double d = this->d_;
  const double ks = this->ks_;

  // Get coefficients
  const double kr = this->kr_;

  // Get stress history terms
  const double p_atm = this->p_atm_;
  const double p_max = this->p_max_;

  // Get loading type
  const bool proportional = this->proportional_;

  // Get distance from stress reversal to current, max pre-stress, and dilation
  // surfaces
  const double rho = this->rho_;
  const double rho_bar = this->rho_bar_;
  const double rho_d = this->rho_d_;

  // Get current surfaces
  const double Rf = this->Rf_;
  const double Rm = this->Rm_;

  // Useful rho ratio
  const double rho_ratio = std::min(std::max((rho / rho_bar), 1.E-5), 1.);

  // Compute w
  double w = 0;
  if (proportional) {

    // Compute wm
    const double wm = (1. / kr) * std::pow((mean_p / p_max), a) *
                      std::pow((Rm / Rf), b) * ((Rd - Rm) / (Rf - Rm));

    // Limit p_max over p_atm ratio
    const double p_max_atm = std::max((p_max / p_atm), 1.);

    // Compute wr
    const double wr = std::pow((Rm / Rf), d) * std::pow(p_max_atm, ks) *
                      (rho_d - rho) * std::pow(rho, m);

    // Combine wm and wr
    const double w_factor = (rho_ratio < 1.) ? std::pow(rho_ratio, 30.) : 1.0;
    w = wr * (1. - w_factor) + wm * w_factor;

  } else {
    w = (1. / kr) * ((Rd - R) / (Rf - Rm));
  }

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

//! Set non-proportional loading value
template <unsigned Tdim>
double mpm::BoundSurfPlasticity<Tdim>::set_loading_dir(const Vector6d& dev_dr,
                                                       const Vector6d& dev_r,
                                                       const bool& loading) {

  // Get params
  const bool proportional = this->proportional_;

  // Loading direction
  double cos_a = -1.;

  // If proportional loading
  if (proportional) {
    if (loading) cos_a_ = 1.;
  } else {
    // Non-proportional loading
    cos_a_ = frobenius_prod(dev_dr, dev_r) / frobenius_norm(dev_dr) /
             frobenius_norm(dev_r);
  }

  return cos_a;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::BoundSurfPlasticity<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // SIGN CONVENTION
  // ---------------------------------------------------------------------------
  // NAME       | SIGN          | NOTES
  // ---------------------------------------------------------------------------
  // dstrain    | tension +     | passed vector, engineering shear strain
  // stress     | tension +     | passed vector
  // sigma      | tension +     | mutable stress vector
  // dsigma_e   | tension +     | elastic portion of stress increase
  // dsigma_ep  | tension +     | elasto-plastic total of stress increase
  // mean_p     | compression + |
  // mean_dp    | compression + |
  // ---------------------------------------------------------------------------

  // Stress-based initial values
  if (first_loop_) compute_first_loop(stress);

  // Save stress and strain as muteable vectors
  Vector6d sigma = stress;

  // Mean pressure
  double mean_p = check_low(-1. * mpm::materials::p(sigma));

  // Check minimum allowalbe mean pressure
  if (mean_p < p_min_) mean_p = p_min_;

  // Compute dilation surface
  const double Rd = compute_Rd(mean_p);

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

  // Normal to yield surface (incremental stress dependent)
  Vector6d n_bar = compute_n_bar(dev_r, Rd);

  // Deviatoric strain (engineering shear) increment
  const double vol_dstrain = -1. * trace(dstrain);
  Vector6d dev_dstrain = -1. * dstrain;
  for (unsigned i = 0; i < 3; ++i) dev_dstrain(i) -= (vol_dstrain / 3.);

  // dev_dr
  Vector6d dev_dr = Vector6d::Zero();

  // Check loading status
  const bool loading = set_loading(dev_dstrain, n_bar);

  // Either loading or unloading
  if (loading) {
    // Loading continues
    // Updated maximum mean pressure
    if (mean_p > p_max_) p_max_ = mean_p;

    // Set elasto-plastic stiffness tensor
    const Matrix6x6 Dep =
        this->compute_elasto_plastic_tensor(sigma, dstrain, ptr, state_vars);

    // Elasto-plastic stress increment
    const Vector6d dsigma_ep = Dep * dstrain;

    // Update stress
    sigma += dsigma_ep;

    // Change in mean pressure
    const double mean_dp = (-1. / 3.) * trace(dsigma_ep);
    mean_p = check_low(-1. * mpm::materials::p(sigma));

    // Incremental deviatoric stress ratio vector
    dev_dr = (-dsigma_ep + sigma * mean_dp / mean_p) / mean_p;

  } else {
    // Unloading begins
    // Reset projection center
    alpha_ = -1. * sigma / mean_p;
    for (unsigned i; i < 3; ++i) alpha_(i) -= 1.;

    // Set elastic stiffness tensor
    const Matrix6x6 De = this->compute_elastic_tensor(sigma, state_vars);

    // Elastic stress increment
    const Vector6d dsigma_e = De * dstrain;

    // Update stress
    sigma += dsigma_e;
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

  // Compute loading direction cos(alpha) for next step
  this->cos_a_ = set_loading_dir(dev_dr, dev_r, loading);

  // Check if past maximum allowable R
  if (R > R_max_) {
    // Update stress
    sigma *= (R_max_ / R);
    for (unsigned i = 0; i < 3; ++i) sigma(i) -= (new_p * (1. - (R_max_ / R)));
  }

  // Return updated stress
  return (sigma);
}
