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
    beta_ = material_properties.at("beta").template get <double>();
    s_ = material_properties.at("s").template get <double>();
    reference_strain_ = material_properties.at("reference_stress").template get <double>();
    p1_ = material_properties.at("p1").template get <double>();
    p2_ = material_properties.at("p2").template get <double>();
    p3_ = material_properties.at("p3").template get <double>();

    // Set elastic tensor
    this->compute_elastic_tensor();
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
  // Shear modulus
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

  Matrix6x6 de = Matrix6x6::Zero();
  // clang-format off
  // compute elasticityTensor
  de(0,0)=a1;    de(0,1)=a2;    de(0,2)=a2;    de(0,3)=0;    de(0,4)=0;    de(0,5)=0;
  de(1,0)=a2;    de(1,1)=a1;    de(1,2)=a2;    de(1,3)=0;    de(1,4)=0;    de(1,5)=0;
  de(2,0)=a2;    de(2,1)=a2;    de(2,2)=a1;    de(2,3)=0;    de(2,4)=0;    de(2,5)=0;
  de(3,0)= 0;    de(3,1)= 0;    de(3,2)= 0;    de(3,3)=G;    de(3,4)=0;    de(3,5)=0;
  de(4,0)= 0;    de(4,1)= 0;    de(4,2)= 0;    de(4,3)=0;    de(4,4)=G;    de(4,5)=0;
  de(5,0)= 0;    de(5,1)= 0;    de(5,2)= 0;    de(5,3)=0;    de(5,4)=0;    de(5,5)=G;
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

    if ((strain.cwiseEqual(0)).all()) {
        // Determine directionality of initial strain
        Vector6d init_dir_ = dstrain.cwiseSign();
        bool init_loading = true;
    } else {
        Vector6d new_dir_ = dstrain.cwiseSign();
        if ((new_dir_.array() != init_dir_.array()).any()) {
            bool init_loading_ = false;
            Vector6d reversal_strain_ = strain;
            Vector6d reversal_stress_ = stress;
            if ((maximum_strain_.cwiseEqual(0)).all()) {
                Vector6d maximum_strain_ = strain;
                double secant_modulus_ = stress(0)/stress(0);
            } else {
                if ((strain.array().abs() > maximum_strain_.array().abs()).any()) {
                    Vector6d maximum_strain_ = strain;
                    double secant_modulus_ = stress(0)/stress(0);
                }
            }
        }
    }

    double shear_modulus = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));
    if (init_loading_) {
        Vector6d calc_stress_ = (shear_modulus*strain.array()*(1+beta_*pow(strain.array() / reference_strain_,s_).inverse())).matrix();
    } else {
        // Strain Vectors and subtracting reversal strain
        Vector6d strain_dif = strain - reversal_strain_;

        // Calculating Reduction Factor
        double reduction_factor = (p1_-p2_*pow(1.-secant_modulus_/shear_modulus,p3_));

        // Calculating Second Term
        Vector6d term2 = (shear_modulus*strain_dif.array()*(1.+beta_*pow(reference_strain_*strain_dif.array()/2,s_)).inverse()).matrix();

        // Calculating Third Term
        Vector6d term3 = (shear_modulus*strain_dif.array()*(1.+beta_*pow(maximum_strain_.array()-reference_strain_,s_)).inverse()).matrix();

        Vector6d calc_stress_ = reduction_factor*(term2+term3)+term3+reversal_stress_;
    }
    return calc_stress_;
}


//! Compute consistent tangent matrix
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6>
    mpm::Hyperbolic<Tdim>::compute_consistent_tangent_matrix(
        const Vector6d& stress, const Vector6d& prev_stress,
        const Vector6d& dstrain, const ParticleBase<Tdim>* ptr,
        mpm::dense_map* state_vars) {
  const Matrix6x6 de = this->compute_elastic_tensor();
  return de;
}