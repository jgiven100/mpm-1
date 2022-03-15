#ifndef MPM_MATERIAL_BOUNDING_SURFACE_PLASTICITY_H_
#define MPM_MATERIAL_BOUNDING_SURFACE_PLASTICITY_H_

#include <cmath>
#include <iomanip>

#include "Eigen/Dense"

#include "infinitesimal_elasto_plastic.h"

namespace mpm {

namespace boundsurfplasticity {
//! Failure state
enum class FailureState { Elastic = 0, Yield = 1 };
}  // namespace boundsurfplasticity

//! BoundSurfPlasticity class
//! \brief Bounding Surface Plasticity material model
//! \details Bounding Surface Plasticity material model
//! \tparam Tdim Dimension
template <unsigned Tdim>
class BoundSurfPlasticity : public InfinitesimalElastoPlastic<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id and material properties
  //! \param[in] material_properties Material properties
  BoundSurfPlasticity(unsigned id, const Json& material_properties);

  //! Destructor
  ~BoundSurfPlasticity() override = default;

  //! Delete copy constructor
  BoundSurfPlasticity(const BoundSurfPlasticity&) = delete;

  //! Delete assignement operator
  BoundSurfPlasticity& operator=(const BoundSurfPlasticity&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! State variables
  std::vector<std::string> state_variables() const override;

  //! Initialise material
  //! \brief Function that initialise material to be called at the beginning of
  //! time step
  void initialise(mpm::dense_map* state_vars) override {
    (*state_vars).at("yield_state") = 0;
  };

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval sigma Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

  //! Compute strain
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_strain(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars);

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Inline ternary function to check negative or zero numbers
  inline double check_low(double val) {
    return (val > 1.0e-15 ? val : 1.0e-15);
  }

  //! Compute elastic tensor
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval de Constitutive relations matrix
  Eigen::Matrix<double, 6, 6> compute_elastic_tensor(
      const Vector6d& stress, mpm::dense_map* state_vars);

  //! Compute constitutive relations matrix for elasto-plastic material
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] hardening Boolean to consider hardening, default=true. If
  //! perfect-plastic tensor is needed pass false
  //! \retval Dep Constitutive relations matrix
  Matrix6x6 compute_elasto_plastic_tensor(const Vector6d& stress,
                                          const Vector6d& dstrain,
                                          const ParticleBase<Tdim>* ptr,
                                          mpm::dense_map* state_vars,
                                          bool hardening = true) override;

  //! Compute fabric dilatancy term Cz
  //! \param[in] stress Stress
  //! \param[in] dstrain Incremental strain
  //! \param[in] G Shear modulus
  //! \param[in] K Bulk modulus
  //! \retval Cz Fabric dilatancy term
  double compute_fabric_dilatancy(const Vector6d& stress,
                                  const Vector6d& dstrain, const double& G,
                                  const double& K);

  //! Compute first loop (stress-based) functions
  //! \param[in] stress Stress
  void compute_first_loop(const Vector6d& stress);

  //! Compute plastic shear modulus coefficient Hr
  //! \param[in] Cz fabric dilatancy term
  //! \param[in] G Shear modulus
  //! \param[in] R Deviatoric stress ratio invariant
  //! \param[in] m Plastic shear modulus power
  //! \retval Hr Plastic shear modulus coefficient
  double compute_Hr(const double& Cz, const double& G, const double& R,
                    const double& m);

  //! Compute unit normal to yield surface (incremental strain)
  //! \param[in] dev_r Deviatoric stress ratio vector
  //! \param[in] dstrain Strain increment vector
  //! \retval n Unit normal to yield surface (incremental strain)
  Vector6d compute_n(const Vector6d& dev_r, const Vector6d& dstrain);

  //! Compute unit normal to yield surface (incremental stress)
  //! \param[in] dev_r Deviatoric stress ratio vector
  //! \param[in] Rd Dilation surface
  //! \retval n_bar Unit normal to yield surface (incremental stress)
  Vector6d compute_n_bar(const Vector6d& dev_r, const double& Rd);

  //! Compute dilation surface
  //! \param[in] mean_p Current mean stress
  //! \retval Rd Deviatoric stress invariant defining the dilation surface
  double compute_Rd(const double& mean_p);

  //! Compute volumetric plastic incremental strain
  //! \param[in] stress Stress
  //! \param[in] dstrain Incremental strain
  //! \param[in] G Shear modulus
  //! \param[in] K Bulk modulus
  //! \retval dstrain_vp Volumetric plastic incremental strain
  double compute_volumetric_plastic_dstrain(const Vector6d& stress,
                                            const Vector6d& dstrain,
                                            const double& G, const double& K);

  //! Compute plastic bulk modulus coefficient w
  //! \param[in] mean_p Current mean stress
  //! \param[in] R Current deviatoric stress invariant
  //! \param[in] m Plastic modulus power
  //! \param[in] Rd Stress invariant for current dilation surface
  //! \retval w Plastic bulk modulus coefficient
  double compute_w(const double& mean_p, const double& R, const double& m,
                   const double& Rd);

  //!  Convert strain from engineering to true shear strain
  //! \param[in] dstrain Incremental strain
  //! \retval dstrain_t incremental strain with true shear strain values
  Vector6d convert_strain(const Vector6d& dstrain);

  //! Compute Frobenius norm
  //! \param[in] vec_A vector form of 2D matrix (i.e., Voigt notation)
  //! \retval norm Frobenius norm of 1 matrix
  double frobenius_norm(const Vector6d& vec_A);

  //! Compute Frobenius inner product
  //! \param[in] vec_A vector form of 2D matrix (i.e., Voigt notation)
  //! \param[in] vec_B vector form of 2D matrix (i.e., Voigt notation)
  //! \retval product Frobenius innner product of 2 matrices
  double frobenius_prod(const Vector6d& vec_A, const Vector6d& vec_B);

  //! Set loading bool
  //! \param[in] dev_dstrain Deviatoric strain (engineering shear) increment
  //! \param[in] n_bar Normal to yield surface (incremental stress dependent)
  //! \retval loading Bool for loading status
  bool set_loading(const Vector6d& dev_dstrain, const Vector6d& n_bar);

  //! Set loading direction
  //! \param[in] dev_dr Deviatoric stress ratio vector increment
  //! \param[in] dev_r Deviatoric stress ratio vector
  //! \param[in] loading Bool for loading or unloading
  //! \retval cos_a Loading direction value
  double set_loading_dir(const Vector6d& dev_dr, const Vector6d& dev_r,
                         const bool& loading);

  //! Inline function to scale
  inline double scale(const double& factor, const double& Dr) {
    return (2. - factor + 2. * (factor - 1.) * Dr);
  }

  //! Inline function to get trace of 3x3 tensor (Voigt notation)
  inline double trace(const Vector6d& vec) {
    return (vec(0) + vec(1) + vec(2));
  }

  //! Model parameter: wm function a power
  double a_pow_{std::numeric_limits<double>::max()};
  //! Previous stress reversal point
  Eigen::Matrix<double, 6, 1> alpha_{Vector6d::Zero()};
  //! Model parameter: wm function b power
  double b_pow_{std::numeric_limits<double>::max()};
  //! Loading direction representation cos(alpha)
  double cos_a_{std::numeric_limits<double>::max()};
  //! Input parameter: dz function cz coefficient
  double cz_{std::numeric_limits<double>::max()};
  //! Input parameter: wr function d power
  double d0_{std::numeric_limits<double>::max()};
  //! Model parameter: wr function d power
  double d_{std::numeric_limits<double>::max()};
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! CSL void ratio intercept
  double e0_{std::numeric_limits<double>::max()};
  //! Critical void ratio
  double ec_{std::numeric_limits<double>::max()};
  //! Initialization bool
  bool first_loop_{true};
  //! Ratio fp = Rp / Rf
  double fp_{std::numeric_limits<double>::max()};
  //! Friction angle
  double friction_{std::numeric_limits<double>::max()};
  //! Shear modulus model parameter
  double G0_{std::numeric_limits<double>::max()};
  //! Input parameter: Hr function hr coefficient
  double hr0_{std::numeric_limits<double>::max()};
  //! Model parameter: Hr function hr coefficient
  double hr_{std::numeric_limits<double>::max()};
  //! Model parameter: m function ka power
  double ka_{std::numeric_limits<double>::max()};
  //! Input parameter: wm function kr coefficient
  double kr0_{std::numeric_limits<double>::max()};
  //! Model parameter: wm function kr coefficient
  double kr_{std::numeric_limits<double>::max()};
  //! Model parameter: wr function ks power
  double ks_{std::numeric_limits<double>::max()};
  //! Critial state slope
  double lambda_{std::numeric_limits<double>::max()};
  //! Last stress state
  Eigen::Matrix<double, 6, 1> last_sigma_{Vector6d::Zero()};
  //! Reference (atmospheric) pressure
  double p_atm_{std::numeric_limits<double>::max()};
  //! Max mean pressure
  double p_max_{std::numeric_limits<double>::max()};
  //! Minimum allowable mean pressure
  double p_min_{std::numeric_limits<double>::max()};
  //! Critical state mean pressure
  double pc_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_{std::numeric_limits<double>::max()};
  //! Proportional loading input flag for cos(alpha) calculations
  bool proportional_{true};
  //! Maximum allowable R
  double R_max_{std::numeric_limits<double>::max()};
  //! Input parameter: failure surface
  double Rf0_{std::numeric_limits<double>::max()};
  //! Failure surface
  double Rf_{std::numeric_limits<double>::max()};
  //! Maximum pre-stress surface
  double Rm_{std::numeric_limits<double>::max()};
  //! Relative density (decimal)
  double relative_density_{std::numeric_limits<double>::max()};
  //! Distance from reversal point to current state
  double rho_{std::numeric_limits<double>::max()};
  //! Distance from reversal point to yield surface
  double rho_bar_{std::numeric_limits<double>::max()};
  //! Distance from reversal point to dilation surface
  double rho_d_{std::numeric_limits<double>::max()};
  //! Default tolerance
  double tolerance_{std::numeric_limits<double>::epsilon()};
  //! Model parameter: Cz function z coefficient
  double z_{std::numeric_limits<double>::max()};
  //! Input parameter: dz function zm coefficient
  double zm_{std::numeric_limits<double>::max()};
  //! Failure state map
  std::map<int, mpm::boundsurfplasticity::FailureState> yield_type_ = {
      {0, mpm::boundsurfplasticity::FailureState::Elastic},
      {1, mpm::boundsurfplasticity::FailureState::Yield}};

};  // Bounding Surface Plasticity class
}  // namespace mpm

#include "bounding_surface_plasticity.tcc"

#endif  // MPM_MATERIAL_BOUNDING_SURFACE_PLASTICITY_H_
