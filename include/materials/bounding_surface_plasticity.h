#ifndef MPM_MATERIAL_BOUNDING_SURFACE_PLASTICITY_H_
#define MPM_MATERIAL_BOUNDING_SURFACE_PLASTICITY_H_

#include <cmath>
#include <iomanip>
#include <iostream>

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
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Compute elastic tensor
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  Eigen::Matrix<double, 6, 6> compute_elastic_tensor(
      const Vector6d& stress, mpm::dense_map* state_vars);

  //! Compute constitutive relations matrix for elasto-plastic material
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] hardening Boolean to consider hardening, default=true. If
  //! perfect-plastic tensor is needed pass false
  //! \retval dmatrix Constitutive relations mattrix
  Matrix6x6 compute_elasto_plastic_tensor(const Vector6d& stress,
                                          const Vector6d& dstrain,
                                          const ParticleBase<Tdim>* ptr,
                                          mpm::dense_map* state_vars,
                                          bool hardening = true) override;

  //! Inline ternary function to check negative or zero numbers
  inline double check_low(double val) {
    return (val > 1.0e-15 ? val : 1.0e-15);
  }

  //! Inline ternary function to check number not greater than one
  inline double check_one(double val) { return (val < 1.0 ? val : 1.0); }

  //! Compute Frobenius norm
  //! \param[in] vec_A vector form of 2D matrix (i.e., Voigt notation)
  //! \retval Frobenius norm of 1 matrix
  double frobenius_norm(const Vector6d& vec_A);

  //! Compute Frobenius inner product
  //! \param[in] vec_A vector form of 2D matrix (i.e., Voigt notation)
  //! \param[in] vec_B vector form of 2D matrix (i.e., Voigt notation)
  //! \retval Frobenius innner product of 2 matrices
  double frobenius_prod(const Vector6d& vec_A, const Vector6d& vec_B);

  //! Inline function to scale
  inline double scale(double factor, double Dr) {
    return (2. - factor + 2. * (factor - 1.) * Dr);
  }

  //! Model parameter: wm function a power
  double a_pow_{std::numeric_limits<double>::max()};
  //! Previous stress reversal point
  double alpha_s_{std::numeric_limits<double>::max()};
  //! Previous stress reversal point
  Eigen::Matrix<double, 6, 1> alpha_v_{Vector6d::Zero()};
  //! Model parameter: wm function b power
  double b_pow_{std::numeric_limits<double>::max()};
  //! Input parameter: wr function d power
  double d0_{std::numeric_limits<double>::max()};
  //! Model parameter: wr function d power
  double d_{std::numeric_limits<double>::max()};
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Cumulative plastic deviatoric strain
  double dep_{std::numeric_limits<double>::max()};
  //! Critical void ratio
  double ec_{std::numeric_limits<double>::max()};
  //! CSL void ratio intercept
  double e0_{std::numeric_limits<double>::max()};
  //! Input parameter: Ci function eta coefficient
  double eta0_{std::numeric_limits<double>::max()};
  //! Model parameter: Ci function eta coefficient
  double eta_{std::numeric_limits<double>::max()};
  //! Initialization bool
  bool first_loop_{true};
  //! Ratio fp = Rp / Rf
  double fp_{std::numeric_limits<double>::max()};
  //! Friction angle
  double friction_{std::numeric_limits<double>::max()};
  //! Shear modulus model parameter
  double G0_{std::numeric_limits<double>::max()};
  //! Input parameter: Ci and Cg functions gm coefficient
  double gm0_{std::numeric_limits<double>::max()};
  //! Model parameter: Ci and Cg functions gm coefficient
  double gm_{std::numeric_limits<double>::max()};
  //! Input parameter: Hr function hr coefficient
  double hr0_{std::numeric_limits<double>::max()};
  //! Model parameter: Hr function hr coefficient
  double hr_{std::numeric_limits<double>::max()};
  //! Model parameter: m function ka power
  double ka_{std::numeric_limits<double>::max()};
  //! Model parameter: Ci and Cg function ke power
  double ke_{std::numeric_limits<double>::max()};
  //! Input parameter: wm function kr coefficient
  double kr0_{std::numeric_limits<double>::max()};
  //! Model parameter: wm function kr coefficient
  double kr_{std::numeric_limits<double>::max()};
  //! Model parameter: wr function ks power
  double ks_{std::numeric_limits<double>::max()};
  //! Critial state slope
  double lambda_{std::numeric_limits<double>::max()};
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
  //! Maximum allowable R
  double R_max_{std::numeric_limits<double>::max()};
  //! Dilation surface
  double Rd_{std::numeric_limits<double>::max()};
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
  //! Failure state map
  std::map<int, mpm::boundsurfplasticity::FailureState> yield_type_ = {
      {0, mpm::boundsurfplasticity::FailureState::Elastic},
      {1, mpm::boundsurfplasticity::FailureState::Yield}};

};  // Bounding Surface Plasticity class
}  // namespace mpm

#include "bounding_surface_plasticity.tcc"

#endif  // MPM_MATERIAL_BOUNDING_SURFACE_PLASTICITY_H_
