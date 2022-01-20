#ifndef MPM_MATERIAL_BOUNDING_SURFACE_PLASTICITY_H_
#define MPM_MATERIAL_BOUNDING_SURFACE_PLASTICITY_H_

#include <cmath>
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

  //! Compute stress invariants (p, q and M_theta)
  //! \param[in] stress Stress
  //! \param[in|out] p Mean stress
  //! \param[in|out] q Deviatoric stress
  void compute_stress_invariants(const Vector6d& stress, double* p, double* q);

  //! Compute state variables (void ratio, p_image, e_image, etc)
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] yield_type Yild type (elastic or yield)
  //! \retval status of computation of stress invariants
  void compute_state_variables(
      const Vector6d& stress, const Vector6d& dstrain,
      mpm::dense_map* state_vars,
      mpm::boundsurfplasticity::FailureState yield_type);

  //! Compute yield function and yield state
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] stress Stress
  //! \retval yield_type Yield type (elastic or yield)
  mpm::boundsurfplasticity::FailureState compute_yield_state(
      double* yield_function, const Vector6d& stress,
      mpm::dense_map* state_vars);

  //! Inline ternary function to check negative or zero numbers
  inline double check_low(double val) {
    return (val > 1.0e-15 ? val : 1.0e-15);
  }

  //! Inline ternary function to check number not greater than one
  inline double check_one(double val) { return (val < 1.0 ? val : 1.0); }

  //! Inline function to scale
  inline double scale(double factor, double Dr) {return (2. - factor + 2. * (factor - 1.) * Dr); }

  // USER INPUTS ///////////////////////////////////////////////////////////////

  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Friction angle
  double friction_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_{std::numeric_limits<double>::max()};
  //! Relative density (decimal)
  double relative_density_{std::numeric_limits<double>::max()};
  //! Initial plastic shear modulus parameter
  double hr0_{std::numeric_limits<double>::max()};
  //! Initial plastic bulk modulus parameter
  double kr0_{std::numeric_limits<double>::max()};
  //! Initial plastic bulk modulus power
  double d0_{std::numeric_limits<double>::max()};
  //! Initial UNKNOWN gamma parameter
  double gamma0_{std::numeric_limits<double>::max()};
  //! Initial UNKNOWN eta (maybe ita) parameter
  double eta0_{std::numeric_limits<double>::max()};
  //! p min
  double pmin_{std::numeric_limits<double>::max()};
  //! Initial maximum stress ratio
  double Rm0_{std::numeric_limits<double>::max()};
  //! Initial stress ratio
  double R0_{std::numeric_limits<double>::max()};
  //! Critical mean pressure
  double pc_{std::numeric_limits<double>::max()};
  //! Reference pressure
  double p_atm_{std::numeric_limits<double>::max()};
  //! Initial mean pressure
  double p0_{std::numeric_limits<double>::max()};
  //! fp ratio
  double fp_{std::numeric_limits<double>::max()};



  // DEFINED BASED ON USER INPUTS //////////////////////////////////////////////

  //! Initial failure surface
  double Rf0_{std::numeric_limits<double>::max()}; 
  //! Failure surface
  double Rf_{std::numeric_limits<double>::max()};
  //! Plastic shear modulus parameter
  double hr_{std::numeric_limits<double>::max()};
  //! Plastic bulk modulus parameter
  double kr_{std::numeric_limits<double>::max()};
  //! Plastic bulk modulus power
  double d_{std::numeric_limits<double>::max()};
  //! UNKNOWN gamma parameter
  double gamma_{std::numeric_limits<double>::max()};
  //! UNKNOWN eta (maybe ita) parameter
  double eta_{std::numeric_limits<double>::max()};
  //! Plastic shear modulus power #1
  double ke_{std::numeric_limits<double>::max()};
  //! Plastic shear modulus power #2
  double ka_{std::numeric_limits<double>::max()};
  //! Plastic bulk modulus power
  double ks_{std::numeric_limits<double>::max()};
  //! R max
  double Rmax_{std::numeric_limits<double>::max()};


  //! Default a parameter
  double a_{std::numeric_limits<double>::max()};
  //! Default b parameter
  double b_{std::numeric_limits<double>::max()};
  
  //! Initial Rp ratio
  double Rp0_{std::numeric_limits<double>::max()};
  //! Rp90 ratio
  double Rp90_{std::numeric_limits<double>::max()};


  // //! Minimum void ratio
  // double lambda_{std::numeric_limits<double>::max()};
  // //! Kappa swelling volumetric
  // double kappa_{std::numeric_limits<double>::max()};
  // //! Gamma void ratio at reference pressure
  // double gamma_{std::numeric_limits<double>::max()};


  // //! Initial void ratio
  // double void_ratio_initial_{std::numeric_limits<double>::max()};

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
