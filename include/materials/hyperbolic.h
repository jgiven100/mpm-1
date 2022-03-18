#ifndef MPM_MATERIAL_HYPERBOLIC_H_
#define MPM_MATERIAL_HYPERBOLIC_H_

#include <cmath>

#include "Eigen/Dense"

#include "infinitesimal_elasto_plastic.h"

namespace mpm {

namespace hyperbolic {
//! Failure state
enum class FailureState { Elastic = 0, Yield = 1 };
}  // namespace hyperbolic

//! Hyperbolic class
//! \brief Hyperbolic Rule material model
//! \details Hyperbolic class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Hyperbolic : public InfinitesimalElastoPlastic<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  //! \param[in] material_properties Material properties
  Hyperbolic(unsigned id, const Json& material_properties);

  //! Destructor
  ~Hyperbolic() override = default;

  //! Delete copy constructor
  Hyperbolic(const Hyperbolic&) = delete;

  //! Delete assignement operator
  Hyperbolic& operator=(const Hyperbolic&) = delete;

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

  //! Compute consistent tangent matrix
  //! \param[in] stress Updated stress
  //! \param[in] prev_stress Stress at the current step
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval dmatrix Constitutive relations mattrix
  Matrix6x6 compute_consistent_tangent_matrix(
      const Vector6d& stress, const Vector6d& prev_stress,
      const Vector6d& dstrain, const ParticleBase<Tdim>* ptr,
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
  Eigen::Matrix<double, 6, 6> compute_elastic_tensor();

  //! Compute constitutive relations matrix for elasto-plastic material
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] hardening Boolean to consider hardening, default=true. If
  //! perfect-plastic tensor is needed pass false
  //! \retval dmatrix Constitutive relations mattrix
  Eigen::Matrix<double, 6, 6> compute_elasto_plastic_tensor(
      const Vector6d& stress, const Vector6d& dstrain,
      const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars,
      bool hardening = true) override;

 private:
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Bulk modulus
  double bulk_modulus_{std::numeric_limits<double>::max()};
  //! Compressional Wave Velocity
  double vp_{std::numeric_limits<double>::max()};
  //! Shear Wave Velocity
  double vs_{std::numeric_limits<double>::max()};
};  // Hyperbolic class
}  // namespace mpm

#include "hyperbolic.tcc"

#endif  // MPM_MATERIAL_HYPERBOLIC_H_
