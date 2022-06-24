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

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Compute elastic tensor
  Eigen::Matrix<double, 6, 6> compute_tensor(const Vector6d& strain);

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
  //! Beta arameter
  double beta_{std::numeric_limits<double>::max()};
  //! s parameter
  double s_{std::numeric_limits<double>::max()};
  //! Reference strain
  double reference_strain_{std::numeric_limits<double>::max()};
  //! p1 parameter
  double p1_{std::numeric_limits<double>::max()};
  //! p2 parameter
  double p2_{std::numeric_limits<double>::max()};
  //! p3 paramaeter
  double p3_{std::numeric_limits<double>::max()};
  //! Calculated Stress
  Eigen::Matrix<double, 6, 1> calc_stress_;
  //! Reveral Strain
  Eigen::Matrix<double, 6, 1> reversal_strain_;
  //! Reversal Stress
  Eigen::Matrix<double, 6, 1> reversal_stress_;
  //! Maximum Stress
  Eigen::Matrix<double, 6, 1> maximum_strain_;
  //! Initial Strain Direction
  Eigen::Matrix<double, 6, 1> init_dir_;
  //! New Strain Direction
  Eigen::Matrix<double, 6, 1> new_dir_;
  //! Secant modulus
  double secant_modulus_{std::numeric_limits<double>::max()};
  //! Initial loading
  bool init_loading_{true};

};  // Hyperbolic class
}  // namespace mpm

#include "hyperbolic.tcc"

#endif  // MPM_MATERIAL_HYPERBOLIC_H_
