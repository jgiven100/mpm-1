#ifndef MPM_MPM_SCHEME_H_
#define MPM_MPM_SCHEME_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mesh.h"

namespace mpm {

//! MPMScheme class
//! \brief Mpmscheme base class to support different stress update schemes
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMScheme {
 public:
  //! Default constructor with mesh class
  MPMScheme(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt);

  //! Intialize
  virtual inline void initialise();

  //! Compute nodal kinematics - map mass and momentum to nodes
  //! \param[in] phase Phase to smooth pressure
  virtual inline void compute_nodal_kinematics(unsigned phase);

  //! Compute stress and strain
  //! \param[in] phase Phase to smooth pressure
  //! \param[in] pressure_smoothing Enable or disable pressure smoothing
  virtual inline void compute_stress_strain(unsigned phase,
                                            bool pressure_smoothing);

  //! Precompute stress and strain (empty call)
  //! \param[in] phase Phase to smooth pressure
  //! \param[in] pressure_smoothing Enable or disable pressure smoothing
  virtual inline void precompute_stress_strain(unsigned phase,
                                               bool pressure_smoothing) = 0;

  //! Postcompute stress and strain (empty call)
  //! \param[in] phase Phase to smooth postssure
  //! \param[in] postssure_smoothing Enable or disable postssure smoothing
  virtual inline void postcompute_stress_strain(unsigned phase,
                                                bool pressure_smoothing) = 0;

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  //! \param[in] pressure_smoothing Enable or disable pressure smoothing
  virtual inline void pressure_smoothing(unsigned phase);

  //! Compute forces
  //! \param[in] gravity Acceleration due to gravity
  //! \param[in] step Number of step in solver
  //! \param[in] concentrated_nodal_forces Boolean for if a concentrated force
  //! is applied or not
  virtual inline void compute_forces(
      const Eigen::Matrix<double, Tdim, 1>& gravity, unsigned phase,
      unsigned step, bool concentrated_nodal_forces);

  //! Compute acceleration velocity position
  //! \param[in] velocity_update Velocity or acceleration update flag
  //! \param[in] phase Phase of particle
  //! \param[in] damping_type Type of damping
  //! \param[in] damping_factor Value of critical damping
  virtual inline void compute_particle_kinematics(
      bool velocity_update, unsigned phase, const std::string& damping_type,
      double damping_factor);

  //! Postcompute nodal kinematics - map mass and momentum to nodes
  //! \param[in] phase Phase to smooth pressure
  virtual inline void postcompute_nodal_kinematics(unsigned phase) = 0;

  //! Compute particle location
  //! \param[in] locate_particles Flag to enable locate particles, if set to
  //! false, unlocated particles will be removed
  virtual inline void locate_particles(bool locate_particles);

  //! Stress update scheme
  //! \retval scheme Stress update scheme
  virtual inline std::string scheme() const = 0;

  /**
   * \defgroup Implicit Functions dealing with implicit MPM
   */
  /**@{*/
  //! Update nodal kinematics by Newmark scheme
  //! \ingroup Implicit
  //! \param[in] newmark_beta Parameter beta of Newmark scheme
  //! \param[in] newmark_gamma Parameter gamma of Newmark scheme
  //! \param[in] phase Phase to smooth pressure
  virtual inline void update_nodal_kinematics_newmark(unsigned phase,
                                                      double newmark_beta,
                                                      double newmark_gamma) {
    throw std::runtime_error(
        "Calling the base class function (update_nodal_kinematics_newmark) in "
        "MPMScheme:: illegal operation!");
  };

  // Update particle stress, strain and volume
  //! \ingroup Implicit
  virtual inline void update_particle_stress_strain_volume() {
    throw std::runtime_error(
        "Calling the base class function "
        "(update_particle_stress_strain_volume) in "
        "MPMScheme:: illegal operation!");
  }
  /**@}*/

  /**
   * \defgroup XMPM Functions dealing with XMPM solver
   */
  /**@{*/
  //! Assign xmpm
  //! \ingroup XMPM
  //! \param[in] xmpm using XMPM solver
  inline void assign_xmpm(bool xmpm) { xmpm_ = xmpm; }
  /**@}*/

 protected:
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Time increment
  double dt_;
  //! MPI Size
  int mpi_size_ = 1;
  //! MPI rank
  int mpi_rank_ = 0;
  //! XMPM slover
  bool xmpm_ = false;
};  // MPMScheme class
}  // namespace mpm

#include "mpm_scheme.tcc"

#endif  // MPM_MPM_SCHEME_H_
