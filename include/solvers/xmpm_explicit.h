#ifndef MPM_XMPM_EXPLICIT_H_
#define MPM_XMPM_EXPLICIT_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "discontinuity_base.h"
#include "mpm_base.h"

namespace mpm {

//! XMPMExplicit class
//! \brief A class that implements the fully explicit one phase xmpm
//! \details A single-phase explicit XMPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class XMPMExplicit : public MPMBase<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;
  //! Default constructor
  XMPMExplicit(const std::shared_ptr<IO>& io);

  //! Solve
  bool solve() override;

  //! Compute stress strain
  //! \param[in] phase Phase to smooth pressure
  void compute_stress_strain(unsigned phase);

  //! Initialise discontinuity
  void initialise_discontinuity();

 protected:
  // Generate a unique id for the analysis
  using mpm::MPMBase<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMBase<Tdim>::dt_;
  //! Current step
  using mpm::MPMBase<Tdim>::step_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nsteps_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nload_balance_steps_;
  //! Output steps
  using mpm::MPMBase<Tdim>::output_steps_;
  //! A unique ptr to IO object
  using mpm::MPMBase<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPMBase<Tdim>::analysis_;
  //! JSON post-process object
  using mpm::MPMBase<Tdim>::post_process_;
  //! Logger
  using mpm::MPMBase<Tdim>::console_;
  //! MPM Scheme
  using mpm::MPMBase<Tdim>::mpm_scheme_;
  //! Stress update method
  using mpm::MPMBase<Tdim>::stress_update_;
  //! Interface scheme
  using mpm::MPMBase<Tdim>::contact_;
  //! xmpm solver
  using mpm::MPMBase<Tdim>::xmpm_;

#ifdef USE_GRAPH_PARTITIONING
  //! Graph
  using mpm::MPMBase<Tdim>::graph_;
#endif

  //! velocity update
  using mpm::MPMBase<Tdim>::velocity_update_;
  //! Gravity
  using mpm::MPMBase<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPMBase<Tdim>::mesh_;
  //! Materials
  using mpm::MPMBase<Tdim>::materials_;
  //! Node concentrated force
  using mpm::MPMBase<Tdim>::set_node_concentrated_force_;
  //! Damping type
  using mpm::MPMBase<Tdim>::damping_type_;
  //! Damping factor
  using mpm::MPMBase<Tdim>::damping_factor_;
  //! Locate particles
  using mpm::MPMBase<Tdim>::locate_particles_;

 private:
  //! Pressure smoothing
  bool pressure_smoothing_{false};
  //! Interface
  bool interface_{false};
  //! With or without discontinuities
  bool setdiscontinuity_{false};
  //! Discontinuities
  std::shared_ptr<mpm::DiscontinuityBase<Tdim>> discontinuity_;
  //! Describe a discontinuity by mesh
  bool surfacemesh_{false};
  //! Initialize the discontinuity by level set values
  bool particle_levelet_{false};
  //! Proparate or not
  bool propagation_{false};
  //! Initiate or not
  bool initiation_{false};
  //! Compute the nodal level set values by: "shepard" of "mls"
  std::string nodal_levelset_{"shepard"};
  //! Compute the average friction coefficient from the neighbour particles
  bool friction_coef_average_{false};

};  // XMPMExplicit class
}  // namespace mpm

#include "xmpm_explicit.tcc"

#endif  // MPM_XMPM_EXPLICIT_H_
