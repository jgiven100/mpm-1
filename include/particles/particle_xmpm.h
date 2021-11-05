#ifndef MPM_PARTICLE_XMPM_H_
#define MPM_PARTICLE_XMPM_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "cell.h"
#include "logger.h"
#include "particle_base.h"

namespace mpm {

//! ParticleXMPM class
//! \brief ParticleXMPM class derived from particle class, stores the
//! information for XMPMExplicit solver \details ParticleXMPM class: id_ and
//! coordinates. \tparam Tdim Dimension
template <unsigned Tdim>
class ParticleXMPM : public Particle<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  ParticleXMPM(Index id, const VectorDim& coord);

  //! Destructor
  ~ParticleXMPM() override = default;

  //! Delete copy constructor
  ParticleXMPM(const ParticleXMPM<Tdim>&) = delete;

  //! Delete assignment operator
  ParticleXMPM& operator=(const ParticleXMPM<Tdim>&) = delete;

  //! Initialise particle from POD data
  //! \param[in] particle POD data of particle
  //! \retval status Status of reading HDF5 particle
  bool initialise_particle(PODParticle& particle) override;

  //! Return particle data as POD
  //! \retval particle POD of the particle
  std::shared_ptr<void> pod() const override;

  //! Initialise properties
  void initialise() override;
  //! Type of particle
  std::string type() const override {
    return (Tdim == 2) ? "P2DXMPM" : "P3DXMPM";
  }
  //! Map particle mass and momentum to nodes
  void map_mass_momentum_to_nodes() noexcept override;

  //! Map body force
  //! \param[in] pgravity Gravity of a particle
  void map_body_force(const VectorDim& pgravity) noexcept override;

  //! Map traction force
  void map_traction_force() noexcept override;

  //! Map internal force
  inline void map_internal_force() noexcept override;

  //! Map particle levelset to nodes
  void map_levelset_to_nodes() noexcept override;

  //! Map particle frictional_coef to nodes
  //! \param[in] friction_coef of the discontinuity
  void map_friction_coef_to_nodes(
      double discontinuity_friction_coef) noexcept override;

  //! Compute displacement gradient
  //! \param[in] dt Analysis time step
  void inline compute_displacement_gradient(double dt) noexcept override;

  //! to do
  //!   virtual void check_levelset() noexcept override;

  //! Compute updated position of the particle
  //! \param[in] dt Analysis time step
  //! \param[in] velocity_update Update particle velocity from nodal vel
  void compute_updated_position(double dt,
                                bool velocity_update = false) noexcept override;

  //! Compute strain rate
  //! \param[in] dn_dx The spatial gradient of shape function
  //! \param[in] phase Index to indicate phase
  //! \retval strain rate at particle inside a cell
  inline Eigen::Matrix<double, 6, 1> compute_strain_rate(
      const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept;

  //! Map levelset from nodes to particles
  void map_levelset_to_particle();

  //! return levelset values
  //! \retval particle levelset values
  double levelset_phi() { return levelset_phi_; }

  //! compute the minimum eigenvalue of the acoustic tensor
  //! \param[in] the normal direction of the previous discontinuity
  //! \param[in] do the initiation detection loop
  bool minimum_acoustic_tensor(VectorDim& normal_cell, bool initiation);

  //! compute the initiation normal direction
  //! \param[in] initiation normal direction
  void compute_initiation_normal(VectorDim& normal);

 private:
  //! Assign the level set function values
  //! \param[in] phivalue The level set values
  void assign_levelsetphi(double phivalue) { levelset_phi_ = phivalue; };

  //! Return 1 if x > 0, -1 if x < 0 and 0 if x = 0
  //! \param[in] x double value
  inline double sgn(double x) noexcept {
    return (x > 0) ? 1. : ((x < 0) ? -1. : 0);
  }

 private:
  //! particle id
  using Particle<Tdim>::id_;
  //! coordinates
  using Particle<Tdim>::coordinates_;
  //! Reference coordinates (in a cell)
  using Particle<Tdim>::xi_;
  //! Cell
  using Particle<Tdim>::cell_;
  //! Cell id
  using Particle<Tdim>::cell_id_;
  //! Nodes
  using Particle<Tdim>::nodes_;
  //! Status
  using Particle<Tdim>::status_;
  //! Material
  using Particle<Tdim>::material_;
  //! Material id
  using Particle<Tdim>::material_id_;
  //! State variables
  using Particle<Tdim>::state_variables_;
  //! Neighbour particles
  using Particle<Tdim>::neighbours_;
  //! Volumetric mass density (mass / volume)
  using Particle<Tdim>::mass_density_;
  //! Mass
  using Particle<Tdim>::mass_;
  //! Volume
  using Particle<Tdim>::volume_;
  //! Size of particle
  using Particle<Tdim>::size_;
  //! Size of particle in natural coordinates
  using Particle<Tdim>::natural_size_;
  //! Stresses
  using Particle<Tdim>::stress_;
  //! Strains
  using Particle<Tdim>::strain_;
  //! dvolumetric strain
  using Particle<Tdim>::dvolumetric_strain_;
  //! Volumetric strain at centroid
  using Particle<Tdim>::volumetric_strain_centroid_;
  //! Strain rate
  using Particle<Tdim>::strain_rate_;
  //! dstrains
  using Particle<Tdim>::dstrain_;
  //! Velocity
  using Particle<Tdim>::velocity_;
  //! Displacement
  using Particle<Tdim>::displacement_;
  //! Particle velocity constraints
  using Particle<Tdim>::particle_velocity_constraints_;
  //! Set traction
  using Particle<Tdim>::set_traction_;
  //! Surface Traction (given as a stress; force/area)
  using Particle<Tdim>::traction_;
  //! Shape functions
  using Particle<Tdim>::shapefn_;
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! dN/dX at cell centroid
  using Particle<Tdim>::dn_dx_centroid_;
  //! Logger
  using Particle<Tdim>::console_;
  //! Map of scalar properties
  using Particle<Tdim>::scalar_properties_;
  //! Map of vector properties
  using Particle<Tdim>::vector_properties_;
  //! Map of tensor properties
  using Particle<Tdim>::tensor_properties_;
  //! Pack size
  using Particle<Tdim>::pack_size_;
  //! du/dx at particles
  Eigen::Matrix<double, 3, 3> du_dx_;

 private:
  //! level set value： phi for discontinuity
  double levelset_phi_{0.};

  //! the minimum eigenvalue of the acoustic tensor
  double minimum_acoustic_eigenvalue_{1e16};

  //! angle
  double discontinuity_angle_{0.};
};  // ParticleXMPM class
}  // namespace mpm

#include "particle_xmpm.tcc"

#endif  // MPM_PARTICLE_XMPM_H__
