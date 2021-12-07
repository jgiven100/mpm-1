//! Initialise xmpm nodal variables
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise_xmpm() noexcept {
  this->initialise();
  // Specific variables for xmpm
  discontinuity_enrich_ = false;
}

//! Initialise shared pointer to nodal properties pool for discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise_discontinuity_property_handle(
    unsigned prop_id,
    std::shared_ptr<mpm::NodalProperties> property_handle) noexcept {
  // the property handle and the property id is set in the node
  this->property_handle_ = property_handle;
  this->discontinuity_prop_id_ = prop_id;
}

//! Update nodal property at the nodes from particle for discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_discontinuity_property(
    bool update, const std::string& property,
    const Eigen::MatrixXd& property_value, unsigned discontinuity_id,
    unsigned nprops) noexcept {
  // Update property
  node_mutex_.lock();
  property_handle_->update_property(property, discontinuity_prop_id_,
                                    discontinuity_id, property_value, nprops);
  node_mutex_.unlock();
}

//! Assign nodal property at the nodes from particle for discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_discontinuity_property(
    bool update, const std::string& property,
    const Eigen::MatrixXd& property_value, unsigned discontinuity_id,
    unsigned nprops) noexcept {
  // assign property
  node_mutex_.lock();
  property_handle_->assign_property(property, discontinuity_prop_id_,
                                    discontinuity_id, property_value, nprops);
  node_mutex_.unlock();
}

// Return data in the nodal properties map at a specific index
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
Eigen::MatrixXd mpm::Node<Tdim, Tdof, Tnphases>::discontinuity_property(
    const std::string& property, unsigned nprops) noexcept {
  // Const pointer to location of property: node_id * nprops x mat_id
  auto property_value =
      property_handle_->property(property, discontinuity_prop_id_, 0, nprops);

  return property_value;
}

//! Compute momentum for discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_momentum_discontinuity(
    unsigned phase, double dt) noexcept {
  momentum_.col(phase) =
      momentum_.col(phase) +
      (internal_force_.col(phase) + external_force_.col(phase)) * dt;
  if (discontinuity_enrich_) {
    property_handle_->update_property(
        "momenta_enrich", discontinuity_prop_id_, 0,
        (property_handle_->property("internal_force_enrich",
                                    discontinuity_prop_id_, 0, Tdim) +
         property_handle_->property("external_force_enrich",
                                    discontinuity_prop_id_, 0, Tdim)) *
            dt,
        Tdim);
  }
  // Apply velocity constraints, which also sets acceleration to 0,
  // when velocity is set.
  this->apply_velocity_constraints_discontinuity();

  this->self_contact_discontinuity(dt);

  this->apply_velocity_constraints_discontinuity();

  return true;
}

//! Compute momentum for discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_momentum_discontinuity_cundall(
    unsigned phase, double dt, double damping_factor) noexcept {
  const double tolerance = 1.0E-15;

  if (!discontinuity_enrich_) {

    if (mass_.col(phase)(0, 0) > tolerance) {

      auto unbalanced_force =
          this->external_force_.col(phase) + this->internal_force_.col(phase);
      this->external_force_.col(phase) -=
          damping_factor * unbalanced_force.norm() *
          this->momentum_.col(phase).normalized();
    }

  } else {

    // obtain the enriched values of enriched nodes
    double mass_enrich = property_handle_->property(
        "mass_enrich", discontinuity_prop_id_, 0, 1)(0, 0);

    // mass for different sides
    double mass_positive = mass_(phase) + mass_enrich;
    double mass_negative = mass_(phase) - mass_enrich;

    Eigen::Matrix<double, Tdim, 1> momenta_enrich = property_handle_->property(
        "momenta_enrich", discontinuity_prop_id_, 0, Tdim);
    Eigen::Matrix<double, Tdim, 1> unbalanced_force =
        this->external_force_.col(phase) + this->internal_force_.col(phase);
    Eigen::Matrix<double, Tdim, 1> unbalanced_force_enrich =
        property_handle_->property("internal_force_enrich",
                                   discontinuity_prop_id_, 0, Tdim) +
        property_handle_->property("external_force_enrich",
                                   discontinuity_prop_id_, 0, Tdim);
    // neet to be fixed
    Eigen::Matrix<double, Tdim, 1> damp_force_positive =
        -damping_factor * (unbalanced_force + unbalanced_force_enrich).norm() *
        (this->momentum_.col(phase) + momenta_enrich).normalized();

    if (mass_positive < tolerance) damp_force_positive.setZero();

    Eigen::Matrix<double, Tdim, 1> damp_force_negative =
        -damping_factor * (unbalanced_force - unbalanced_force_enrich).norm() *
        (this->momentum_.col(phase) - momenta_enrich).normalized();

    if (mass_negative < tolerance) damp_force_negative.setZero();
    this->external_force_.col(phase) +=
        0.5 * (damp_force_positive + damp_force_negative);

    property_handle_->update_property(
        "external_force_enrich", discontinuity_prop_id_, 0,
        0.5 * (damp_force_positive - damp_force_negative), Tdim);
  }

  compute_momentum_discontinuity(phase, dt);

  return true;
}

//! Apply velocity constraints for discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof,
               Tnphases>::apply_velocity_constraints_discontinuity() {
  // Set velocity constraint
  for (const auto& constraint : this->velocity_constraints_) {
    // Direction value in the constraint (0, Dim * Nphases)
    const unsigned dir = constraint.first;
    // Direction: dir % Tdim (modulus)
    const auto direction = static_cast<unsigned>(dir % Tdim);
    // Phase: Integer value of division (dir / Tdim)
    const auto phase = static_cast<unsigned>(dir / Tdim);

    if (!generic_boundary_constraints_) {
      // Velocity constraints are applied on Cartesian boundaries
      this->momentum_(direction, phase) = this->mass(phase) * constraint.second;

      // Set acceleration to 0 in direction of velocity constraint
      this->internal_force_(direction, phase) = 0;
      this->external_force_(direction, phase) = 0;
      if (discontinuity_enrich_) {
        property_handle_->assign_property(
            "momenta_enrich", discontinuity_prop_id_ * Tdim + direction, 0,
            property_handle_->property("mass_enrich", discontinuity_prop_id_, 0,
                                       1) *
                constraint.second,
            1);
        Eigen::Matrix<double, 1, 1> zero_force;
        zero_force.setZero();
        property_handle_->assign_property(
            "internal_force_enrich", discontinuity_prop_id_ * Tdim + direction,
            0, zero_force, 1);
        property_handle_->assign_property(
            "external_force_enrich", discontinuity_prop_id_ * Tdim + direction,
            0, zero_force, 1);
      }
    } else {  // need to be done
    }
  }
}
//! Apply self-contact of the discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::self_contact_discontinuity(
    double dt) noexcept {
  // need to be fixed
  if (!discontinuity_enrich_) return;

  double contact_distance = property_handle_->property(
      "contact_distance", discontinuity_prop_id_, 0, 1)(0, 0);

  Eigen::Matrix<double, Tdim, 1> normal_vector = property_handle_->property(
      "normal_unit_vectors_discontinuity", discontinuity_prop_id_, 0, Tdim);

  if (contact_distance >= 0) return;
  //  single phase for solid
  unsigned phase = 0;
  const double tolerance = 1.0E-15;

  // obtain the enriched values of enriched nodes
  double mass_enrich = property_handle_->property(
      "mass_enrich", discontinuity_prop_id_, 0, 1)(0, 0);
  Eigen::Matrix<double, Tdim, 1> momenta_enrich = property_handle_->property(
      "momenta_enrich", discontinuity_prop_id_, 0, Tdim);

  // mass for different sides
  double mass_positive = mass_(phase) + mass_enrich;
  double mass_negative = mass_(phase) - mass_enrich;

  if (mass_positive < tolerance || mass_negative < tolerance) return;

  // velocity for different sides
  auto velocity_positive =
      (momentum_.col(phase) + momenta_enrich) / mass_positive;
  auto velocity_negative =
      (momentum_.col(phase) - momenta_enrich) / mass_negative;

  // relative normal velocity
  if ((velocity_positive - velocity_negative).col(phase).dot(normal_vector) >=
      0)
    return;

  // the contact momentum, force vector for sticking contact
  auto momentum_contact =
      (mass_enrich * momentum_.col(phase) - mass_(phase) * momenta_enrich) /
      mass_(phase);
  auto force_contact = momentum_contact / dt;

  // friction_coef < 0: move together without slide
  double friction_coef = property_handle_->property(
      "friction_coef", discontinuity_prop_id_, 0, 1)(0, 0);

  if (friction_coef < 0) {
    property_handle_->update_property("momenta_enrich", discontinuity_prop_id_,
                                      0, momentum_contact.col(phase), Tdim);
    property_handle_->update_property("external_force_enrich",
                                      discontinuity_prop_id_, 0,
                                      force_contact.col(phase), Tdim);
  } else {
    // the contact momentum, force value for sticking contact at normal
    // direction
    double momentum_contact_norm =
        momentum_contact.col(phase).dot(normal_vector);
    double force_contact_norm = momentum_contact_norm / dt;

    // the cohesion at nodes
    double cohesion = property_handle_->property(
        "cohesion", discontinuity_prop_id_, 0, 1)(0, 0);

    // the cohesion at nodes
    double cohesion_area = property_handle_->property(
        "cohesion_area", discontinuity_prop_id_, 0, 1)(0, 0);

    double max_friction_force =
        friction_coef * abs(force_contact_norm) + 2 * cohesion * cohesion_area;

    // the contact momentum, force vector for sticking contact at tangential
    // direction
    auto momentum_tangential =
        momentum_contact.col(phase) - momentum_contact_norm * normal_vector;
    auto force_tangential = momentum_tangential / dt;

    // the friction force magnitude
    double force_tangential_value = force_tangential.norm();

    double force_friction = force_tangential_value < max_friction_force
                                ? force_tangential_value
                                : max_friction_force;

    // adjust the momentum and force
    property_handle_->update_property(
        "momenta_enrich", discontinuity_prop_id_, 0,
        momentum_contact_norm * normal_vector +
            force_friction * force_tangential.col(phase).normalized() * dt,
        Tdim);
    property_handle_->update_property(
        "external_force_enrich", discontinuity_prop_id_, 0,
        force_contact_norm * normal_vector +
            force_friction * force_tangential.col(phase).normalized(),
        Tdim);
  }
}

//! Add a cell id
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::add_cell_id(Index id) noexcept {
  cells_.emplace_back(id);
}