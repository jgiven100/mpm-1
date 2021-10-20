//! Constructor
template <unsigned Tdim>
mpm::DiscontinuityBase<Tdim>::DiscontinuityBase(
    const Json& discontinuity_props) {

  friction_coef_ = 0;

  std::string logger = "discontinuity";
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

  try {
    // assign friction_coef_ if it's given in input file
    if (discontinuity_props.contains("friction_coefficient"))
      friction_coef_ =
          discontinuity_props.at("friction_coefficient").template get<double>();
    // assign cohesion_ if it's given in input file
    if (discontinuity_props.contains("cohesion"))
      cohesion_ = discontinuity_props.at("cohesion").template get<double>();
    // assign contact_distance_ if it's given in input file
    if (discontinuity_props.contains("contact_distance"))
      contact_distance_ =
          discontinuity_props.at("contact_distance").template get<double>();

    // assign move direction if it's given in input file
    if (discontinuity_props.contains("move_direction"))
      move_direction_ =
          discontinuity_props.at("move_direction").template get<int>();
    // assign width if it's given in input file
    if (discontinuity_props.contains("width"))
      width_ = discontinuity_props.at("width").template get<double>();

    // assign maximum_pdstrain if it's given in input file
    if (discontinuity_props.contains("maximum_pdstrain"))
      maximum_pdstrain_ =
          discontinuity_props.at("maximum_pdstrain").template get<double>();

  } catch (Json::exception& except) {
    console_->error("discontinuity parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! create points from file
template <unsigned Tdim>
bool mpm::DiscontinuityBase<Tdim>::create_points(
    const std::vector<VectorDim>& coordinates) {
  bool status = true;
  try {
    // Check if point coordinates is empty
    if (coordinates.empty())
      throw std::runtime_error("List of coordinates is empty");
    // Iterate over all coordinates
    for (const auto& point_coordinates : coordinates) {

      // Add point
      mpm::discontinuity_point<Tdim> point(point_coordinates);

      points_.emplace_back(point);  //
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
//! Locate points in a cell
template <unsigned Tdim>
void mpm::DiscontinuityBase<Tdim>::locate_discontinuity_mesh(
    const Vector<Cell<Tdim>>& cells,
    const Map<Cell<Tdim>>& map_cells) noexcept {
  for (auto& point : this->points_)
    point.locate_discontinuity_mesh(cells, map_cells);
  // for (auto& point : this->points_) point.assign_node_enrich(map_cells);
}

// Compute updated position of the particle
template <unsigned Tdim>
void mpm::DiscontinuityBase<Tdim>::compute_updated_position(
    double dt) noexcept {
  for (auto& point : this->points_)
    point.compute_updated_position(dt, move_direction_);
}

// Compute updated position of the particle
template <unsigned Tdim>
void mpm::DiscontinuityBase<Tdim>::compute_shapefn() noexcept {
  for (auto& point : this->points_) point.compute_shapefn();
}

//! Insert new point
template <unsigned Tdim>
void mpm::DiscontinuityBase<Tdim>::insert_particles(
    VectorDim& coordinates, const Vector<Cell<Tdim>>& cells,
    const Map<Cell<Tdim>>& map_cells) {
  for (auto& point : this->points_) {
    point.assign_tip(false);
  }
  // Add point
  mpm::discontinuity_point<Tdim> point(coordinates);
  point.locate_discontinuity_mesh(cells, map_cells);
  point.compute_shapefn();
  points_.emplace_back(point);
}

//! Output mark points
template <unsigned Tdim>
void mpm::DiscontinuityBase<Tdim>::output_markpoints(int step) {

  std::ofstream path("markpoints.txt", std::ios::app);
  path << step << ":" << std::endl;
  for (auto& point : this->points_) {
    path << point.coordinates()[0] << "    " << point.coordinates()[1] << "    "
         << point.coordinates()[2] << std::endl;
  }
  path << std::endl;
  path.close();
}