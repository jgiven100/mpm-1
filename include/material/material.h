#ifndef MPM_MATERIAL_H_
#define MPM_MATERIAL_H_

#include <limits>

#include "Eigen/Dense"

namespace mpm {

// Material class
//! \brief Base class that stores the information about materials
//! \details Material class stresses and strains
class Material {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  // Constructor with id
  Material(unsigned id);

  //! Destructor
  virtual ~Material(){};

  //! Delete copy constructor
  Material(const Material&) = delete;

  //! Delete assignement operator
  Material& operator=(const Material&) = delete;

  //! Return id of the material
  unsigned id() const { return id_; }

 protected:
  //! material id
  unsigned id_{std::numeric_limits<unsigned>::max()};
};  // Material class
}  // mpm namespace

#include "material.tcc"

#endif  // MPM_MATERIAL_H_
