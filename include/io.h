#ifndef MPM_IO_H_
#define MPM_IO_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>

#include <boost/filesystem.hpp>

#include "tclap/CmdLine.h"
//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;
//! Speed Logger
#include "spdlog/spdlog.h"

//! \brief Input/Output handler
class IO {
 public:
  //! Constructor with argc and argv
  //! \param[in] argc Number of input arguments
  //! \param[in] argv Input arguments
  IO(int argc, char** argv);

  //! Return input mesh file name
  std::string file_name(const std::string& file);

 private:
  //! Working directory
  std::string working_dir_;
  //! Input file name
  std::string input_file_{"mpm.json"};
  //! Input JSON object
  Json json_;
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};

#endif  // MPM_IO_H_
