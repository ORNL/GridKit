# [[
#  EMT frequency-domain parameter models and the Overhead aggregate.
# ]]

gridkit_add_library(
  emt_parameter_views
  INTERFACE_TARGET
  HEADERS Views.hpp
  LINK_LIBRARIES INTERFACE GridKit::emt_core)

add_subdirectory(Geometry/Conductor)
add_subdirectory(Geometry/Tower)
add_subdirectory(Geometry/Path)
add_subdirectory(Effects/GeometricInductance)
add_subdirectory(Effects/SkinEffect)
add_subdirectory(Effects/Carson)
add_subdirectory(Effects/SeriesImpedance)
add_subdirectory(Effects/ShuntPotential)
add_subdirectory(Effects/ShuntAdmittance)
add_subdirectory(Effects/Reduction)
add_subdirectory(Modal)
add_subdirectory(Response/Gamma)
add_subdirectory(Response/Yc)
add_subdirectory(Response/Zc)

gridkit_add_library(
  emt_overhead
  INTERFACE_TARGET
  HEADERS Overhead.hpp
          OverheadData.hpp
          OverheadDataJSONParser.hpp
          OverheadImpl.hpp
  LINK_LIBRARIES
    INTERFACE
    GridKit::emt_core
    INTERFACE
    GridKit::emt_tower
    INTERFACE
    GridKit::emt_path
    INTERFACE
    GridKit::emt_geometric_inductance
    INTERFACE
    GridKit::emt_skin_effect
    INTERFACE
    GridKit::emt_carson
    INTERFACE
    GridKit::emt_series_impedance
    INTERFACE
    GridKit::emt_shunt_potential
    INTERFACE
    GridKit::emt_shunt_admittance
    INTERFACE
    GridKit::emt_reduction
    INTERFACE
    GridKit::emt_gamma
    INTERFACE
    GridKit::emt_modal
    INTERFACE
    GridKit::emt_yc
    INTERFACE
    GridKit::emt_zc
    INTERFACE
    GridKit::emt_conductor
    INTERFACE
    GridKit::dense_vector)
