gridkit_add_library(
  emt_pwm
  INTERFACE_TARGET
  HEADERS Pwm.hpp PwmData.hpp
  LINK_LIBRARIES INTERFACE GridKit::emt_core)
