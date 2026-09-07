cmake_minimum_required(VERSION 3.24)

project(
  AstroCFA
  VERSION 0.1.0
  DESCRIPTION "CFA-aware, measurement-constrained RAW development for astrophotography"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Qt6 REQUIRED COMPONENTS Widgets)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(LibRaw QUIET IMPORTED_TARGET libraw)
  pkg_check_modules(LibTiff QUIET IMPORTED_TARGET libtiff-4)
  pkg_check_modules(LibJpeg QUIET IMPORTED_TARGET libjpeg)
endif()

set(ASTROCFA_CORE_DEPS)
if(TARGET PkgConfig::LibRaw)
  list(APPEND ASTROCFA_CORE_DEPS PkgConfig::LibRaw)
else()
  find_package(unofficial-libraw CONFIG REQUIRED)
  list(APPEND ASTROCFA_CORE_DEPS unofficial::libraw::libraw)
endif()

if(TARGET PkgConfig::LibTiff)
  list(APPEND ASTROCFA_CORE_DEPS PkgConfig::LibTiff)
else()
  find_package(TIFF REQUIRED)
  list(APPEND ASTROCFA_CORE_DEPS TIFF::TIFF)
endif()

if(TARGET PkgConfig::LibJpeg)
  list(APPEND ASTROCFA_CORE_DEPS PkgConfig::LibJpeg)
else()
  find_package(JPEG REQUIRED)
  list(APPEND ASTROCFA_CORE_DEPS JPEG::JPEG)
endif()

add_library(astrocfa_core
  src/background_model.cpp
  src/calibration.cpp
  src/cfa.cpp
  src/connected_components.cpp
  src/diagnostic_maps.cpp
  src/demosaic.cpp
  src/drizzle.cpp
  src/frequency_cfa.cpp
  src/image_writer.cpp
  src/joint_reconstruction.cpp
  src/luminance_proxy.cpp
  src/multiframe_benchmark.cpp
  src/noise_model.cpp
  src/output_transform.cpp
  src/psf_estimator.cpp
  src/raw_color.cpp
  src/raw_inspector.cpp
  src/raw_loader.cpp
  src/reconstruction_metrics.cpp
  src/registration.cpp
  src/star_detector.cpp
  src/synthetic_astro_scene.cpp
)

target_include_directories(astrocfa_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(astrocfa_core PUBLIC cxx_std_23)
target_link_libraries(astrocfa_core PUBLIC ${ASTROCFA_CORE_DEPS})

add_executable(astrocfa-nogui
  src/nogui_main.cpp
)

target_compile_features(astrocfa-nogui PRIVATE cxx_std_23)
target_link_libraries(astrocfa-nogui PRIVATE astrocfa_core)

add_executable(astrocfa
  src/gui_main.cpp
)

target_compile_features(astrocfa PRIVATE cxx_std_23)
target_link_libraries(astrocfa PRIVATE Qt6::Widgets astrocfa_core)

add_executable(astrocfa_cfa_tests
  tests/cfa_tests.cpp
)

target_compile_features(astrocfa_cfa_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_cfa_tests PRIVATE astrocfa_core)

add_executable(astrocfa_calibration_tests
  tests/calibration_tests.cpp
)

target_compile_features(astrocfa_calibration_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_calibration_tests PRIVATE astrocfa_core)

add_executable(astrocfa_component_tests
  tests/component_tests.cpp
)

target_compile_features(astrocfa_component_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_component_tests PRIVATE astrocfa_core)

add_executable(astrocfa_star_detector_tests
  tests/star_detector_tests.cpp
)

target_compile_features(astrocfa_star_detector_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_star_detector_tests PRIVATE astrocfa_core)

add_executable(astrocfa_demosaic_tests
  tests/demosaic_tests.cpp
)

target_compile_features(astrocfa_demosaic_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_demosaic_tests PRIVATE astrocfa_core)

add_executable(astrocfa_drizzle_tests
  tests/drizzle_tests.cpp
)

target_compile_features(astrocfa_drizzle_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_drizzle_tests PRIVATE astrocfa_core)

add_executable(astrocfa_registration_tests
  tests/registration_tests.cpp
)

target_compile_features(astrocfa_registration_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_registration_tests PRIVATE astrocfa_core)

add_executable(astrocfa_frequency_cfa_tests
  tests/frequency_cfa_tests.cpp
)

target_compile_features(astrocfa_frequency_cfa_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_frequency_cfa_tests PRIVATE astrocfa_core)

add_executable(astrocfa_image_writer_tests
  tests/image_writer_tests.cpp
)

target_compile_features(astrocfa_image_writer_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_image_writer_tests PRIVATE astrocfa_core)

add_executable(astrocfa_output_transform_tests
  tests/output_transform_tests.cpp
)

target_compile_features(astrocfa_output_transform_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_output_transform_tests PRIVATE astrocfa_core)

add_executable(astrocfa_diagnostic_maps_tests
  tests/diagnostic_maps_tests.cpp
)

target_compile_features(astrocfa_diagnostic_maps_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_diagnostic_maps_tests PRIVATE astrocfa_core)

add_executable(astrocfa_synthetic_benchmark_tests
  tests/synthetic_benchmark_tests.cpp
)

target_compile_features(astrocfa_synthetic_benchmark_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_synthetic_benchmark_tests PRIVATE astrocfa_core)

add_executable(astrocfa_joint_reconstruction_tests
  tests/joint_reconstruction_tests.cpp
)

target_compile_features(astrocfa_joint_reconstruction_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_joint_reconstruction_tests PRIVATE astrocfa_core)

add_executable(astrocfa_raw_color_tests
  tests/raw_color_tests.cpp
)

target_compile_features(astrocfa_raw_color_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_raw_color_tests PRIVATE astrocfa_core)

add_executable(astrocfa_background_model_tests
  tests/background_model_tests.cpp
)

target_compile_features(astrocfa_background_model_tests PRIVATE cxx_std_23)
target_link_libraries(astrocfa_background_model_tests PRIVATE astrocfa_core)

if(MSVC)
  target_compile_options(astrocfa_core PRIVATE /W4)
  target_compile_options(astrocfa-nogui PRIVATE /W4)
  target_compile_options(astrocfa PRIVATE /W4)
  target_compile_options(astrocfa_cfa_tests PRIVATE /W4)
  target_compile_options(astrocfa_calibration_tests PRIVATE /W4)
  target_compile_options(astrocfa_component_tests PRIVATE /W4)
  target_compile_options(astrocfa_star_detector_tests PRIVATE /W4)
  target_compile_options(astrocfa_demosaic_tests PRIVATE /W4)
  target_compile_options(astrocfa_drizzle_tests PRIVATE /W4)
  target_compile_options(astrocfa_registration_tests PRIVATE /W4)
  target_compile_options(astrocfa_frequency_cfa_tests PRIVATE /W4)
  target_compile_options(astrocfa_image_writer_tests PRIVATE /W4)
  target_compile_options(astrocfa_output_transform_tests PRIVATE /W4)
  target_compile_options(astrocfa_diagnostic_maps_tests PRIVATE /W4)
  target_compile_options(astrocfa_synthetic_benchmark_tests PRIVATE /W4)
  target_compile_options(astrocfa_joint_reconstruction_tests PRIVATE /W4)
  target_compile_options(astrocfa_raw_color_tests PRIVATE /W4)
  target_compile_options(astrocfa_background_model_tests PRIVATE /W4)
else()
  target_compile_options(astrocfa_core PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa-nogui PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_cfa_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_calibration_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_component_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_star_detector_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_demosaic_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_drizzle_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_registration_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_frequency_cfa_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_image_writer_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_output_transform_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_diagnostic_maps_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_synthetic_benchmark_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_joint_reconstruction_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_raw_color_tests PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(astrocfa_background_model_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()

add_test(NAME astrocfa_nogui_help COMMAND astrocfa-nogui --help)
add_test(NAME astrocfa_nogui_version COMMAND astrocfa-nogui --version)
add_test(NAME astrocfa_nogui_inspect_requires_input COMMAND astrocfa-nogui inspect)
set_tests_properties(astrocfa_nogui_inspect_requires_input PROPERTIES WILL_FAIL TRUE)
add_test(NAME astrocfa_cfa_tests COMMAND astrocfa_cfa_tests)
add_test(NAME astrocfa_calibration_tests COMMAND astrocfa_calibration_tests)
add_test(NAME astrocfa_component_tests COMMAND astrocfa_component_tests)
add_test(NAME astrocfa_star_detector_tests COMMAND astrocfa_star_detector_tests)
add_test(NAME astrocfa_demosaic_tests COMMAND astrocfa_demosaic_tests)
add_test(NAME astrocfa_drizzle_tests COMMAND astrocfa_drizzle_tests)
add_test(NAME astrocfa_registration_tests COMMAND astrocfa_registration_tests)
add_test(NAME astrocfa_frequency_cfa_tests COMMAND astrocfa_frequency_cfa_tests)
add_test(NAME astrocfa_image_writer_tests COMMAND astrocfa_image_writer_tests)
add_test(NAME astrocfa_output_transform_tests COMMAND astrocfa_output_transform_tests)
add_test(NAME astrocfa_diagnostic_maps_tests COMMAND astrocfa_diagnostic_maps_tests)
add_test(NAME astrocfa_synthetic_benchmark_tests COMMAND astrocfa_synthetic_benchmark_tests)
add_test(NAME astrocfa_joint_reconstruction_tests COMMAND astrocfa_joint_reconstruction_tests)
add_test(NAME astrocfa_raw_color_tests COMMAND astrocfa_raw_color_tests)
add_test(NAME astrocfa_background_model_tests COMMAND astrocfa_background_model_tests)
