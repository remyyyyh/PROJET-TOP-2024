add_library(stencil SHARED stencil/config.c stencil/comm_handler.c stencil/mesh.c stencil/init.c stencil/solve.c)
target_include_directories(stencil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(stencil PUBLIC m)

add_library(utils SHARED chrono.c)
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(stencil::stencil ALIAS stencil)
add_library(stencil::utils ALIAS utils)
