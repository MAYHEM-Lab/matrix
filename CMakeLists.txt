cmake_minimum_required(VERSION 3.13)
project(matrix)

find_package(lapacke REQUIRED)
find_package(distributions REQUIRED)
find_package(woof REQUIRED)

add_library(mioarray mioarray.c)
target_link_libraries(mioarray PUBLIC mio)

add_library(polyco polyco.c)
target_link_libraries(polyco PUBLIC mio)

add_library(mioregress mioregress.c)
target_link_libraries(mioregress PUBLIC mioarray distributions)

install(TARGETS mioarray mioregress EXPORT matrix DESTINATION lib)
install(EXPORT matrix DESTINATION lib/cmake/matrix)

file(GLOB matrix_headers "*.h")

install(FILES ${matrix_headers} DESTINATION include)
install(FILES matrix-config.cmake DESTINATION lib/cmake/matrix)
