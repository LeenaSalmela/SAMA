add_library(libdai alldai.cpp bp.cpp bipgraph.cpp daialg.cpp exceptions.cpp factor.cpp factorgraph.cpp graph.cpp jtree.cpp mf.cpp parallelbp.cpp util.cpp) 
find_package(Boost 1.71 REQUIRED)
target_include_directories(libdai PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS} ${GMP_DIR})
target_link_libraries(libdai ${Boost_LIBRARIES} -L${GMP_DIR}/lib)
