add_executable(sama main.cpp settings.cpp dbgraph.cpp dsnode.cpp
               kmernpp.cpp util.cpp debug.cpp coverage.cpp alignment.cpp
               graphaln.cpp refcomp.cpp correctgraph.cpp crfmult.cpp
               readaln.cpp threshold.cpp)

target_link_libraries(sama readfile kmer pgm pthread libdai gmp)

if (ZLIB_FOUND)
        target_link_libraries(sama ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)

install(TARGETS sama RUNTIME DESTINATION bin)
add_subdirectory(readfile)
add_subdirectory(kmer)
add_subdirectory(pgm)
add_subdirectory(libdai)

