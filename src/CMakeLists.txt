# enable C++11 support
find_package (CXX11)

# search SeqAn library    
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# enable SeqAn dependencies
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")
add_definitions (${SEQAN_DEFINITIONS})
include_directories (${SEQAN_INCLUDE_DIRS})

# register our target
add_executable (bam_realigner
                bam_realigner.cpp
                bam_realigner_app.cpp
                bam_realigner_app.h
                bam_realigner_options.h
                bam_realigner_options.cpp
                realigner_step.h
                realigner_step.cpp)
target_link_libraries (bam_realigner ${SEQAN_LIBRARIES})
