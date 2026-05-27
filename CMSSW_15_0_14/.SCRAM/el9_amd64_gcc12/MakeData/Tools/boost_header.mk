ALL_TOOLS      += boost_header
boost_header_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/boost/1.80.0-71a9f99e945ddf5b85afd9761a632291/include
boost_header_EX_USE := sockets root_cxxdefaults
boost_header_EX_FLAGS_CPPDEFINES  := -DBOOST_SPIRIT_THREADSAFE -DPHOENIX_THREADSAFE -DBOOST_MATH_DISABLE_STD_FPCLASSIFY -DBOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX
boost_header_EX_FLAGS_CXXFLAGS  := -Wno-error=unused-variable
boost_header_EX_FLAGS_SYSTEM_INCLUDE  := 1

