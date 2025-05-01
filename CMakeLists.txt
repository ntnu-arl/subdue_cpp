cmake_minimum_required(VERSION 3.0.2)
project(subdue)

find_package(catkin_simple REQUIRED)
find_package(GSL REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_definitions(-Wall -O3 -Werror=return-type) # -O3

#############
# LIBRARIES #
#############
cs_add_library(${PROJECT_NAME}
  src/compress.c
  src/discover.c
  src/dot.c
  src/evaluate.c
  src/extend.c
  src/graphmatch.c
  src/graphops.c
  src/labels.c
  src/sgiso.c
  src/subops.c
  src/test.c
  src/utility.c
  src/avl.c
  src/gendata.c
  src/incboundary.c
  src/inccomp.c
  src/incextend.c
  src/incgraphops.c
  src/incutil.c
  src/subdue.c
)

include_directories(include ${catkin_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} GSL::gsl GSL::gslcblas)

# catkin_package(
#   INCLUDE_DIRS include ${catkin_INCLUDE_DIRS}
# )

# subs2dot: subs2dot_main.o $(OBJS)
# 	$(CC) $(LDFLAGS) -o subs2dot subs2dot_main.o $(OBJS) $(LDLIBS) $(GSLLIBS)

cs_add_executable(subdue_node src/main.c)
target_link_libraries(subdue_node ${PROJECT_NAME} GSL::gsl GSL::gslcblas)

cs_add_executable(subs2dot src/subs2dot_main.c)
target_link_libraries(subs2dot ${PROJECT_NAME} GSL::gsl GSL::gslcblas)

cs_add_executable(subdue_cpp_node src/main.cpp)
target_link_libraries(subdue_cpp_node ${PROJECT_NAME} GSL::gsl GSL::gslcblas)

##########
# EXPORT #
##########
cs_install()
cs_export(LIBRARIES)
