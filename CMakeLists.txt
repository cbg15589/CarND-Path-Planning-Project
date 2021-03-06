project(Path_Planning)
cmake_minimum_required (VERSION 3.5)


# compiler flags.
add_definitions(-std=c++11)
set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

file(GLOB cpp "src/*.cpp")
file(GLOB headers "include/*.h*")
set(sources ${cpp} ${headers})

include_directories(include)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

    set(vcpkg_dir "C:/vcpkg")
    set(CMAKE_TOOLCHAIN_FILE "${vcpkg_dir}./scripts/buildsystems/vcpkg.cmake")

    include_directories("${vcpkg_dir}./installed/x86-windows/include")
    link_directories("${vcpkg_dir}./installed/x86-windows/lib")

endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    include_directories(/usr/local/include)
    include_directories(/usr/local/opt/openssl/include)
    link_directories(/usr/local/lib)
    link_directories(/usr/local/opt/openssl/lib)
    link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

add_executable(path_planning ${sources})

# linking (linux)
if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(path_planning  z ssl uv uWS)
endif (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# linking (windows)
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(path_planning zlib ssleay32 libuv uWS)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")