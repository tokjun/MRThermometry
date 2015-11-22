cmake_minimum_required(VERSION 2.8.9)

project(MRThermometry)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://www.example.com/Slicer/Extensions/MRThermometry")
set(EXTENSION_CATEGORY "IGT")
set(EXTENSION_CONTRIBUTORS "Junichi Tokuda (Brigham and Women's Hospital)")
set(EXTENSION_DESCRIPTION "MR Thermometry using proton resonant frequency shift.")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/MRThermometry.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/MRThermometry/Screenshots/1.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(PRFThermometry)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
