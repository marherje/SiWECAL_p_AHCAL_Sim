# Install script for directory: /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so.0.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so.0.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17/lib64:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/ilcutil/v01-06-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES
    "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build/lib/libLCIO2BuildProcessor.so.0.1.0"
    "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build/lib/libLCIO2BuildProcessor.so.0.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so.0.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so.0.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17/lib64:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/ilcutil/v01-06-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
           NEW_RPATH "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17/lib64:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/ilcutil/v01-06-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/cvmfs/sft.cern.ch/lcg/releases/binutils/2.30-e5b21/x86_64-centos7/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so"
         RPATH "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17/lib64:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/ilcutil/v01-06-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build/lib/libLCIO2BuildProcessor.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so"
         OLD_RPATH "/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17/lib64:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/ilcutil/v01-06-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17/lib64:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/ilcutil/v01-06-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/cvmfs/sft.cern.ch/lcg/releases/binutils/2.30-e5b21/x86_64-centos7/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libLCIO2BuildProcessor.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
