Fitting your model
----------------------------------

Code
~~~~

.. code-block:: console

  $ make check

.. code-block:: cmake

   cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

   project(robestApp CXX)

   set(CMAKE_CXX_STANDARD 11)
   
   include_directories(${PROJECT_SOURCE_DIR}/example)

   add_executable(robest-line-fitting 
   example/robest-line-fitting.cpp
   example/LineFitting/LineFitting.hpp
   example/LineFitting/LineFitting.cpp)

Explanation
~~~~~~~~~~~
