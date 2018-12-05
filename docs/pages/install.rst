==================================
Install
==================================

There is no the installation procedure, you just need to clone the repository or download a zip file, and then add the necessary file to your project.

- Cloning repository for Linux systems:

.. code-block:: terminale

    git clone --recursive https://github.com/avkudr/robest.git

- Downloading .zip file for Windows:

1. On GitHub, navigate to the main page of the **robest** repository.

2. Under the repository name, click Clone or download. And choose download zip. 

.. image:: images/clone-repo-clone-url-button.png


Add to your project
___________________

cmake   
~~~~~

A small example of filling a *CMAKE* file and including a **robest** library for a simple project:

.. code-block:: cmake

    cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

    project(robestApp CXX)

    set(CMAKE_CXX_STANDARD 11)

    include_directories(${PROJECT_SOURCE_DIR}/example)

    add_executable(robest-line-fitting
    example/robest-line-fitting.cpp
    example/LineFitting/LineFitting.hpp
    example/LineFitting/LineFitting.cpp)

In this case, the *robestApp* project uses the line fitting class, which requires the presence of the source file of **robest** library.
