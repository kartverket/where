---
title: Where
subtitle: A new geodetic software being developed at the Norwegian Mapping Authority
author: Team Where
date: November 27, 2017
---

# A short history

The **Where** activity was started in the fall of 2015 with the goal of
building software that can analyse (and combine) data for _VLBI_, _SLR_, _GNSS_
(and _DORIS_).

+ **Where** builds on ideas and experiences from the _Geosat_ software

+ The **Where-team** consists of four researchers at the Norwegian Mapping
  Authority (NMA):

    + Michael Dähnn (GNSS)
    + Ingrid Fausk (SLR)
    + Geir Arne Hjelle (VLBI/SLR/GNSS)
    + Ann-Silje Kirkvik (VLBI)

+ NMA is currently building / expanding the observatory at Ny-Ålesund


# Current status

+ All models from the _IERS Conventions_ are implemented for _VLBI_ and _SLR_

    + Many of the models can be reused for the other techniques

    + In 2015 we participated in a _VLBI Analysis Software Comparison Campaign_
      (_VASCC_) organized by Onsala with _Geosat_. Later we have benchmarked
      **Where** against these results (more later)

    + We have implemented an orbit integrator for _SLR_ satellites, and are
      working on improving it

+ For _GNSS_ we can handle data from _GPS_ and _Galileo_, and perform several
  types of analysis. At the moment we are focusing on supporting a project for
  doing independent monitoring of the _Galileo_-system.


# Technology

The **Where** software is mainly being written in _Python_

+ Cross-platform: Runs on Linux, Mac, Windows

+ Solid, flexible and fast libraries like `numpy`, `astropy`, `matplotlib` and
  `scipy` are available

+ We use a **HDF5**-based format for storing data while the program is running

+ A companion tool **There** is used for visualizing and analyzing the results
  from **Where**

+ _Python_ has effective interfaces to _C_ and _Fortran_ code, and we use the
  **Sofa** and **IERS** software libraries directly


# Implementation

**Where** uses the idea of a _pipeline_: Data are handled in several steps,
where the result of one step is input to the next step.

+ **Read**: Read observation data from file

+ **Edit**: Clean the data, remove observations that should be ignored

+ **Orbit**: Calculate / read satellite orbits

+ **Calculate**: Calculate theoretical delay for all models

    + Includes both site position and delay models
    + Output: _prefit residual_ = _obs - calc_

+ **Estimate**: Estimate troposphere, clock, positions, EOP etc

    + Output: _postfit residual_ = _obs - (calc + est)_

+ **Write**: Output to Sinex and others


# Future plans

At the moment, the highest priorities for **Where** are

+ finishing the _VLBI_ analysis

    + Setting up operational VLBI-analysis at NMA

    + Becoming an IVS analysis center

+ finishing the _SLR_ analysis

    + Cleaning up the orbit integrator
    
    + Speeding up the analysis

+ support different kinds of _GNSS_-analysis dependent on other projects

    + Monitoring the quality of broadcast orbits

    + Assessing the quality of different ionosphere models


# Future plans (continued)

At the moment, the highest priorities for **Where** are

+ exposing parts of **Where** as a library that can be reused for other
  projects

    + Handling of different data formats, time scales, reference frames,
      satellite orbits, etc

+ preparing for the first release of **Where**

    + We plan to release the software as an open source project available to
      anyone interested


# Quick demonstration
