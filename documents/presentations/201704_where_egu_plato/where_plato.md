---
title: Where
subtitle: A new geodetic software being developed at the Norwegian Mapping Authority
author: Geir Arne Hjelle and Eirik Mysen
date: April 27, 2017
---

# A short history

The **Where** project was started in the fall of 2015 with the goal of building
software that can analyse and combine data for _VLBI_, _SLR_, _GNSS_ and _DORIS_.

+ **Where** builds on ideas and experiences from the _GEOSAT_ software

+ The **Where-team** consists of five researchers at the Norwegian Mapping
  Authority (NMA):

    + Michael Dähnn (GPS)
    + Ingrid Fausk (SLR)
    + Geir Arne Hjelle (VLBI, GPS)
    + Ann-Silje Kirkvik (VLBI)
    + Eirik Mysen (VLBI)

+ NMA is currently building / expanding the observatory at Ny-Ålesund

# The new Ny-Ålesund observatory

![The twin telescopes will open next year, SLR to come later](nyal_antennas.jpg)

# Current status

+ All models from the _IERS Conventions_ are implemented for _VLBI_ and _SLR_

    + The _VLBI_ delay models are consistent with other softwares to the
      millimeter-level (tested against the _VASCC_ data presented last year at IVS)

    + Many of the models can be reused for the other techniques

+ Some work remains on the orbit integrator for _SLR_ satellites

+ We are currently working on Precise Point Positioning (PPP) for _GPS_

+ _DORIS_ is put on hold until we have resources to focus on it


# Technology

The **Where** software is mainly being written in _Python_

+ Solid, flexible and fast libraries like `numpy`, `astropy`, `matplotlib` and
  `scipy` are available

+ We use a **HDF5**-based format for storing data while the program is running

+ _Python_ has effective interfaces to _C_ and _Fortran_ code, and we can use
  the **Sofa** and **IERS** software libraries directly


# Technology -- plans

![The planned architecture](../../../monitor/figure/code_structure.pdf)


# Future plans

At the moment, the highest priorities for **Where** are

+ finishing the _VLBI_ analysis

    + Resolving the last issues in estimation

    + Doing operational testing

+ finishing the _SLR_ analysis

    + The orbit integrator needs some more work

+ finishing PPP for _GPS_ and starting to look at _Galileo_ and possibly
  _Glonass_

    + Orbit integration for GNSS-satellites
