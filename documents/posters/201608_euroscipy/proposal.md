---
title: EuroSciPy 2016 -- Proposal
author: Geir Arne Hjelle
date: May 6, 2016
---

# High precision positioning using Python

## Abstract

Using data from quasars, lasers and satellites the geodetic community is monitoring the everchanging surface of the
Earth. Traditionally, this has been done using Fortran code. At Kartverket we are currently using Python to develop
flexible, solid and fast software that can do this analysis, while also taking advantage of the experience of the old
softwares.

## Description

The presentation will consist of

+ an introduction to reference frames and geodetic techniques,
+ a description of the problem we are solving: combining huge datasets from different sources, and how Python's
  flexibility is helping,
+ experiences in using Python's scientific libraries and interfacing with legacy C and Fortran code,
+ a demo of the software and supporting tools.

## Detailed abstract

A combination of data from quasars (VLBI), lasers (SLR) and satellites (GPS and DORIS) is used to develop a coordinate
system for the Earth, which is necessary to monitor the changes of the Earth's surface: Land rising after the last ice
age, sea level increases and the tectonic plates moving several centimeters every year. Analysis from several research
institutes around the world is combined in order to establish this coordinate system, known as the International
Terrestrial Reference Frame (ITRF). Traditionally, most of this analysis has been done using Fortran code. At
Kartverket (the Norwegian Mapping Authority) we are currently developing new software in Python that will take part in
this international collaboration.

This presentation will show how Python is helping us effectively build flexible, solid and fast code for doing very
complex physical modeling. We are able to make use of well known data science packages like `numpy`, `scipy` and
`pandas`, in addition to lesser known, more specialized modules like `astropy` and `jplephem`. Furthermore, Python's
capabilities in interfacing with C and Fortran code makes it easy to take advantage of models the geodetic community
has made available.

More detailed, the presentation will consist of

+ a short introduction to the International Terrestrial Reference Frame (ITRF) and the different geodetic techniques
  (VLBI, SLR, GPS, and DORIS),

+ a short description of the problem we are solving: Combination of huge datasets from different sources, and how
  Python's flexibility is helping us write clean and well-organized code,

+ experiences in using `numpy`, `scipy` and `pandas`, as well as the more specialized `astropy` and `jplephem` packages
  for analysing geodetic data,

+ experiences in interfacing Python code with well established C and Fortran libraries already in use in the geodetic
  community,

+ a short demo of the software, including some of the tools we have built to ease the analysis and reporting of data.

## Tags

reference frame, geodetic analysis, quasars, lasers, satellites, python, fortran, astropy, jplephem, demo
