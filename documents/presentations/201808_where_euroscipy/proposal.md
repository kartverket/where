---
title: EuroSciPy 2018 - Proposal
author:
    - Geir Arne Hjelle
    - Michael Dähnn
    - Ingrid Fausk
    - Ann-Silje Kirkvik
date: May 13, 2018
---

# Listening to Quasars and Shooting Satellites With Lasers

## Abstract

Quasars, lasers and satellites are all used to keep track of the Earth as it tumbles through space. By combining observations of far away objects, some billions of light-years away, we are able to monitor the centimeter changes at the surface of our planet. Python plays an increasing part in this analysis.

## Description

Remote sensing from satellites is used to monitor changes of the Earth's surface: land rising after the last ice age, sea level increases and continents drifting apart. To be able to this, a common and accurate coordinate system for the Earth is needed. This system is called the International Terrestrial Reference Frame (ITRF) and is constructed and maintained through the joint effort of many research institutes around the world.

At the Norwegian Mapping Authority we are taking part in this analysis with our own Python software called Where and its support library, Midgard. Mainly aimed at the geodetic community, we hope that some of the tools will be useful also for a wider audience. Both Where and Midgard are released as open source packages.

This presentation will give a very brief introduction to the International Terrestrial Reference Frame (ITRF) and the geodetic techniques used to construct it. Next, we will dive into how we have been able to leverage the Python data analysis stack, including `numpy`, `pandas`, `scipy`, `astropy`, `h5py`, and `cython` in our analysis. We will end with some of the experiences we have made in building high-performant analysis pipelines in Python.

## Abstract as Tweet

How Python is used to combine data from quasars, lasers and satellite to monitor the Earth's ever-changing surface.

## Biography

Geir Arne works as a researcher/data scientist at the Norwegian Mapping Authority, trying to answer the big question: Where are we?

Playing with numbers since he was little, Geir Arne followed his interest and ended up studying and teaching mathematics in Trondheim, Barcelona and St. Louis. After an interlude pondering energy markets, he currently masquerades as a geodesist. Outside the office he teaches kids how to code and writes tutorials at realpython.com.

Other interests include board games, photography, hiking, square roots and uphill skiing. 
