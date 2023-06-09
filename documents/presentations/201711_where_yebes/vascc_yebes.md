---
title: Benchmarking Where VLBI
subtitle: The VLBI Analysis Software Comparison Campaign (VASCC)
author: Team Where
date: November 27, 2017
---

# The VLBI Analysis Software Comparison Campaign

+ Done in 2015/2016 by Grzegorz Klopotek at Chalmers University (Sweden)

+ 11 different software packages compared

+ We took part using our legacy software _GEOSAT_

+ Have redone the analysis using **Where**


# The VLBI Analysis Software Comparison Campaign

+ Two networks of stations:

    + 4 stations on the southern hemisphere (SH)
    + 5 stations on the northern hemisphere (NH)

+ Virtual observations of one radio source for each network are scheduled every
  minute for 15 days

    + 129600 observations for SH
    + 216000 observations for NH

+ A leap second is introduced in the middle of the sessions (midnight June 30th 2015)


# The VLBI Analysis Software Comparison Campaign

\begin{align*}
     \tau =& \tau_{\text{geom}} + \tau_{\text{grav}} + \Delta\tau_{\text{tropo}}
                     + \Delta\tau_{\text{axis}} - \Delta\tau_{\text{therm}} \\
                 +& \tau_{\text{clock}} - \Delta\tau_{\text{cable}} + \tau_{\text{iono}} ,
\end{align*}

+ The following terms are included in the campaign:

    + Geometric and gravitational delay from the IERS 2010 Conventions
    + Troposphere (hydrostatic delay with GMF mapping function)
    + Thermal deformations (with constant temperatures)
    + Axis offset (as described by Nothnagel)

+ Delays due to cable calibration, the ionosphere and clocks are ignored.

+ Site displacement models are also compared.


# Results

+ 6 of the 11 softwares participating got sub-millimeter agreements.

+ Among these 6 softwares,

    + the biggest RMS of differences was 0.71 mm, and
    + the smallest RMS of differences was 0.21 mm.

+ The best agreement was found between _C5++_ and _VieVS_.


# C5++ vs VieVS, Northern

\begin{figure}[p]
  \includegraphics[width=\columnwidth]{figure/c5++_vs_vievs_NH}
  \caption{The RMS of all differences is 0.21~mm.}
\end{figure}

# C5++ vs VieVS, Southern

\begin{figure}[p]
  \includegraphics[width=\columnwidth]{figure/c5++_vs_vievs_SH}
  \caption{The RMS of all differences is 0.39~mm.}
\end{figure}


# Benchmarking Where

+ We have carried out the same analysis for **Where**, and compared to _C5++_
  and _VieVS_.

+ This has been a great exercise, both for finding bugs and trusting our
  analysis.

+ On the models tested in the campaign, **Where** gives very agreement with both
  _C5++_ and _VieVS_.

# Where vs C5++, Northern

\begin{figure}[p]
  \includegraphics[width=\columnwidth]{figure/where_vs_c5++_NH}
  \caption{The RMS of all differences is 0.49~mm.}
\end{figure}

# Where vs C5++, Southern

\begin{figure}[p]
  \includegraphics[width=\columnwidth]{figure/where_vs_c5++_SH}
  \caption{The RMS of all differences is 0.18~mm.}
\end{figure}

# Where vs VieVS, Northern

\begin{figure}[p]
  \includegraphics[width=\columnwidth]{figure/where_vs_vievs_NH}
  \caption{The RMS of all differences is 0.44~mm.}
\end{figure}

# Where vs VieVS, Southern

\begin{figure}[p]
  \includegraphics[width=\columnwidth]{figure/where_vs_vievs_SH}
  \caption{The RMS of all differences is 0.43~mm.}
\end{figure}


