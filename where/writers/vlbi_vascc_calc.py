"""Write a list of calculated observations in the VASCC format

Description:
------------

Write a list of calculated observations in a format specified by the VLBI Analysis Software Comparison Campaign (VASCC)
2015 document [1]_. The calc-variable is taken from the given :class:`~where.data.dataset.Dataset` (assumed to be
meter) and converted to seconds before written to the textfile. The following columns of information are stored:

* Observation number, counter from 1 to number of observations
* Date of observation (UTC), from ``time`` field
* Time of observation (UTC), from ``time`` field
* Name of source, from ``source`` field
* Name of station 1, from ``station_1`` field
* Name of station 2, from ``station_2`` field
* Calculated delay in seconds, converted from ``calc`` field

Example:
--------

.. highlight:: none

The following example shows an example of the first 10 lines of the output::

     1 2015/03/14 00:00:10.00 1424-418 HOBART12 YARRA12M -4.40034976520716E-03
     2 2015/03/14 00:01:18.00 1749+096 HOBART12 YARRA12M -8.35102410861801E-03
     3 2015/03/14 00:02:20.00 2145+067 HOBART12 YARRA12M +7.08658005475699E-04
     4 2015/03/14 00:03:37.00 1958-179 HOBART12 YARRA12M -1.96341421641893E-03
     5 2015/03/14 00:05:25.00 0208-512 HOBART12 YARRA12M +9.11774803487017E-03
     6 2015/03/14 00:06:48.00 1057-797 HOBART12 YARRA12M +3.37566164295730E-03
     7 2015/03/14 00:08:31.00 2255-282 HOBART12 YARRA12M +5.34818723943018E-03
     8 2015/03/14 00:10:19.00 NRAO530  HOBART12 YARRA12M -7.29513514203987E-03
     9 2015/03/14 00:10:19.00 NRAO530  KATH12M  YARRA12M -2.67561482753287E-03
    10 2015/03/14 00:12:05.00 CTA102   HOBART12 YARRA12M +1.79053750343479E-03


References:
-----------

.. [1] VLBI analysis software comparison campaign 2015, Grzegorz Klopotek.
       http://www.hobiger.org/VASCC2015/VASCC_INFO.pdf

"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from midgard.math.constant import constant
from where.lib import files


@plugins.register
def vascc_calc(dset):
    """Write a list of calculated observations in the VASCC format.

    Args:
        dset:  Dataset, data for a model run.
    """
    with files.open("output_vascc_calc", file_vars=dict(session=dset.dataset_name, **dset.vars), mode="wt") as fid:
        for obs, (time, src, sta_1, sta_2, calc) in enumerate(
            dset.values("time", "source", "station_1", "station_2", "calc"), start=1
        ):
            time_str = time.utc.datetime.strftime("%Y/%m/%d %H:%M:%S.%f")[:22]
            fid.write(
                "{:6d} {:>22s} {:<8s} {:<8s} {:<8s} {:+16.14E}\n".format(
                    obs, time_str, src, sta_1, sta_2, calc / constant.c
                )
            )
