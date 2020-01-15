"""Write UERE analysis results

Description:
------------



"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config


@plugins.register
def uere_writer(dset):
    """Write UERE analysis results

    Args:
        dset:       Dataset, a dataset containing the data.
    """

    with config.files.open("output_uere", file_vars=dset.vars, mode="wt") as fid:

        # Write header
        fid.write(
            "#{:>19s}{:>14s}{:>5s}{:>16s}{:>16s}{:>16s}\n"
            "".format("YYYY:MM:DD:hh:mm:ss", "MJD", "SAT", "SISRE", "UEE", "UERE")
        )
        fid.write("#{:>54s}{:>16s}{:>16s}\n" "".format("[m]", "[m]", "[m]"))
        fid.write("#{}\n".format("_" * 88))

        # Loop over all observations
        for idx in range(0, dset.num_obs):

            #
            # Write SISRE analysis results
            #

            # #YYYY:MM:DD:hh:mm:ss           MJD  SAT           SISRE             UEE            UERE
            # #                                                   [m]             [m]             [m]
            # #______________________________________________________________________________________
            #  2017:10:28:00:00:00  58054.000000  G01         0.86957        -1.26053        -0.51241
            #  2017:10:28:00:00:00  58054.000000  G02        -0.27013         0.25391        -0.39711
            # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+--
            fid.write(
                " {:>19s}{:14.6f}{:>5s}{:16.5f}{:16.5f}{:16.5f}\n"
                "".format(
                    dset.time.gps.datetime[idx].strftime("%Y:%m:%d:%H:%M:%S"),
                    dset.time.gps.mjd[idx],
                    dset.satellite[idx],
                    dset.sisre[idx],
                    dset.uee[idx],
                    dset.uere[idx],
                )
            )
