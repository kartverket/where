#!/usr/bin/env python3
"""Compare Where and gLAB analysis

Description:
------------
1. Run gLAB GNSS analysis
cd ~/SOURCE/GNSS/gLAB_3.0.0_Linux

./gLAB_linux  -input:cfg gLAB_broadcast_spp.cfg
./gLAB_linux  -input:cfg gLAB_broadcast.cfg
./gLAB_linux  -input:cfg gLAB_precise.cfg

grep -E 'MODEL.*C1C' broadcast_spp.out > ~/where/src/analysis/compare_where_glab/20180201_glab_model_brdc_c1c_spp.out
grep -E 'MODEL.*C1C' broadcast.out > ~/where/src/analysis/compare_where_glab/20180201_glab_model_brdc_c1c.out
grep -E 'MODEL.*C1C' precise.out > ~/where/src/analysis/compare_where_glab/20180201_glab_model_precise_c1c.out
 
grep -E 'FILTER   20' broadcast_spp.out > ~/where/src/analysis/compare_where_glab/20180201_glab_filter_brdc_c1c_spp.out
grep -E 'FILTER   20' broadcast.out > ~/where/src/analysis/compare_where_glab/20180201_glab_filter_brdc_c1c.out
grep -E 'FILTER   20' precise.out > ~/where/src/analysis/compare_where_glab/20180201_glab_filter_precise_c1c.out

grep -E 'OUTPUT 20' broadcast_spp.out > ~/where/src/analysis/compare_where_glab/20180201_glab_output_brdc_c1c_spp.out
grep -E 'OUTPUT 20' broadcast.out > ~/where/src/analysis/compare_where_glab/20180201_glab_output_brdc_c1c.out
grep -E 'OUTPUT 20' precise.out > ~/where/src/analysis/compare_where_glab/20180201_glab_output_precise_c1c.out

#NOTE: gLAB configuration file can be generated via gLAB_GUI.py.


2. Run Where GNSS analysis

where 2016 3 1 -g

GENERAL CONFIGURATION:
[gnss]
sampling_rate                  = 300
systems                        = G
use_mixed_brdc_file            = True
brdc_block_nearest_to          = toe_positive
calc_models                    = gnss_satellite_clock, gnss_relativistic_clock, gnss_range, troposphere_radio
                                 #or for precise orbits in addition: gnss_satellite_phase_center_offset
estimate_stochastic            = gnss_rcv_clock
estimate_constant              = gnss_site_pos
apriori_orbit                  = broadcast    #or precise

[gnss_rcv_clock]
knot_interval                  = 1
process_noise                  = 9e10
apriori_stdev                  = 1
apriori_rate_stdev             = 1
unit                           = second
display_unit                   =

[gnss_site_pos]
fix_stations                   =
knot_interval                  =
process_noise                  =
apriori_stdev                  = 1
apriori_rate_stdev             = 1
unit                           = meter
display_unit                   =

[troposphere_radio]
mapping_function               = gpt2w


CONFIGURATION files.conf:
[gnss_rinex_obs]
filename        = {$STATION}00NOR_R_{$yyyy}{$doy}0000_01D_300s_MO.rnx{$gz}
directory       = {$path_data}/obs/gnss/{$yyyy}/{$doy}

[gnss_rinex_nav_M]
filename        = brdm{$doy}0.{$yy}p{$gz}
directory       = {$path_data}/obs/orb/brdc/{$yyyy}

[gnss_orbit_sp3]
filename        = com{$gpsweek}{$dow}.sp3
directory       = {$path_data}/obs/orb/igs/{$gpsweek}

[gnss_rnx_clk]
filename        = com{$gpsweek}{$dow}.clk
directory       = {$path_data}/obs/orb/igs/{$gpsweek}

[gnss_antex]
filename        = igs14.atx
directory       = {$path_data}/apriori/antenna


3. Read gLAB output file and Where Dataset

4. Generate gLAB Dataset

5. Compare Where and gLAB Dataset

"""
# Standard library imports
import subprocess

# Where imports
from where import data
from where import parsers
from where.lib import log

TECH = "gnss"


def _difference_where_glab(dwhere, dglab, stage):

    if stage == "model":
        difference_by = ["time", "satellite"]
    else:
        difference_by = ["time"]

    if dwhere.num_obs == 0 or dglab.num_obs == 0:
        log.warn(
            f"Nothing to compare. Number of observations are zero at least for one dataset "
            f"(dwhere: {dwhere.num_obs}, dglab: {dglab.num_obs})."
        )

    ddiff = dwhere.difference_with(dglab, difference_by=difference_by)
    ddiff.write_as(tech=TECH, stage="dwhrglb")


def _read_glab_file(inpath, rundate, pipeline, station):
    """Generate gLAB dataset

    Read gLAB output file.
    """

    # Read file
    p = parsers.parse_file(parser_name="glab_output", file_path=inpath)

    # Get dataset
    dset = p.write_to_dataset()

    # Write dataset
    dset.write_as(pipeline=pipeline, stage="glab")

    return dset


def _generate_where_dset(rundate, pipeline, station, stage):
    """Generate Where dataset

    Start Where analysis with Where and read Dataset.
    """
    program = "where"
    system = "G"

    # Start processing with Where
    where_args = [
        str(rundate.year),
        str(rundate.month),
        str(rundate.day),
        "--" + pipeline,
        "--station=" + station,
        "--systems=" + system,
        "--calc_models=gnss_satellite_clock, gnss_range, gnss_relativistic_clock, gnss_ionosphere, troposphere_radio, gnss_total_group_delay",
        "--sampling_rate=300",
        "--estimate_epochwise=True",
        "--elevation:cut_off=5",
        "--profile=cnes_lnav_l1",
        "-N",
        "-D",
        "-T",
    ]
    log.info("Start {}: {} {}".format(pipeline.upper(), program, " ".join(where_args)))
    process = subprocess.run([program] + where_args)
    if process.returncode:
        log.error("Where {} failed with error code {} ({})", pipeline.upper(), process.returncode, " ".join(process.args))

    # Read and write GNSS Dataset
    dset = data.Dataset(rundate=rundate, pipeline=pipeline, stage=stage, station=station, label="None", id="")
    dset.write_as(pipeline=pipeline, stage="where")

    return dset


############################################################################################
#
#
#                                MAIN PROGRAM
#
#
############################################################################################
if __name__ == "__main__":

    log.init(log_level="info")
    year = 2018
    month = 2
    day = 1
    rundate = date(year, month, day)
    station = "stas"
    pipeline ="gnss"

    # Read gLAB output file
    stage = "output"  # "filter", "output"
    where_stage = {"model": "calculate", "filter": "estimate", "output": "write"}
    inpath = "{:%Y%m%d}_glab_{}_brdc_c1c_spp.out".format(rundate, stage)
    # inpath = "{:%Y%m%d}_glab_{}_brdc_c1c.out".format(rundate, stage)
    # inpath = '{:%Y%m%d}_glab_{}_precise_c1c.out'.format(rundate, stage)

    # Where dataset
    dwhere = _generate_where_dset(rundate, pipeline, station, where_stage[stage])

    # gLAB dataset
    dglab = _read_glab_file(inpath, rundate, pipeline, station)

    # Determine difference between the Where and gLAB dataset
    _difference_where_glab(dwhere, dglab, stage)
