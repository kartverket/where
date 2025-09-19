"""Write RINEX navigation file

Description:
------------
Write data in the Rinex navigation file format 3.03 (see :cite:`rinex3`).

If the writer is called by RINEX_NAV analysis, then the given Dataset includes already necessary apriori broadcast
ephemeris information. If the writer is called from another analysis (e.g. SISRE), then the Dataset does not include
apriori broadcast ephemeris. In this case the given Dataset has to be overwritten with a apriori broadcast ephemeris
Dataset.


"""
# Standard library imports
from datetime import datetime, timedelta
from datetime import time as dt_time

# Midgard imports
from midgard.dev import plugins

# Where imports
import where
from where import apriori
from where.lib import config, log, util

# TODO: SYSTEM_TIME_OFFSET_TO_GPS_SECOND & SYSTEM_TIME_OFFSET_TO_GPS_WEEK should be placed in constans.conf

# The constant shows the 'second' relation between the GNSS time systems to the GPS time scale. Galileo (E),
# QZSS (I) and IRNSS (J) uses the same second as the GPS time scale, but BeiDou (C) BDT time system is 14 seconds behind
# GPS time.
SYSTEM_TIME_OFFSET_TO_GPS_SECOND = dict(C=14, E=0, I=0, J=0)

# The constant shows the relation between the GNSS week to the GPS week. Galileo (E), IRNSS (I) and QZSS (J) week
# corresponds to GPS week, whereas BeiDou (C) week starts at GPS week 1356.
SYSTEM_TIME_OFFSET_TO_GPS_WEEK = dict(C=1356, E=0, I=0, J=0)


@plugins.register
def rinex3_nav(dset: "Dataset"):
    """Write RINEX navigation file

    Args:
        dset:  Dataset, a dataset containing the data.
    """

    # Overwrite Dataset. This is necessary if the writer is called from a analysis (e.g. SISRE) with does not include
    # broadcast ephemeris information.
    # TODO: Is that the best solution?
    if dset.vars["pipeline"] != "rinex_nav":
        brdc = apriori.get(
            "orbit",
            rundate=dset.analysis["rundate"],
            system=tuple(dset.unique("system")),
            station=dset.vars["station"],
            apriori_orbit="broadcast",
        )
        meta = brdc.dset_edit.meta[dset.analysis["rundate"].strftime("%Y-%m-%d")]
        data = brdc.dset_edit  # TODO: Another possibility: brdc.dset_raw

    else:
        meta = dset.meta[dset.analysis["rundate"].strftime("%Y-%m-%d")]
        data = dset  # TODO: Another possibility: brdc.dset_raw

    sat_sys_definition = dict(
        G="GPS", R="GLONASS", E="Galileo", J="QZSS", C="BDS", I="IRNSS", S="SBAS Payload", M="Mixed"
    )
    rinex_version = "3.03"

    file_path = config.files.path("output_rinex3_nav", file_vars={**dset.vars, **dset.analysis})
    log.info(f"Write file {file_path}.")

    with config.files.open_path(file_path, create_dirs=True, mode="wt") as fid:

        #
        # Write RINEX navigation header
        #

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #      3.03           N: GNSS NAV DATA    E: GALILEO          RINEX VERSION / TYPE
        file_type = "N: GNSS NAV DATA"
        sat_sys = set(dset.system).pop() if len(set(dset.system)) == 1 else "M"
        fid.write(
            "{:>9s}{:11s}{:20s}{:20s}RINEX VERSION / TYPE\n"
            "".format(rinex_version, "", file_type, sat_sys + ": " + sat_sys_definition[sat_sys])
        )

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # CCRINEXN V1.6.0 UX  CDDIS               19990903 152236 UTC     PGM / RUN BY / DATE
        pgm = "where " + where.__version__
        run_by = util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else ""
        file_created = datetime.utcnow().strftime("%Y%m%d %H%M%S") + " UTC"
        fid.write("{:20s}{:20s}{:20s}PGM / RUN BY / DATE\n".format(pgm, run_by, file_created))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # IGS BROADCAST EPHEMERIS FILE                                COMMENT
        # TODO fid.write('{:60s}COMMENT\n'.format(line))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # BDSA    .1397E-07   .0000E+00  -.5960E-07   .5960E-07       IONOSPHERIC CORR
        # BDSB    .1106E+06  -.3277E+05  -.2621E+06   .1966E+06       IONOSPHERIC CORR
        if "iono_para" in meta:
            for type_, val in sorted(meta["iono_para"].items()):
                fid.write(
                    "{:4s} {:>12.4e}{:>12.4e}{:>12.4e}{:>12.4e} {:1s} {:2s}  IONOSPHERIC CORR\n"
                    "".format(
                        type_,
                        val["para"][0],
                        val["para"][1],
                        val["para"][2],
                        val["para"][3],
                        val["time_mark"],
                        val["sv_id"],
                    )
                )

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # BDUT -5.5879354477e-09-0.000000000e+00     14 1886          TIME SYSTEM CORR
        # GAUT  0.0000000000e+00 0.000000000e+00 172800 1886          TIME SYSTEM CORR
        if "time_sys_corr" in meta:
            for type_, val in sorted(meta["time_sys_corr"].items()):
                fid.write(
                    "{:4s} {:>17.10e}{:>16.9e}{:>7d}{:>5d}{:10s}TIME SYSTEM CORR\n"
                    "".format(type_, val["a0"], val["a1"], val["t"], val["w"], "")
                )

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #     16    17  1851     3                                    LEAP SECONDS
        if "leap_seconds" in meta:
            fid.write(
                "{:>6s}{:>6s}{:>6s}{:>6s}{:3s}{:33s}LEAP SECONDS\n"
                "".format(
                    meta["leap_seconds"]["leap_seconds"],
                    meta["leap_seconds"]["future_past_leap_seconds"],
                    meta["leap_seconds"]["week"],
                    meta["leap_seconds"]["week_day"],
                    meta["leap_seconds"]["time_sys"],
                    "",
                )
            )

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #                                                             END OF HEADER
        fid.write("{:60s}END OF HEADER\n".format(""))

        #
        # Write RINEX navigation data
        #

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # E11 2016 02 28 22 00 00  .654120231047E-04  .109707798401E-10  .000000000000E+00
        #       .400000000000E+01  .129375000000E+02  .319691887879E-08 -.292934515480E+01
        #       .460073351860E-06  .329698785208E-03  .683590769768E-05  .544061308098E+04
        #       .792000000000E+05  .391155481339E-07 -.125937621871E+01  .316649675369E-07
        #       .967916388734E+00  .197406250000E+03 -.653087089047E+00 -.571166648526E-08
        #       .276797244002E-09  .257000000000E+03  .188600000000E+04  .000000000000E+00
        #      -.100000000000E+01  .000000000000E+00 -.249128788710E-07 -.225845724344E-07
        #       .798850000000E+05  .000000000000E+00  .000000000000E+00  .000000000000E+00

        for idx in range(0, data.num_obs):

            # Remove observation epochs, which does not fit in the given time period
            # TODO: Is the time handling ok. Especially for BeiDou from day to day or week to week?
            rundate = datetime.combine(data.analysis["rundate"], dt_time.min)
            if data.time.gps.datetime[idx] < rundate or data.time.gps.datetime[idx] >= (rundate + timedelta(days=1)):
                continue

            # TODO:
            #        for drow in data.get_rows():
            #            fid.write('{d.sat:2d} {d.time:%Y%m%d %H%M%S} {d.inc0:13.4f}'.format(d=drow))
            #            fid.write('  {d.hurra:14.10f} ...'.format(d=drow))

            if data.system[idx] in ["R", "S"]:
                log.warning("Writing of RINEX navigation message is not implemented for GLONASS and SBAS satellites.")
                continue

            time = _time_system_correction(data, idx)
            gnss_data = _get_fields_based_on_system(data, idx)

            # BROADCAST ORBIT - 1
            fid.write(
                "{:3s} {:>4d} {:>2s} {:>2s} {:>2s} {:>2s} {:>2s}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(
                    data.satellite[idx],
                    time["toc"].year,
                    str(time["toc"].month).zfill(2),
                    str(time["toc"].day).zfill(2),
                    str(time["toc"].hour).zfill(2),
                    str(time["toc"].minute).zfill(2),
                    str(time["toc"].second).zfill(2),
                    data.sat_clock_bias[idx],
                    data.sat_clock_drift[idx],
                    data.sat_clock_drift_rate[idx],
                )
            )

            # BROADCAST ORBIT - 2
            fid.write(
                "    {:19.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.iode[idx], data.crs[idx], data.delta_n[idx], data.m0[idx])
            )

            # BROADCAST ORBIT - 3
            fid.write(
                "    {:19.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.cuc[idx], data.e[idx], data.cus[idx], data.sqrt_a[idx])
            )

            # BROADCAST ORBIT - 4
            fid.write(
                "    {:19.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(time["toe"], data.cic[idx], data.Omega[idx], data.cis[idx])
            )

            # BROADCAST ORBIT - 5
            fid.write(
                "    {:19.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.i0[idx], data.crc[idx], data.omega[idx], data.Omega_dot[idx])
            )

            # BROADCAST ORBIT - 6
            fid.write(
                "    {:19.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.idot[idx], gnss_data["data_info"], time["week"], gnss_data["l2p_flag"])
            )

            # BROADCAST ORBIT - 7
            fid.write(
                "    {:19.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(
                    data.sv_accuracy[idx], data.sv_health[idx], gnss_data["tgd_bgd"], gnss_data["iodc_groupdelay"]
                )
            )

            # BROADCAST ORBIT - 8
            fid.write(
                "    {:19.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(time["transmission_time"], gnss_data["interval"], 0.0, 0.0)
            )


def _time_system_correction(data, idx):
    """Apply correction to given time system for getting time scale defined in RINEX convention (:cite:`rinex3`)

    Following relationship are given between GNSS time scale (either BeiDou, Galileo, IRNSS or QZSS)
    :math:`t_{GNSS}` and GPS time scale :math:`t_{GPS}` (see Section 2.1.4 in :cite:`teunissen2017`):
    .. math::
          t_{GPS}  = t_{GNSS} + \\Delta t

    The time offset :math:`\\Delta t` is 0 s for Galileo, IRNSS and QZSS and for BeiDou 14 s. All these time scales
    are related to the International Atomic Time (TAI) by a certain time offset. In addition the GNSS week number
    is different depending on the GNSS. Galileo, IRNSS and QZSS are referring to the same GPS week, whereas BeiDou
    week starts at GPS week 1356.

    In this routine the given navigation epoch (time of clock (toc)), time of ephemeris (toe) and the broadcast
    ephemeris transmission time will be transformed to RINEX convention time scale (see :cite:`rinex3`) for example GPS
    or BeiDou time scale.
    """
    time = dict()
    # TODO hjegei: Introduction of BDS time scale in Time module
    if data.system[idx] == "C":
        time["toc"] = data.time.gps.datetime[idx] - timedelta(seconds=SYSTEM_TIME_OFFSET_TO_GPS_SECOND.get("C", 0))
    elif data.system[idx] in "EGIJ":
        time["toc"] = data.time.gps.datetime[idx]

    for key in ["toe", "transmission_time"]:
        if data.system[idx] == "C":
            time[key] = data[key].gps.gps_ws.seconds[idx] - SYSTEM_TIME_OFFSET_TO_GPS_SECOND.get("C", 0)
        elif data.system[idx] in "EGIJ":
            time[key] = data[key].gps.gps_ws.seconds[idx]

    if data.system[idx] == "C":
        time["week"] = data.gnss_week[idx] - SYSTEM_TIME_OFFSET_TO_GPS_WEEK.get("C", 0)
    elif data.system[idx] in "EGIJ":
        time["week"] = data.gnss_week[idx]

    # TODO: week rollover -> Should it be applied?

    return time


def _get_fields_based_on_system(data, idx):
    """Get general GNSS fields from GNSS specific ones

    Several fields in the RINEX navigation message have a different meaning depending on the GNSS. The GNSS specfic
    fields are converted to the general fields like 'data_info', 'interval', 'iodc_groupdelay', 'l2p_flag' and
    'tgd_bgd'. In the following table the relationship between the general fields and the GNSS dependent fields
    are shown:

    ===================== ==================== =================================================================
     General field         GNSS field           Description
    ===================== ==================== =================================================================
    data_info                                   Depending on GNSS this field has different meaning:
                           E: data_source        - Galileo: Data source information about the broadcast
                                                   ephemeris block, that means if the ephemeris block is based
                                                   on FNAV or INAV navigation message.
                           G: codes_l2           - GPS: Codes on L2 channel. Indication which codes are used
                                                   on L2 channel (P code, C/A code). See section 20.3.3.3.1.2
                                                   in :cite:`is-gps-200h`).
                           J: codes_l2           - QZSS: Codes on L2 channel. Indication if either C/A- or P-
                                                   code is used on L2 channel (0: spare, 1: P-code, 2: L1C/A
                                                   code). See section 4.1.2.7 in :cite:`is-qzss-pnt-001`.
                           C, I, R: None         - BeiDou, IRNSS and GLONASS: not used

    interval                                    Interval indicates either the curve-fit interval for GPS or QZSS
                                                ephemeris or for BeiDou the extrapolation interval for clock
                                                correction parameters:
                                                BeiDou to the clock correction parameters:
                           C: age_of_clock_corr   - BeiDou: Age of data, clock (AODC) is the extrapolated interval
                                                    of clock correction parameters. It indicates the time
                                                    difference between the reference epoch of clock correction
                                                    parameters and the last observation epoch for extrapolating
                                                    clock correction parameters. Meaning of AODC
                                                     < 25  Age of the satellite clock correction parameters in hours
                                                       25  Age of the satellite clock correction parameters is two days
                                                       ...
                                                    See section 5.2.4.9 in :cite:`bds-sis-icd`.
                           G: fit_interval        - GPS: Indicates the curve-fit interval used by the GPS Control
                                                    Segment in determining the ephemeris parameters, which is
                                                    given in HOURS (see section 6.6 in :cite:`rinex2`).
                           J: fit_interval        - QZSS: Fit interval is given as flag (see section 4.1.2.7 in
                                                    :cite:`is-qzss-pnt-001`)
                                                          0 - 2 hours
                                                          1 - more than 2 hours
                                                          blank - not known
                           E, I, R: None          - Galileo, IRNSS and GLONASS: not used

    iodc_groupdelay                             Clock issue of data or group delay depending on GNSS:
                           C: tgd_b2_b3           - BeiDou: B2/B3 TGD2
                           E: bgd_e1_e5b          - Galileo: E1-E5b BGD (see section 5.1.5 in 
                                                    :cite:`galileo-os-sis-icd`)
                           G: iodc                - GPS: IODC (Clock issue of data indicates changes
                                                    (set equal to IODE))
                           J: iodc                - QZSS: IODC
                           I, R: None             - IRNSS and GLONASS: not used

    l2p_flag                                    L2 P-code data flag is only used by GPS and QZSS:
                           G: l2p_flag            - GPS: When bit 1 of word four is a "1", it shall indicate that
                                                    the NAV data stream was commanded OFF on the P-code of the L2
                                                    channel (see section 20.3.3.3.1.6 in :cite:`is-gps-200h`).
                           J: l2p_flag            - QZSS: L2P data flag set to 1 since QZSS does not track L2P.
                                                    See section 4.1.2.7 in :cite:`is-qzss-pnt-001`.
                           C, E, I, R: None       - BeiDou, Galileo, IRNSS and GLONASS: not used

    tgd_bgd                                     Total group delay (TGD) or broadcast group delay (BGD) for
                                                Galileo:
                           C: tgd_b1_b3           - BeiDou: B1/B3 TGD1
                           E: bgd_e1_e5a          - Galileo: E1-E5a BGD (see section 5.1.5 in 
                                                    :cite:`galileo-os-sis-icd`)
                           G: tgd                 - GPS: TGD (:math:`L_1 - L_2` delay correction term. See
                                                    section 20.3.3.3.3.2 in :cite:`is-gps-200h`.)
                           I: tgd                 - IRNSS: TGD
                           J: tgd                 - QZSS: TGD
                           R: None                - GLONASS: not used

    Args:
        data (Dataset):   Where Dataset with navigation information
        idx (int):        Index of current navigation epoch


    Returns:
        dict:   Dictionary with general fields valid for given navigation epoch
    """
    gnss_data = dict()

    if data.system[idx] == "C":
        gnss_data["data_info"] = 0
        gnss_data["interval"] = data.age_of_clock_corr[idx]
        gnss_data["iodc_groupdelay"] = data.tgd_b2_b3[idx]
        gnss_data["l2p_flag"] = 0
        gnss_data["tgd_bgd"] = data.tgd_b1_b3[idx]

    elif data.system[idx] == "E":
        gnss_data["data_info"] = data.data_source[idx]
        gnss_data["interval"] = 0
        gnss_data["iodc_groupdelay"] = data.bgd_e1_e5b[idx]
        gnss_data["l2p_flag"] = 0
        gnss_data["tgd_bgd"] = data.bgd_e1_e5a[idx]

    elif data.system[idx] in "GJ":
        gnss_data["data_info"] = data.codes_l2[idx]
        gnss_data["interval"] = data.fit_interval[idx]
        gnss_data["iodc_groupdelay"] = data.iodc[idx]
        gnss_data["l2p_flag"] = data.l2p_flag[idx]
        gnss_data["tgd_bgd"] = data.tgd[idx]

    elif data.system[idx] == "I":
        gnss_data["data_info"] = 0
        gnss_data["interval"] = 0
        gnss_data["iodc_groupdelay"] = 0
        gnss_data["l2p_flag"] = 0
        gnss_data["tgd_bgd"] = data.tgd[idx]

    elif data.system[idx] == "R":
        gnss_data["data_info"] = 0
        gnss_data["interval"] = 0
        gnss_data["iodc_groupdelay"] = 0
        gnss_data["l2p_flag"] = 0
        gnss_data["tgd_bgd"] = 0

    return gnss_data
