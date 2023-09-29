"""Write RINEX observations  in Rinex format 3.03

Description:
------------
Write data in the RINEX observations file format (see :cite:`rinex3`).


"""

# Standard library imports
from datetime import datetime
from math import isnan

# Midgard imports
from midgard.dev import plugins

# Where imports
import where
from where.lib import config
from where.lib import log


@plugins.register
def rinex3_obs(dset):
    """Write RINEX observations in Rinex format 3.03

    Args:
        dset:       Dataset, a dataset containing the data.
    """

    # Initialze variables
    meta = dset.meta
    version = "3.03"
    program = "Where v{}".format(where.__version__)
    run_by = "NMA"
    date = datetime.utcnow()
    time_sys = "GPS"  # TODO: So far only GPS time system can be handled by Where.
    file_created = "{:15s} {:3s}".format(date.strftime("%Y%m%d %H%M%S"), "UTC")
    pos_x = dset.site_pos.trs.x[0]
    pos_y = dset.site_pos.trs.y[0]
    pos_z = dset.site_pos.trs.z[0]

    cfg_sampling_rate = config.tech.sampling_rate.float
    num_satellites = len(dset.unique("satellite"))

    if meta["file_type"] == "O":
        file_type = "OBSERVATION DATA"

    if meta["interval"] <= float(cfg_sampling_rate):
        sampling_rate = cfg_sampling_rate
    else:
        sampling_rate = meta["interval"]
    dset.vars["sampling_rate"] = str(int(sampling_rate))  # Used as placeholder for determination of output file name

    with config.files.open("output_rinex3_obs", file_vars=dset.vars, mode="wt") as fid:

        # ================================
        #  Write RINEX observation header
        # ================================
        #
        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #      3.02           OBSERVATION DATA    M (MIXED)           RINEX VERSION / TYPE
        fid.write("{:>9s}{:11s}{:20s}{:20s}RINEX VERSION / TYPE\n".format(version, "", file_type, meta["sat_sys"]))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # MAKERINEX 2.0.20023 BKG/GOWETTZELL      2016-03-02 00:20    PGM / RUN BY / DATE
        fid.write("{:20s}{:20s}{:20s}PGM / RUN BY / DATE\n".format(program, run_by, file_created))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # G = GPS R = GLONASS E = GALILEO S = GEO M = MIXED           COMMENT

        # Add comments related applying HAS code bias corrections
        if "has_corrected_obstypes" in meta:
            if not "comment" in meta.keys():
                meta["comment"] = list()
                
            meta["comment"].extend([
                # ----+----1----+----2----+----3----+----4----+----5----+----6
                  "",
                  "HAS corrections are applied for code observations by Where",
                  "software. The corrected observations are listed for a given",
                  "GNSS and used HAS code bias signal (after convention",
                  "<sys> <sig>: <obs_codes>):",
            ])
            for sys in sorted(meta["has_corrected_obstypes"].keys()):
                for signal in sorted(meta["has_corrected_obstypes"][sys].keys()):
                    line = f"{sys} {signal+':':6s} {' '.join(meta['has_corrected_obstypes'][sys][signal])}"
                    meta["comment"].append(line)

        if "comment" in meta:
            for line in meta["comment"]:
                fid.write("{:60s}COMMENT\n".format(line))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # stas                                                        MARKER NAME
        fid.write("{:60s}MARKER NAME\n".format(meta["marker_name"]))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # 66008M005                                                   MARKER NUMBER
        if "marker_number" in meta:
            fid.write("{:60s}MARKER NUMBER\n".format(meta["marker_number"]))

        if "marker_type" in meta:
            fid.write("{:60s}MARKER TYPE\n".format(meta["marker_type"]))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # SATREF              Norwegian Mapping Authority             OBSERVER / AGENCY
        fid.write("{:20s}{:40s}OBSERVER / AGENCY\n".format(meta["observer"], meta["agency"]))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # 3008040             SEPT POLARX4        2.9.0               REC # / TYPE / VERS
        fid.write(
            "{:20s}{:20s}{:20s}REC # / TYPE / VERS\n"
            "".format(meta["receiver_number"], meta["receiver_type"], meta["receiver_version"])
        )

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # CR620012101         ASH701945C_M    SCIS                    ANT # / TYPE
        fid.write("{:20s}{:40s}ANT # / TYPE\n".format(meta["antenna_number"], meta["antenna_type"]))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #   3275756.7623   321111.1395  5445046.6477                  APPROX POSITION XYZ
        fid.write("{:>14.4f}{:>14.4f}{:>14.4f}{:18s}APPROX POSITION XYZ\n".format(pos_x, pos_y, pos_z, ""))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #         0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N
        fid.write(
            "{:>14.4f}{:>14.4f}{:>14.4f}{:18s}ANTENNA: DELTA H/E/N\n"
            "".format(meta["antenna_height"], meta["antenna_east"], meta["antenna_north"], "")
        )

        if "ant_vehicle_x" in meta:
            fid.write(
                "{:>14.4f}{:>14.4f}{:>14.4f}{:18s}ANTENNA: DELTA X/Y/Z\n"
                "".format(meta["ant_vehicle_x"], meta["ant_vehicle_y"], meta["ant_vehicle_z"], "")
            )

        # TODO: ANTENNA:PHASECENTER
        # TODO: ANTENNA:B.SIGHT XYZ
        # TODO: ANTENNA:ZERODIR AZI
        # TODO: ANTENNA:ZERODIR XYZ
        # TODO: CENTER OF MASS: XYZ

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # G   26 C1C C1P L1C L1P D1C D1P S1C S1P C2P C2W C2S C2L C2X  SYS / # / OBS TYPES
        #        L2P L2W L2S L2L L2X D2P D2W D2S D2L D2X S2P S2W S2S  SYS / # / OBS TYPES
        # R   16 C1C C1P L1C L1P D1C D1P S1C S1P C2C C2P L2C L2P D2C  SYS / # / OBS TYPES
        #        D2P S2C S2P                                          SYS / # / OBS TYPES
        for sys in sorted(meta["obstypes"]):
            obstypes = meta["obstypes"][sys].copy()
            num_lines = int(len(obstypes) / 13) + 1
            for line in range(0, num_lines):
                num_obstypes = len(obstypes)
                num_obstypes_str = str(num_obstypes) if line == 0 else ""
                spaces = "  " if meta["version"].startswith("2") else " "
                if num_obstypes <= 13:
                    fid.write(
                        "{:1s}{:>5s} {:53s}SYS / # / OBS TYPES\n".format(sys, num_obstypes_str, spaces.join(obstypes))
                    )
                    if num_obstypes == 13:
                        break
                else:
                    fid.write(
                        "{:1s}{:>5s} {:53s}SYS / # / OBS TYPES\n".format(
                            sys, num_obstypes_str, spaces.join(obstypes[0:13])
                        )
                    )
                    del obstypes[0:13]

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # DBHZ                                                        SIGNAL STRENGTH UNIT
        if "signal_strength_unit" in meta:
            fid.write("{:60s}SIGNAL STRENGTH UNIT\n".format(meta["signal_strength_unit"]))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #      1.000                                                  INTERVAL
        if "interval" in meta:
            fid.write("{:>10.3f}{:50s}INTERVAL\n".format(sampling_rate, ""))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #   2016    03    01    00    00   00.0000000     GPS         TIME OF FIRST OBS
        if not meta["time_sys"] == "GPS":
            log.fatal(f"Time system {meta['time_sys']!r} is not implemented so far in Where")
        d = dset.time.gps.datetime[0]
        fid.write(
            "{:>6d}{:>6d}{:>6d}{:>6d}{:>6d}{:>13.7f}{:>8s}{:9s}TIME OF FIRST OBS\n"
            "".format(d.year, d.month, d.day, d.hour, d.minute, d.second, time_sys, "")
        )

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #   2016    03    01    23    59   59.0000000     GPS         TIME OF LAST OBS
        if "time_last_obs" in meta:
            d = dset.time.gps.datetime[-1]
            fid.write(
                "{:>6d}{:>6d}{:>6d}{:>6d}{:>6d}{:>13.7f}{:>8s}{:9s}TIME OF LAST OBS\n"
                "".format(d.year, d.month, d.day, d.hour, d.minute, d.second, time_sys, "")
            )

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #      0                                                      RCV CLOCK OFFS APPL
        if "rcv_clk_offset_flag" in meta:
            fid.write("{:>6s}{:54s}RCV CLOCK OFFS APPL\n".format(meta["rcv_clk_offset_flag"], ""))

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # G APPL_DCB          xyz.uvw.abc//pub/dcb_gps.dat            SYS / DCBS APPLIED
        if "dcbs_applied" in meta:
            for sys in sorted(meta["dcbs_applied"]):
                if sys in meta["obstypes"]:
                    fid.write(
                        "{:1s} {:17s} {:40s}SYS / DCBS APPLIED\n"
                        "".format(sys, meta["dcbs_applied"][sys]["prg"], meta["dcbs_applied"][sys]["url"])
                    )

        if "pcvs_applied" in meta:
            for sys in sorted(meta["pcvs_applied"]):
                if sys in meta["obstypes"]:
                    fid.write(
                        "{:1s} {:17s} {:40s}SYS / PCVS APPLIED\n"
                        "".format(sys, meta["pcvs_applied"][sys]["prg"], meta["pcvs_applied"][sys]["url"])
                    )
        # TODO: SYS / SCALE FACTOR

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # G L1C  0.00000  12 G01 G02 G03 G04 G05 G06 G07 G08 G09 G10  SYS / PHASE SHIFT
        #                    G11 G12                                  SYS / PHASE SHIFT
        # G L1W  0.00000                                              SYS / PHASE SHIFT
        if "phase_shift" in meta:
            num_sat_limit = 10
            for sys, obstypes in sorted(meta["phase_shift"].items()):

                if sys not in meta["obstypes"]:
                    continue

                if not obstypes:
                    # Note: Phase corrections are unknown.
                    fid.write("{:1s}{:59s}SYS / PHASE SHIFT\n".format(sys, ""))
                    continue

                for type_ in obstypes:
                    if type_ in meta["obstypes"][sys]:
                        # TODO: Remove unused satellites
                        sats = meta["phase_shift"][sys][type_]["sat"].copy()
                        num_lines = int(len(sats) / num_sat_limit) + 1
                        for line in range(0, num_lines):
                            num_sats = len(sats)
                            if line == 0:
                                num_sats_str = str(num_sats) if num_sats > 0 else ""
                                correction = meta["phase_shift"][sys][type_]["corr"]
                                correction = f"{float(correction):>8.5f}" if correction else f"{correction:>8s}"
                                phase_shift_str = "{:1s} {:>3s} {}{:>4s}" "".format(
                                    sys, 
                                    type_, 
                                    correction, 
                                    num_sats_str,
                                )
                            else:
                                phase_shift_str = ""

                            if num_sats <= num_sat_limit:
                                fid.write("{:18s} {:41s}SYS / PHASE SHIFT\n".format(phase_shift_str, " ".join(sats)))
                            else:
                                fid.write(
                                    "{:18s} {:41s}SYS / PHASE SHIFT\n".format(
                                        phase_shift_str, " ".join(sats[0:num_sat_limit])
                                    )
                                )
                                del sats[0:num_sat_limit]

        # TODO: WAVELENGTH FACT L1/2  -> given only for RINEX 2.11, but could be of interest in RINEX file

        if "R" in meta["obstypes"]:
            # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
            #  22 R01  1 R02 -4 R03  5 R04  6 R05  1 R06 -4 R07  5 R08  6 GLONASS SLOT / FRQ #
            #     R09 -6 R10 -7 R11  0 R13 -2 R14 -7 R15  0 R17  4 R18 -3 GLONASS SLOT / FRQ #
            #     R19  3 R20  2 R21  4 R22 -3 R23  3 R24  2               GLONASS SLOT / FRQ #
            # TODO: Remove unused satellites from 'GLONASS SLOT / FRQ #'
            if "glonass_slot" in meta:
                num_sat = len(meta["glonass_slot"])
                glonass_slots = dict(meta["glonass_slot"])
                num_lines = int(num_sat / 8) + 1
                for idx in range(0, num_lines):
                    line = "{:>3d}".format(num_sat) if idx == 0 else "   "
                    for num, (slot, bias) in enumerate(sorted(glonass_slots.items())):
                        if num == 8:
                            break
                        line = line + " {:3s} {:>2d}".format(slot, bias)
                        del glonass_slots[slot]
                    fid.write(line.ljust(60) + "GLONASS SLOT / FRQ #\n")

            # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
            #  C1C  -10.000 C1P  -10.123 C2C  -10.432 C2P  -10.634        GLONASS COD/PHS/BIS
            line = ""
            if "glonass_bias" in meta:
                for type_, bias in sorted(meta["glonass_bias"].items()):
                    if type_ in meta["obstypes"]["R"]:
                        line = line + " {:3s} {:8.3f}".format(type_, float(bias))
            fid.write(line.ljust(60) + "GLONASS COD/PHS/BIS\n")

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #     16    17  1851     3                                    LEAP SECONDS
        #
        # NOTE: Entries 'future_past_leap_seconds', 'week', 'week_day' and 'time_sys' are not given in RINEX version
        #       2.11.
        if "leap_seconds" in meta:
            if meta["version"].startswith("2"):
                fid.write("{:>6d}{:54s}LEAP SECONDS\n".format(int(meta["leap_seconds"]["leap_seconds"]), ""))
            else:
                fid.write(
                    "{:>6d}{:>6s}{:>6s}{:>6s}{:3s}{:33s}LEAP SECONDS\n"
                    "".format(
                        int(meta["leap_seconds"]["leap_seconds"]),
                        meta["leap_seconds"]["future_past_leap_seconds"],
                        meta["leap_seconds"]["week"],
                        meta["leap_seconds"]["week_day"],
                        meta["leap_seconds"]["time_sys"],
                        "",
                    )
                )

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #     71                                                      # OF SATELLITES
        fid.write("{:>6d}{:54s}# OF SATELLITES\n".format(num_satellites, ""))
        # TODO: PRN / # OF OBS

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #                                                             END OF HEADER
        fid.write("{:60s}END OF HEADER\n".format(""))

        # ================================
        #  Write RINEX observation data
        # ================================
        #
        epoch_prev = dset.time.gps.datetime[0]
        first_obs_in_epoch = True
        obs_epoch_cache = dict()

        # Loop over all observations
        for idx in range(0, dset.num_obs):

            # Write epoch (reading of observations from epoch 'd_prev' in 'obs_epoch_cache' is finished)
            epoch = dset.time.gps.datetime[idx]
            if epoch_prev != epoch:
                num_sat = idx - idx_epoch_start  # TODO: idx_epoch_start is not defined
                _write_epoch(dset, fid, obs_epoch_cache, idx, num_sat, epoch_prev)
                first_obs_in_epoch = True

            if first_obs_in_epoch is True:
                obs_epoch_cache = dict()
                idx_epoch_start = idx
                first_obs_in_epoch = False

            # Save observations for a given epoch in obs_epoch_cache
            #
            # NOTE: The caching is mainly necessary to determine the number of satellites for an epoch and to be
            #       flexible in what kind of order the observation types should be written. The order of the
            #       observation types for a given GNSS is defined via dset.meta['obstypes'] variable.
            if dset.satellite[idx] in obs_epoch_cache:
                log.fatal(f"Satellite {dset.satellite[idx]} occurs twice in epoch {dset.time.gps.datetime[idx]}")

            for type_ in dset.meta["obstypes"][dset.system[idx]]:
                obs = " " if isnan(dset.obs[type_][idx]) else f"{dset.obs[type_][idx]:>14.3f}"
                lli = " " if isnan(dset.lli[type_][idx]) else str(int(dset.lli[type_][idx]))
                snr = " " if isnan(dset.snr[type_][idx]) else str(int(dset.snr[type_][idx]))
                obs_epoch_cache.setdefault(dset.satellite[idx], list()).append({"obs": obs, "lli": lli, "snr": snr})
            epoch_prev = epoch

        # Write last epoch
        num_sat = (idx + 1) - idx_epoch_start
        _write_epoch(dset, fid, obs_epoch_cache, idx, num_sat, epoch_prev)


def _write_epoch(dset, fid, obs_epoch_cache, idx, num_sat, epoch):
    """Write current RINEX observation file epoch

    Args:
        dset:            Dataset, a dataset containing the data.
        fid:             File object
        obs_epoch_cache: Cache with observations for given epoch
        idx:             Current observation index
        num_sat:         Number of satellites
        epoch:           Current epoch given as Datetime object
    """
    # Write epoch record
    if "rcv_clk_offset_flag" in dset.meta:
        if dset.meta["rcv_clk_offset_flag"] == "0":
            rcv_clk_offset = "{:15s}".format("")  # Blank if receiver clock offset is not given.
    else:
        rcv_clk_offset = "" if isnan(dset.rcv_clk_offset[idx]) else "{:>15.12f}".format(dset.rcv_clk_offset[idx])
    fid.write(
        "> {:3d}{:>3d}{:>3d}{:>3d}{:>3d}{:>11.7f}{:>3d}{:3d}{:6s}{:15s}\n"
        "".format(
            epoch.year,
            epoch.month,
            epoch.day,
            epoch.hour,
            epoch.minute,
            epoch.second,
            int(dset.epoch_flag[idx]),
            num_sat,
            "",
            rcv_clk_offset,
        )
    )

    # Write observation record
    for sat in sorted(obs_epoch_cache):
        fid.write("{:3s}".format(sat))
        for kk in range(0, len(obs_epoch_cache[sat])):
            if obs_epoch_cache[sat][kk]["obs"] == 0.0:
                fid.write("{:>16s}".format(""))  # Write blank if no observation is given.
            else:
                fid.write(
                    "{:>14s}{:1s}{:1s}".format(
                        obs_epoch_cache[sat][kk]["obs"],
                        obs_epoch_cache[sat][kk]["lli"],
                        obs_epoch_cache[sat][kk]["snr"],
                    )
                )
        fid.write("\n")
