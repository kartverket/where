"""A parser for reading data from igs.snx file based on IGS sitelog files in SINEX format

Example:
--------

    from where import parsers
    p = parsers.parse_file(parser_name='gnss_sinex_igs', file_path='igs.snx')
    data = p.as_dict()

Description:
------------

Reads station information (e.g. approximated station coordinates, receiver and antenna type, station eccentricities,
...) igs.snx file in SINEX format.


"""
# Standard library imports
from datetime import datetime

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_sinex import SinexParser


@plugins.register
class IgsSnxParser(SinexParser):
    """A parser for reading data from igs.snx file based on IGS sitelog files in SINEX format

            site - Site dictionary, whereby keys are the site identifiers and values are a site entry
                   dictionary with the keys 'site_antenna', 'site_eccentricity', 'site_id' and 'site_receiver'. The
                   site dictionary has following strucuture:

                      self.site[site] = { 'site_antenna':          [],  # SITE/ANTENNA SINEX block information
                                          'site_eccentricity':     [],  # SITE/ECCENTRICITY block information
                                          'site_id':               {},  # SITE/ID block information
                                          'site_receiver':         [],  # SITE/RECEIVER block information }

                   with the site entry dictionary entries

                      site_antenna[ii]      = { 'point_code':         point_code,
                                                'soln':               soln,
                                                'obs_code':           obs_code,
                                                'start_time':         start_time,
                                                'end_time':           end_time,
                                                'antenna_type':       antenna_type,
                                                'radome_type':        radome_type,
                                                'serial_number':      serial_number }

                      site_eccentricity[ii] = { 'point_code':         point_code,
                                                'soln':               soln,
                                                'obs_code':           obs_code,
                                                'start_time':         start_time,
                                                'end_time':           end_time,
                                                'reference_system':   reference_system,
                                                'vector_1':           vector_1,
                                                'vector_2':           vector_2,
                                                'vector_3':           vector_3,
                                                'vector_type':        UNE }

                      site_id               = { 'point_code':         point_code,
                                                'domes':              domes,
                                                'marker':             marker,
                                                'obs_code':           obs_code,
                                                'description':        description,
                                                'approx_lon':         approx_lon,
                                                'approx_lat':         approx_lat,
                                                'approx_height':      approx_height }

                      site_receiver[ii]     = { 'point_code':         point_code,
                                                'soln':               soln,
                                                'obs_code':           obs_code,
                                                'start_time':         start_time,
                                                'end_time':           end_time,
                                                'receiver_type':      receiver_type,
                                                'serial_number':      serial_number,
                                                'firmware':           firmware }

                   The counter 'ii' ranges from 0 to n and depends on how many antenna type, receiver type and
                   antenna monument changes were done at each site. Note also, that time entries (e.g. start_time,
                   end_time) are given in Modified Julian Date. If the time is defined as 00:000:00000 in the SINEX
                   file, then the value is saved as 'None' in the Sinex class.

    """

    def setup_parser(self):
        return (
            self.site_id,
            self.site_receiver,
            self.site_antenna,
            self.site_gps_phase_center,
            self.site_eccentricity,
            self.satellite_id,
            self.satellite_phase_center,
        )

    def parse_site_id(self, data):
        """Parse SITE/ID SINEX block
        """
        for d in data:
            site_key = d["site_code"]
            self.data.setdefault(site_key, dict())
            self.data[site_key].setdefault("site_id", dict())
            self.data[site_key]["site_id"] = {n: d[n] for n in d.dtype.names}

    def parse_site_antenna(self, data):
        """Parse SITE/ANTENNA SINEX block
        """
        for d in data:
            site_key = d["site_code"]
            # TODO_hjegei: How to remove d['site_code'] from d?
            add_dict = {n: d[n] for n in d.dtype.names}  # Generate dictionary with all SINEX field entries
            add_dict["antenna_type"], add_dict["radome_type"] = d["antenna_type"].split()
            self.data.setdefault(site_key, dict())
            self.data[site_key].setdefault("site_antenna", list())
            self.data[site_key]["site_antenna"].append(add_dict)

    def parse_site_receiver(self, data):
        """Parse SITE/RECEIVER SINEX block
        """
        for d in data:
            site_key = d["site_code"]
            # TODO_hjegei: How to remove d['site_code'] from d?
            self.data.setdefault(site_key, dict())
            self.data[site_key].setdefault("site_receiver", list())
            self.data[site_key]["site_receiver"].append({n: d[n] for n in d.dtype.names})

    def parse_site_eccentricity(self, data):
        """Parse SITE/ECCENTRICITY SINEX block
        """
        for d in data:
            site_key = d["site_code"]
            # TODO_hjegei: How to remove d['site_code'] from d?
            self.data.setdefault(site_key, dict())
            self.data[site_key].setdefault("site_eccentricity", list())
            self.data[site_key]["site_eccentricity"].append({n: d[n] for n in d.dtype.names})

    # TODO: Improve handling of SITE/GPS_PHASE_CENTER, SATELLITE/ID and SATELLITE/PHASE_CENTER
    #       SINEX block.
