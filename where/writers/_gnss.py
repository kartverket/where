"""Common functions used for plotting GNSS analysis results

Description:
------------
TODO

"""
# Standard liberay imports
from datetime import datetime
from typing import List

# External liberary imports
import numpy as np


# Where imports
from where import apriori


def get_grc_csv_header() -> List[str]:
    """Get header line of GRC CSV format

    Returns:
        Header line of GRC CSV format
    """
    return [
            "Service Line",
            "Service Category",
            "Business Service",
            "Batch",
            "Satellite",
            "PRN",
            "Slot",
            "GSS Site",
            "Station",
            "Service",
            "Type",
            "Mode",
            "Target",
            "Unit",
            "Month/Year",
            "Result",
    ]


def get_grc_csv_row(
            kpi: str, 
            mode: str, 
            date: str, 
            result: float, 
            station: str="",
            satellite: str="",
) -> List[str]:
    """Create row in GRC CSV format in dependency of given arguments

    Args:
        kpi:        KPI name (e.g. hpe, vpe)
        mode:       Signal combination (e.g. e1, e1e5b)
        date:       Date in format month-year (e.g. July-2021)
        result:     KPI value
        station:    Station name
        satellite:  Satellite name

    Returns:
        Line in GRC CSV format
    """
    batch = ""
    prn = ""
    slot = ""

    batch_def = { 
            "GALILEO-1": "IOV",
            "GALILEO-2": "FOC",
    }

    business_def = {
            "hpe": "FOM-OS-0015 Horizontal Positioning Service Accuracy per Station over Month",
            "site_vel_3d": "FOM-OS-0021 3D Velocity Service Accuracy per Station over Month",
            "sisre": "SDD-OS-0012 SIS Ranging Accuracy over All Satellites over Month",
            "sisre_sat": "SDD-OS-0013 SIS Ranging Accuracy per Satellite over Month",
            "vpe": "FOM-OS-0016 Vertical Positioning Service Accuracy per Station over Month",
    }

    category_def = {
            "hpe": "Position Domain",
            "site_vel_3d": "Position Domain",
            "sisre": "Ranging Domain",
            "sisre_sat": "Ranging Domain",
            "vpe": "Position Domain",
    } 

    # TODO: Better solution? dset.meta["obstypes"] handling has to be improved.
    mode_def = {
            "e1": "E1",
            "e1e5a": "E1/E5a",
            "e1e5b": "E1/E5b",
            "l1":    "L1",
            "l1l2":  "L1/L2",
    }

    slot_def = {
            "GSAT101": "B05",
            "GSAT102": "B06",
            "GSAT103": "C04",
            "GSAT104": "C05",
            "GSAT201": "EXT01",
            "GSAT202": "EXT02",
            "GSAT203": "B08",
            "GSAT204": "B03",
            "GSAT205": "A08",
            "GSAT206": "A05",
            "GSAT207": "C06",
            "GSAT208": "C07",
            "GSAT209": "C02",
            "GSAT210": "A02",
            "GSAT211": "A06",
            "GSAT212": "C08",
            "GSAT213": "C03",
            "GSAT214": "C01",
            "GSAT215": "A03",
            "GSAT216": "A07",
            "GSAT217": "A04",
            "GSAT218": "A01",
            "GSAT219": "B04",
            "GSAT220": "B01",
            "GSAT221": "B02",
            "GSAT222": "B07",
    }

    station_def = {
            "brux": "Brussels (Belgium)",
            "cpvg": "Cap-Vert (Cabo Verde)",
            "koug": "Kourou (French Guiana)",
            "hofs": "Hoefn (Iceland)",
            "hons": "Honningsvag (Norway)",
            "krss": "Kristiansand (Norway)",
            "mas1": "Maspalomas (Spain)",
            "nabd": "Ny Alesund (Norway)",
            "nklg": "N' Koltang (Gabon)",
            "vegs": "Vega (Norway)",
    }

    target_def = {
            "hpe": "",
            "site_vel_3d": "",
            "sisre": 2,
            "sisre_sat": "7",
            "vpe": "",
    }

    type_ = "Single" if len(mode) == 2 else "Dual" # TODO: Better solution?

    unit_def = {
            "hpe": "m",
            "site_vel_3d": "m/s",
            "sisre": "m",
            "sisre_sat": "m",
            "vpe": "m",
    }

    mode_to_write = mode_def[mode] if mode in mode_def.keys() else mode
    station_to_write = station_def[station] if station in station_def.keys() else ""

    if satellite:
        atx = apriori.get("gnss_antenna_correction")
        used_date = datetime.strptime(f"{date}-01", "%y-%b-%d")
        sat_info = atx.get_satellite_info(satellite, used_date)

        batch = batch_def[sat_info["sat_type"]] if sat_info["sat_type"] in batch_def.keys() else ""
        prn = sat_info["sat_code"].replace("E", "GSAT") if sat_info["sat_code"].startswith("E") else sat_info["sat_code"]
        slot = slot_def[prn] if prn in slot_def.keys() else ""


    return [
        "Open Service", # Service Line
        category_def[kpi], # Service Category
        business_def[kpi], # Business Service
        batch, # Batch
        satellite, # Satellite
        prn, # PRN
        slot, # Slot
        station_to_write, # GSS Site
        station.upper(), # Station
        "OS", # Service
        type_, # Type
        mode_to_write, # Mode
        target_def[kpi], # Target
        unit_def[kpi], # Unit
        date, # Month/Year
        f"{result:.5f}", # Result
    ]


