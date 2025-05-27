"""Common functions used for writing of GNSS pipeline analysis results

Description:
------------
TODO

"""
# Standard liberay imports
from datetime import datetime
from typing import List

# Where imports
from where import apriori


def get_grc_csv_header() -> List[str]:
    """Get header line of GRC CSV format

    Returns:
        Header line of GRC CSV format
    """
    return [
            "Constellation",
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
            "Date",
            "Result",
    ]


def get_grc_csv_row(
            constellation: str, 
            kpi: str, 
            mode: str, 
            date: str, 
            result: float, 
            station: str="",
            satellite: str="",
) -> List[str]:
    """Create row in GRC CSV format in dependency of given arguments

    Args:
        constellation:  Constellation name (e.g. Galileo, GPS)
        kpi:            KPI name (e.g. hpe, vpe)
        mode:           Signal combination (e.g. e1, e1e5b)
        date:           Date in format year-month (e.g. 2021-Jul)
        result:         KPI value
        station:        Station name
        satellite:      Satellite name

    Returns:
        Line in GRC CSV format as list
    """
    batch = ""
    svn = ""
    slot = ""

    batch_def = { 
            "GALILEO-1": "IOV",
            "GALILEO-2": "FOC",
            "BLOCK IIIA": "III",
            "BLOCK IIF": "IIF",
            "BLOCK IIR-A": "IIR",
            "BLOCK IIR-B": "IIR",
            "BLOCK IIR-M": "IIR-M",
    }

    business_def = {
            "Galileo": {
                "hpe": "FOM-OS-07: Horizontal Positioning Service Accuracy per Station over Month",
                "site_vel_3d": "FOM-OS-11: 3D Velocity Service Accuracy per Station over Month",
                "sisre": "SDD-OS-09: SIS Ranging Accuracy at GA over All Satellites over Month",
                "sisre_sat": "SDD-OS-10: SIS Ranging Accuracy at GA per Satellite over Month",
                "vpe": "FOM-OS-08: Vertical Positioning Service Accuracy per Station over Month",
            },
            "GPS": {
                "hpe": "FOM-OS-07: Horizontal Positioning Service Accuracy per Station over Month",
                "site_vel_3d": "FOM-OS-11: 3D Velocity Service Accuracy per Station over Month",
                "sisre": "SPS-OS-02: SIS Ranging Accuracy over All Satellites over Month",
                "sisre_sat": "SPS-OS-03: SIS Ranging Accuracy per Satellite over Month",
                "vpe": "FOM-OS-08: Vertical Positioning Service Accuracy per Station over Month",
            },
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
            "SVN41": "B03",
            "SVN43": "F06",
            "SVN44": "B03",
            "SVN45": "D03",
            "SVN48": "A04",
            "SVN50": "E03",
            "SVN51": "E04",
            "SVN52": "A02",
            "SVN53": "C04",
            "SVN55": "F02",
            "SVN56": "B01",
            "SVN57": "C01",
            "SVN58": "B04",
            "SVN59": "C05",
            "SVN61": "D01",
            "SVN62": "B02",
            "SVN63": "D02",
            "SVN64": "A03",
            "SVN65": "A01",
            "SVN66": "C02",
            "SVN67": "D04",
            "SVN68": "F03",
            "SVN69": "E01",
            "SVN70": "F01",
            "SVN71": "B05",
            "SVN72": "C03",
            "SVN73": "E02",
            "SVN74": "F04",
            "SVN75": "D06",
            "SVN76": "E05",
            "SVN77": "B06",
            "SVN78": "D05",
            "SVN79": "A06",
            "SVN80": "D02",
            "GSAT0101": "B05",
            "GSAT0102": "B06",
            "GSAT0103": "C04",
            "GSAT0104": "C14",
            "GSAT0201": "EXT01",
            "GSAT0202": "EXT02",
            "GSAT0203": "B08",
            "GSAT0204": "B14",
            "GSAT0205": "A08",
            "GSAT0206": "A05",
            "GSAT0207": "C06",
            "GSAT0208": "C07",
            "GSAT0209": "C02",
            "GSAT0210": "A02",
            "GSAT0211": "A06",
            "GSAT0212": "C08",
            "GSAT0213": "C03",
            "GSAT0214": "C01",
            "GSAT0215": "A03",
            "GSAT0216": "A07",
            "GSAT0217": "A04",
            "GSAT0218": "A01",
            "GSAT0219": "B04",
            "GSAT0220": "B01",
            "GSAT0221": "B02",
            "GSAT0222": "B07",
            "GSAT0223": "B03",
            "GSAT0224": "B15",
            "GSAT0225": "C05",
            "GSAT0226": "A02",
            "GSAT0227": "C12",
            "GSAT0232": "A17",
    }

    station_def = {
            "altc": "Alta (Norway)",
            "brux": "Brussels (Belgium)",
            "cpvg": "Cap-Vert (Cabo Verde)",
            "koug": "Kourou (French Guiana)",
            "hofs": "Hoefn (Iceland)",
            "hons": "Honningsvag (Norway)",
            "janm": "Jan Mayen (Norway)",
            "krss": "Kristiansand (Norway)",
            "mas1": "Maspalomas (Spain)",
            "nabd": "Ny Alesund (Norway)",
            "nklg": "N' Koltang (Gabon)",
            "vegs": "Vega (Norway)",
    }

    target_def = {
            "hpe": "",
            "site_vel_3d": "",
            "sisre": "2",
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
        used_date = datetime.strptime(f"{date}-01", "%Y-%b-%d")
        sat_info = atx.get_satellite_info(satellite, used_date)

        batch = batch_def[sat_info["sat_type"]] if sat_info["sat_type"] in batch_def.keys() else ""

        if sat_info["sat_code"].startswith("E"):
            svn = sat_info["sat_code"].replace("E", "GSAT0")
        elif sat_info["sat_code"].startswith("G"):
            svn = sat_info["sat_code"].replace("G0", "SVN")
        else:
            svn = sat_info["sat_code"]
        slot = slot_def[svn] if svn in slot_def.keys() else ""

    return [
        constellation, # Constellation
        "Open Service", # Service Line
        category_def[kpi], # Service Category
        business_def[constellation][kpi], # Business Service
        batch, # Batch
        svn, # Satellite
        satellite, # PRN
        slot, # Slot
        station_to_write, # GSS Site
        station.upper(), # Station
        "OS", # Service
        type_, # Type
        mode_to_write, # Mode
        target_def[kpi], # Target
        unit_def[kpi], # Unit
        date, # Month/Year
        f"{result:.3f}", # Result
    ]


