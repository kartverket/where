import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import functools

from midgard import parsers


def make_plots(sta_data):
    # trop_tot plot
    names = {"trop_tot": "Total troposphere delay [mm]",
            "trop_gn_tot": "Gradient north [mm]",
            "trop_ge_tot": "Gradient east [mm]"}
    
    def plot(field):
        fig = plt.figure(figsize=(12, 8), dpi=150)
        for sta in sta_data.keys():
            x = sta_data[sta]["epoch"]
            y = sta_data[sta][field]
            plt.scatter(x, y, marker=".", label=sta)
            plt.ylabel(names[field])
        plt.legend()
        plt.savefig(f"img/GNSS_{field}.png", bbox_inches='tight')
    
    for n in names.keys():
        plot(n)


def do_statistics(sta_data):
    field_dict = {}
    for sta in sta_data.keys():
        for field in sta_data[sta].keys():
            
            field_dict.setdefault(field, {}).update({sta:sta_data[sta][field]})

    
    for field, value in field_dict.items():
        print(f"{field}")
        df = pd.DataFrame(value)
        print(df.corr())
    

def main():
    path = "/home/kirann/Downloads/bern_trop"
    stations = ["NABG", "NABD", "NYA1", "NYAL"]
    sta_data = {s:{} for s in stations}
    
    # 2020 was a leap year
    for doy in range(1, 367):
        filename = f"{path}/F1_20{doy:03}0.TRO"
        
        p = parsers.parse_file("sinex_tro", filename)
        if all(sta in p.data for sta in stations):
            # To avoid arrays of different size only use data when all stations are in the file
            for sta in stations:
                try:
                    for k, v in p.data[sta].items():
                        sta_data[sta].setdefault(k, []).append(v)
                except KeyError:
                    print(f"{sta} missing in {p.file_path}")


    make_plots(sta_data)
    do_statistics(sta_data)

if __name__ == "__main__":
    main()