import importlib
import numpy as np
import os
import sys
import time


def calculate_tax(version, pop_income):
    module_name = "tax_{version}".format(version=version)
    taxes = importlib.import_module(module_name)

    num_pop = len(pop_income)
    mean_income = np.mean(pop_income)
    print("Antall personer:         {:>10d}".format(num_pop))
    print("Gj.snittlig inntekt:     {:>13.2f} NOK".format(mean_income))

    start_time = time.perf_counter()
    pop_tax = taxes.calculate_pop_tax(pop_income)
    running_time = 1e3 * (time.perf_counter() - start_time)
    mean_tax = np.mean(pop_tax)
    print("Gj.snittlig skatt:       {:>13.2f} NOK".format(mean_tax))

    percent_tax = 100 * mean_tax / mean_income
    print("Gj.snittlig skatteprosent:      {:>6.2f} %".format(percent_tax))
    print("Kjoretid:                   {:>10.2f} ms".format(running_time))

    return running_time


def main():
    num_pop = 3100000
    mu = np.log(450000)
    sigma = 1
    np.random.seed(1)

    pop_income = np.random.lognormal(mu, sigma, num_pop)

    versions = sys.argv[1:]
    if not versions or "all" in versions:
        versions = reversed(
            sorted([f[4:].split(".")[0] for f in os.listdir() if f.startswith("tax_") and ".py" in f[-4:]])
        )

    running_times = dict()
    for version in versions:
        print("### {version}".format(version=version))
        running_times[version] = calculate_tax(version, pop_income)

    print("\nVersion             Running Time   Speed-up")
    max_time = running_times[max(running_times, key=running_times.get)]
    for version in reversed(sorted(running_times, key=running_times.get)):
        run_time = running_times[version]
        rel_time = max_time / run_time
        print(
            "{version:<20s} {run_time:8.2f} ms   {rel_time:8.2f}".format(
                version=version, run_time=run_time, rel_time=rel_time
            )
        )


if __name__ == "__main__":
    main()
