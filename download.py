#!/usr/bin/env python3
"""Download external code or files needed by Where

Usage:
------

    python3 download.py library


Description:
------------

This script downloads external libraries used by Where.
"""
from collections import OrderedDict
from configparser import ConfigParser
import pathlib
import subprocess
import sys

import pycurl


class CasedConfigParser(ConfigParser):
    """ConfigParser with case-sensitive keys"""

    def optionxform(self, optionstr):
        return optionstr

    def read(self, path):
        """Replace f-style variables with values from the __vars__ section after the config file is read"""
        super().read(path)

        fmt_vars = {k: v for k, v in self["__vars__"].items()}
        for name, section in self._sections.items():
            if name.startswith("__"):
                continue

            new_section = OrderedDict()
            for key, value in section.items():
                new_section[key.format(**fmt_vars)] = value.format(**fmt_vars)
            self._sections[name] = new_section


def main() -> None:
    """Check command line parameters
    """
    # Read from command line
    opts = [o for o in sys.argv[1:] if o.startswith("-")]
    args = [a for a in sys.argv[1:] if not a.startswith("-")]

    # Show help message
    if "-h" in opts or "--help" in opts or not args:
        help()

    # Check if config for library exists
    library = args[0]
    cfg_path = pathlib.Path.cwd() / "config" / f"download_{library}.conf"
    if not cfg_path.exists():
        help()

    # Download library
    download(cfg_path)


def download(cfg_path: pathlib.Path) -> None:
    """Download files

    Information about the files that should be downloaded are given in a
    configuration file.

    Args:
        cfg_path:  Path to configuration file
    """
    # Read configuration
    cfg = CasedConfigParser()
    cfg.read(cfg_path)

    # Report on progress
    lib_name = cfg["library"].get("name")
    lib_url = cfg["library"].get("source_url")
    target_dir = pathlib.Path(cfg["library"].get("target_dir"))
    print(f"Downloading {lib_name} from '{lib_url}'")

    # Do preprocessing
    process(cfg, "preprocess")

    # Download files
    file_sections = [s for s in cfg.sections() if s.startswith("files")]
    for section in file_sections:
        target_directory = cfg[section].get("__target__", "")
        for file_name, lib_directory in cfg[section].items():
            if file_name.startswith("__"):
                continue
            if not lib_directory.endswith("/"):
                lib_directory += "/"

            target_path = target_dir / target_directory / file_name
            print(f"+ {file_name:<20} -> {target_path}")
            download_file(f"{lib_url}/{lib_directory}{file_name}", target_path)

    # Do postprocessing
    process(cfg, "postprocess")


def download_file(url: str, path: pathlib.Path) -> None:
    """Download a single file from a url, save at given path

    Args:
        url:   URL where file will be downloaded from.
        path:  File path where file will be saved.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, mode="wb") as fid:
        c = pycurl.Curl()
        c.setopt(c.URL, url)
        c.setopt(c.WRITEDATA, fid)
        c.perform()
        response_code = c.getinfo(c.RESPONSE_CODE)
        if not (200 <= response_code <= 299):
            print(f"There was a problem downloading {c.getinfo(c.EFFECTIVE_URL)} ({response_code})")
        c.close()


def process(cfg: ConfigParser, section: str) -> None:
    """Process commands in a section

    Args:
        cfg:      Configuration.
        section:  Name of section in configuration containing commands to be processed.
    """
    if section not in cfg:
        return

    for command, hint in cfg[section].items():
        if command.startswith("__"):
            continue

        print(f"# {command}")
        try:
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print(f"Command '{command}' failed.")
            input(f"Please {hint} manually, and hit enter afterwards")


def help() -> None:

    """Print the help message from the module doc-string.
    """
    print(__doc__)
    raise SystemExit


if __name__ == "__main__":
    main()
