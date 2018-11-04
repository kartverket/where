#!/usr/bin/env python3
"""Set up local configuration files for Where

Usage:
------

    python3 config_wizard.py <local_config_file> [options]

The options can be any of

    -h, --help                - Show this help message
        --ignore-existing     - Do not run the wizard if the config file exists
        --overwrite-existing  - Overwrite the config file if it exists


Description:
------------

This script starts a wizard for creating a local configuration file.
"""
import sys
import pathlib

from midgard.dev.console import color
from midgard.dev.exceptions import MissingConfigurationError
from where.lib import config


def main() -> None:
    """Check command line parameters
    """
    # Read from command line
    opts = [o for o in sys.argv[1:] if o.startswith("-")]
    args = [a for a in sys.argv[1:] if not a.startswith("-")]

    # Show help message
    if "-h" in opts or "--help" in opts or not args:
        help()

    # Check that path for local configuration is valid and does not already exist
    cfg_path = pathlib.Path(args[0])
    if not str(cfg_path).endswith("_local.conf"):
        help(f"Name of configuration file should end with '_local.conf'")
    if cfg_path.exists():
        if "--ignore-existing" in opts:
            raise SystemExit
        elif "--overwrite-existing" not in opts:
            help(f"The file {cfg_path} already exists")

    # Read main configuration
    cfg_name = cfg_path.name[:-11]
    try:
        cfg = read_configuration(cfg_name)
    except MissingConfigurationError as e:
        help(str(e))

    # Set up local config
    print(f"Setting up configuration file {cfg_path}")
    cfg_local = config_wizard(cfg)

    # Save local config to file
    help_path = cfg_path.with_suffix(".template")
    help_text = help_path.read_text() if help_path.exists() else ""
    with open(cfg_path, mode="w") as fid:
        fid.write(help_text)
        fid.write(cfg_local.as_str(width=80, key_width=20) + "\n")


def read_configuration(cfg_name: str) -> config.Configuration:
    """Read a configuration from disk

    """
    cfg = config.Configuration("main")
    file_paths = list(config.config_paths(cfg_name))
    if not file_paths:
        raise MissingConfigurationError(f"Where has no configuration named {cfg_name!r}")
    for file_path in file_paths:
        if not str(file_path).endswith("_local.conf"):
            cfg.update_from_file(file_path)

    return cfg


def config_wizard(cfg: config.Configuration) -> config.Configuration:
    """Set up a local configuration file

    Args:
        cfg:  Main configuration to create a local configuration for.

    Returns:
        Local configuration
    """
    # Local configuration
    cfg_local = config.Configuration("local")

    # Loop through sections and entries of main configuration
    for section in cfg.section_names:
        for key, entry in cfg[section].items():
            if "wizard" not in entry.meta:
                continue

            description = entry.meta["wizard"] if entry.meta["wizard"] else entry.meta.get("help", "")
            print(f"\n+ Enter value for {color.Style.BRIGHT}{section}:{key}{color.Style.NORMAL} ({description})")
            value = input(f"  [{entry.str}] ")
            if value:
                cfg_local.update(section=section, key=key, value=value, source=__name__)
    return cfg_local


def help(error_msg: str = "") -> None:
    """Print the help message from the module doc-string.

    Args:
        error_msg:  Optional extra error message that will be printed below the help message.
    """
    print(__doc__)
    if error_msg:
        print(f"\n{color.Fore.RED}Error: {error_msg}")

    raise SystemExit


if __name__ == "__main__":
    main()
