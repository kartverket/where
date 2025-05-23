#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals

AUTHOR = "Norwegian Mapping Authority"
SITENAME = "Where"
SITEURL = "https://kartverket.github.io/where"


THEME = 'pelican-bootstrap3'
JINJA_ENVIRONMENT = {'extensions': ['jinja2.ext.i18n']}
PLUGIN_PATHS = ['pelican-plugins']
#PLUGINS = ['i18n_subsites', "pin_to_top"]
# pin_to_top does not seem to be maintained any more and has not migrated to new plugins structure
PLUGINS = ["i18n_subsites"] 
BOOTSTRAP_THEME = 'cerulean'  # https://bootswatch.com/  cosmo, sandstone, lumen

MARKDOWN = {
    "extension_configs": {
        "markdown.extensions.toc": {"title": "Table of Contents"},
    }
}

PATH = "content"
STATIC_PATHS = ["extras", "images"]
EXTRA_PATH_METADATA = {"extras/custom.css": {"path": "static/css/custom.css"}}
CUSTOM_CSS = "static/css/custom.css"
FAVICON = "images/favicon.ico"
TIMEZONE = "Europe/Paris"
DEFAULT_LANG = "en"

BANNER = "images/banner.jpg"
BANNER_SUBTITLE = "Version 2.1.3"
# BANNER_ALL_PAGES = True

DISPLAY_PAGES_ON_MENU = True
DISPLAY_CATEGORIES_ON_MENU = True
# BOOTSTRAP_NAVBAR_INVERSE = True

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Blogroll
GITHUB_URL = "https://github.com/kartverket/where"
LINKS = (
    ("Midgard", "https://github.com/kartverket/midgard"),
    # ("Norwegian Mapping Authority", "https://kartverket.no/en/"),
    ("Anaconda Python", "https://www.anaconda.com/download"),
    ("International VLBI Service", "https://ivscc.gsfc.nasa.gov/"),
    ("IERS Conventions 2010", "https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html"),
)

# Social widget
SOCIAL = (
    ("github", "https://github.com/kartverket/where"),
    # ("twitter", "#"),
)


DEFAULT_PAGINATION = 10
USE_PAGER = False
# Uncomment following line if you want document-relative URLs when developing
#RELATIVE_URLS = True
