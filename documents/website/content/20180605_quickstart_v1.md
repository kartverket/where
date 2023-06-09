Title: Quickstart
Date: 2018-6-5
Tags: vlbi, where, there
Status: hidden
Summary: Get started using Where for VLBI analysis.

**Where** is a geodetic software package, used for professional analysis of VLBI data.

This is the quickstart guide for Where version 1 (or older). For the latest version see [quickstart guide]({filename}20230103_quickstart_v2.md).

## Running a VLBI analysis

To do an analysis with Where you only need to specify which analysis
[**pipeline**]({filename}20180604_about.md#pipeline) you want to run
(i.e. VLBI), as well as the **date** and **session** of your data. The most
basic Where command looks something like

    where 2016 12 27 --vlbi --session=XA

This command can be extended by setting any of the many configuration options
available for an analysis. More details about the configuration settings are in
a [separate page]({filename}20180605_where.md#configuration) You can also run
`where --help` for a quick overview.



## Exploring the analysis results

To look at the results of your analysis, you can use the companion tool
**There**.  The `there` command works similarly to `where`, so to see the
results from the analysis earlier you can do:

    there 2016 12 27 --vlbi --session=XA

More details about There can be found on a
[separate page]({filename}20180605_there.md).


### Exploring with Python

If you are comfortable using the Python programming language, you can also
explore the analysis results in a Python interactive shell. Run Where with the
`-I` / `--interactive` option. This will open an interactive Python (`ipython`)
shell the analysis preloaded into variables named `v0`, `v1` and so on:

    where 2016 12 27 --vlbi --session=XA -I
