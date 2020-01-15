"""A dataset for handling time series data

Description:
------------

Mainly based on the Midgard dataset. Some extra fields are added, typically
those that rely on features that require compilation of external libraries,
like earth orientation parameters.

"""

# Add Where fieldtypes before importing dataset
from where.data import fieldtypes  # noqa

# Standard library imports
from datetime import date
import getpass
import os
from typing import Dict, Optional, Union

# Midgard imports
from midgard.data import dataset as mg_dataset
from midgard.data.dataset import Dataset as MgDataset

# Where imports
import where
from where.lib import config


class Dataset(MgDataset):
    """A Midgard dataset customized for use inside of Where"""

    version = f"{MgDataset.version}, Where v{where.__version__}"

    def __init__(
        self,
        num_obs: int = 0,
        rundate: Optional[date] = None,
        pipeline: str = "",
        stage: str = "",
        label: Union[str, int] = 1,
        user: str = "",
        id: str = "",
        **dset_args: str,
    ):
        """Create dataset, and optionally set dataset variables"""
        super().__init__(num_obs=num_obs)
        self.vars.update(self._dset_vars(rundate=rundate, pipeline=pipeline, stage=stage, label=label, **dset_args))
        self.analysis = self._analysis_vars(rundate, user, id)  # Dataset variables that will not be stored on file.

    @classmethod
    def read(
        cls,
        rundate: date,
        pipeline: str,
        stage: str,
        label: Union[str, int],
        user: str = "",
        id: str = "",
        **dset_args: str,
    ) -> "Dataset":
        """Read a dataset from file"""

        dset_vars = cls._dset_vars(rundate=rundate, pipeline=pipeline, stage=stage, label=label, **dset_args)
        if label == "last":
            ids = config.files.glob_variable("dataset", variable="label", pattern=r"[\w]+", file_vars=dset_vars)
            try:
                dset_vars.update(label=max(ids))
            except ValueError:
                # No ids found
                raise ValueError(f"No labels found for date {rundate}, pipeline {pipeline}, stage {stage}")
        analysis_vars = cls._analysis_vars(rundate, user, id)

        file_vars = dset_vars.copy()
        for k, v in analysis_vars.items():
            if k not in file_vars:
                file_vars[k] = v
        file_path = config.files.path("dataset", file_vars=file_vars)

        dset = super().read(file_path)
        dset.vars.update(dset_vars)
        dset.analysis = analysis_vars  # Dataset variables that will not be stored on file.
        return dset

    def delete_stage(self, stage, **kwargs):
        file_vars = dict(stage=stage)
        file_vars.update(kwargs)
        file_paths = config.files.glob_paths("dataset", file_vars=file_vars)
        for file_path in file_paths:
            os.remove(file_path)

    def write(self, write_level: str = None) -> None:
        """Write a dataset to file"""
        file_vars = self.vars.copy()
        for k, v in self.analysis.items():
            if k not in file_vars:
                file_vars[k] = v
        file_path = config.files.path("dataset", file_vars=file_vars)
        write_level = config.tech.write_level.str if write_level is None else write_level
        super().write(file_path, write_level=write_level)

    def write_as(
        self,
        rundate: date = None,
        pipeline: str = None,
        stage: str = None,
        label: Union[str, int] = None,
        write_level: str = None,
        **kwargs: str,
    ) -> None:
        """Write a dataset to file, renaming the given dataset variables"""

        # Construct all variables
        all_vars = self._dset_vars(rundate=rundate, pipeline=pipeline, stage=stage, label=str(label), **kwargs)
        # Pick out the variables that should be updated
        update_vars = {k: v for k, v in all_vars.items() if v}
        self.vars.update(update_vars)

        # Write the dataset
        self.write(write_level=write_level)

    @staticmethod
    def _dset_vars(
        rundate: date, pipeline: str, stage: str, label: Union[str, int], **dset_args: str
    ) -> Dict[str, str]:
        """Create a dataset variable dictionary"""
        rundate_str = rundate.strftime(config.FMT_date) if rundate else ""
        dset_vars = dict(rundate=rundate_str, pipeline=pipeline, stage=stage, label=label, **dset_args)
        return dset_vars

    @staticmethod
    def _analysis_vars(rundate, user="", id=""):
        analysis_vars = dict()
        analysis_vars["rundate"] = rundate  # Rundate as date
        analysis_vars["user"] = user if user else getpass.getuser().lower()
        analysis_vars["id"] = id
        analysis_vars.update(config.date_vars(rundate))
        return analysis_vars


for field_type in fieldtypes.names():
    setattr(Dataset, f"add_{field_type}", mg_dataset._add_field_factory(field_type))
