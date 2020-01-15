# """Framework for calculating models
#
# Description:
# ------------
#
# Each model should be defined in a separate .py-file inside one of the model directories, delay, site, orbit and
# satellite. The function inside the .py-file that should be called need to be decorated with the
# :func:`~midgard.dev.plugins.register` decorator as follows::
#
#     from midgard.dev import plugins
#
#     @plugins.register
#     def vlbi_vacuum_delay(dset):
#         ...
#
# See the corresponding :file:`__init__.py`-files for more information.
#
#
# See also:
# ---------
#
# * :mod:`where.models.delay`
# * :mod:`where.models.site`
# * :mod:`where.models.orbit`
#
# """
# # External library imports
#
# # Where imports
# from where.lib import log
#
#
# def _calculate_model(calculate_func, config_key, dset_in, dset_out, shape, write_levels=None):
#     """Call models and store output in dataset
#
#     If the model output is empty, we still create a dummy field in the table only containing zeros. This is done to
#     assert that the table will always exist after doing a `models.calculate...`-call.
#
#     Args:
#         calculate_func (Function):  The function that calls models.
#         config_key (String):        Key in config with list of models, also table the model output is stored in.
#         dset_in (Dataset):          Dataset to read data from.
#         dset_out (Dataset):         Dataset to store data to.
#         shape (Tuple of int):       Shape of output.
#     """
#     dset_out = dset_in if dset_out is None else dset_out
#     write_levels = dict() if write_levels is None else write_levels
#
#     model_output = calculate_func(config_key, dset_in)
#     if not model_output:
#         # TODO update to dataset3
#         model_name = "{}_zeros".format(config_key)
#         if model_name not in dset_out.fields:
#             dset_out.add_float(
#                 model_name,
#                 table=config_key,
#                 shape=shape,
#                 unit="meter",
#                 write_level=write_levels.get(model_name, "analysis"),
#             )
#
#
#     for model_name, values in sorted(model_output.items()):
#         for multiplier, pos_delta in zip(dset_out.for_each_suffix("station"), values):
#             field_name = f"{config_key}.{model_name}{dset_out.default_field_suffix}"
#             if field_name in dset_out.fields:
#                 dset_out[field_name][:] = values
#             else:
#                 dset_out.add_position_delta(field_name, pos_delta, multiplier=multiplier, write_level=write_levels.get(model_name, "analysis"))
#         log.info(f"Average correction = {dset_out.rms(f'{config_key}.{model_name}'):14.5f} in {model_name} model")
#
#
#         log.info(f"Average correction = {dset_out.rms(model_name):14.5f} in {model_name} model")
