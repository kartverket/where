"""Where library module for opening files

Example:
--------

from where.lib import files
with files.open('eopc04_iau', mode='rt') as fid:
    for line in fid:
        print(line.strip())

Description:
------------

This module handles opening of files. All regular Where files should be registered in the Where file list, and then
opened using files.open. Other files could be opened using files.open_path, although this should be the exception,
rather than the rule.

The files.open and files.open_path functions are both wrappers around the built-in open function, and behave mainly
similar. In particular, they accept all the same keyword arguments (like for instance mode). Furthermore, to make sure
files are properly closed they should normally be used with a context manager as in the example above.

"""

# Standard library imports
import builtins
from contextlib import contextmanager
import gzip
import itertools
import pathlib
import re
import shutil

# External library imports
import h5py
import pycurl

# Midgard imports
from midgard.dev import console

# Where imports
from where.lib import config
from where.lib import log


class URL(str):
    """Simple wrapper around String to have URLs work similar to pathlib.Path

    TODO: Remove when files.path() is removed
    """

    @property
    def name(self):
        """Part of URL after the last /"""
        return self.split("/")[-1]

    def with_name(self, name):
        """Replace part of URL after the last / with a new name

        Args:
            name (String):  New name.

        Return:
            URL:  URL with part after the last / replaced with the new name.
        """
        base_url, slash, _ = self.rpartition("/")
        return self.__class__(f"{base_url}{slash}{name}")

    def exists(self):
        """Check whether the given URL returns a valid document

        Try to download the first byte of the document (avoid downloading a big file if it exists).

        Return:
            Boolean:  True if URL leads to a valid document, False otherwise.
        """
        c = pycurl.Curl()
        c.setopt(c.URL, self)
        c.setopt(c.RANGE, "0-0")
        try:
            c.perform()
        except pycurl.error:
            return False
        return True


@contextmanager
def open(file_key, file_vars=None, create_dirs=False, is_zipped=None, download_missing=True, **kwargs):
    """Open a Where file

    Open a Where file based on file key which is looked up in the Where file list.

    The function automatically handles reading from gzipped files if the filename is specified with the special
    {gz}-ending (including the curly braces) in the file list. In that case, the mode should be specified to be 'rt' if
    the contents of the file should be treated as text. If both a zipped and an unzipped version is available, the
    zipped version is used. This can be overridden by specifying True or False for the is_zipped-parameter.

    This function behaves similar to the built-in open-function, and should typically be used with a context manager as
    follows:

    Example:
        with files.open('eopc04_iau', mode='rt') as fid:
            for line in fid:
                print(line.strip())

    Args:
        file_key:    String that is looked up in the Where file list.
        file_vars:   Dict, used to replace variables in file name and path.
        create_dirs: True or False, if True missing directories are created.
        kwargs:      All keyword arguments are passed on to open_path.

    Returns:
        File object representing the file.
    """
    import sys

    caller = sys._getframe(2)
    func_name = caller.f_code.co_name
    file_name = caller.f_code.co_filename
    line_num = caller.f_lineno
    log.dev(
        f"{file_name} ({line_num}) {func_name}: 'lib.files.open()' is deprecated. Use 'lib.config.files.open()' instead"
    )

    download_missing = download_missing and "r" in kwargs.get("mode", "r")
    file_path = path(file_key, file_vars, is_zipped=is_zipped, download_missing=download_missing)
    kwargs.setdefault("encoding", encoding(file_key))
    try:
        with open_path(
            file_path, description=file_key, create_dirs=create_dirs, is_zipped=is_path_zipped(file_path), **kwargs
        ) as fid:
            yield fid
    except Exception:
        raise


@contextmanager
def open_datafile(file_key, file_vars=None, mode="r", create_dirs=True, write_log=True, **kwargs):
    """Open an HDF5 datafile

    TODO: Remove when switching to dataset 3
    """
    file_path = path(file_key, file_vars)
    if write_log:
        _log_file_open(file_path, file_key, mode)
    if create_dirs:
        file_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with h5py.File(file_path, mode=mode, **kwargs) as datafile:
            yield datafile
    except Exception:
        raise


@contextmanager
def open_path(file_path, description="", mode="rt", create_dirs=False, is_zipped=None, write_log=True, **kwargs):
    """Open a local file

    TODO: Possibly not needed, if needed -> move to a method on midgard.config.files.FileConfiguration

    Open a local file based on file name. This function behaves similar to the built-in open-function, the difference
    is that we do some extra logging. The function should typically be used with a context manager as follows:

    Example:
        with files.open_path('local_output.txt', mode='rt') as fid:
            for line in fid:
                print(line.strip())

    Args:
        file_path (String/Path):   The file_path, should be a full path.
        description (String):      Description used for logging.
        mode (String):             Same as for the built-in open, usually 'rt' or 'wt'.
        create_dirs (Boolean):     True or False, if True missing directories are created.
        is_zipped (Boolean):       True or False, if True the gzip module will be used.
        kwargs:                    All keyword arguments are passed on to the built-in open.

    Returns:
        File object representing the file.
    """
    if write_log:
        _log_file_open(file_path, description, mode)
    if create_dirs:
        file_path.parent.mkdir(parents=True, exist_ok=True)
    is_zipped = is_path_zipped(file_path) if is_zipped is None else is_zipped

    try:
        if is_zipped:
            with gzip.open(file_path, mode=mode, **kwargs) as fid:
                yield fid
        else:
            with builtins.open(file_path, mode=mode, **kwargs) as fid:
                yield fid
    except Exception:
        raise


def encoding(file_key):
    """Look up the encoding for a given file key

    TODO: Remove when files.path() is removed

    Args:
        file_key (String):  Key that is looked up in the Where file list.

    Returns:
        String:  Name of encoding. If encoding is not specified None is returned.
    """
    file_encoding = config.files.get("encoding", section=file_key, default="").str
    return file_encoding or None


def path(file_key, file_vars=None, default=None, is_zipped=None, download_missing=False, use_aliases=True):
    """Construct a filepath for a given file with variables

    If `is_zipped` is None, and the file_path contains `<filename>{gz}`, the file will be assumed to be a gzip-file if
    there exists a file named `<filename>.gz`.

    When setting `use_aliases` to True, the aliases as specified in the files configuration file represent alternative
    filenames. In particular,

        + if directory / file_name exists it is returned
        + otherwise the first directory / alias that exists is returned
        + if none of these exist, directory / file_name is returned

    Args:
        file_key (String):        Key that is looked up in the Where file list.
        file_vars (Dict):         Values used to replace variables in file name and path.
        default (String):         Value to use for variables that are not in file_vars.
        is_zipped (Bool/None):    True, False or None. If True, open with gzip. If None automatically decide.
        download_missing (Bool):  Whether to try to download missing files.
        use_aliases (Bool):       Fall back on aliases if file does not exist.

    Return:
        Path: Full path with replaced variables in file name and path.
    """
    import sys

    caller = sys._getframe(1)
    func_name = caller.f_code.co_name
    file_name = caller.f_code.co_filename
    line_num = caller.f_lineno
    log.dev(
        f"{file_name} ({line_num}) {func_name}: 'lib.files.path()' is deprecated. Use 'lib.config.files.path()' instead"
    )

    file_vars = dict() if file_vars is None else file_vars
    directory = config.files[file_key].directory.replace(default=default, **file_vars).path
    file_name = config.files[file_key].filename.replace(default=default, **file_vars).path
    file_path = _replace_gz(directory / file_name)

    # Check for aliases
    if use_aliases and not path_exists(file_path):
        aliases = config.files.get("aliases", section=file_key, default="").replace(default=default, **file_vars).list
        for alias in aliases:
            aliased_path = _replace_gz(file_path.with_name(alias))
            if path_exists(aliased_path):
                return aliased_path

    # Try to download the file if it is missing
    if download_missing and not path_exists(file_path):
        downloaded_path = download_file(file_key, file_vars)
        if downloaded_path is not None:
            file_path = downloaded_path

    return file_path


def url(file_key, file_vars=None, default=None, is_zipped=None, use_aliases=True):
    """Construct a URL for a given file with variables

    TODO: Remove when files.path() is removed

    If `is_zipped` is None, and the url contains `<filename>{gz}`, the url will be assumed to point to a gzip-file if
    there exists a file named `<filename>.gz` on the server.

    Args:
        file_key (String):        Key that is looked up in the Where file list.
        file_vars (Dict):         Values used to replace variables in file name and path.
        default (String):         Value to use for variables that are not in file_vars.
        is_zipped (Bool/None):    True, False or None. If True, open with gzip. If None automatically decide.

    Return:
        String: Full URL with replaced variables in file name and url.
    """
    file_vars = dict() if file_vars is None else file_vars
    base_url = config.files[file_key].url.replace(default=default, **file_vars).str.rstrip("/")
    file_name = config.files[file_key].filename.replace(default=default, **file_vars).str
    file_url = _replace_gz(URL(f"{base_url}/{file_name}"), is_zipped=is_zipped)

    # Check for aliases
    if use_aliases and not file_url.exists():
        aliases = config.files.get("aliases", section=file_key, default="").replace(default=default, **file_vars).list
        for alias in aliases:
            aliased_url = _replace_gz(file_url.with_name(alias), is_zipped=is_zipped)
            if aliased_url.exists():
                return aliased_url

    return file_url


def _replace_gz(file_path, is_zipped=None):
    """Replace the {gz} pattern with '.gz' or '' depending on whether the file is zipped

    TODO: Remove when files.path() is removed

    If `is_zipped` is None, and the file_path contains `<filename>{gz}`, the file will be assumed to be a gzip-file if
    there exists a file named `<filename>.gz`.

    Args:
        file_path (Path):       Path to a file
        is_zipped (Bool/None):  True, False or None. If True, open with gzip. If None automatically decide.

    Returns:
        Path:  File path with {gz} replaced.
    """
    if "{gz}" not in file_path.name:
        return file_path

    if is_zipped is None:
        is_zipped = file_path.with_name(file_path.name.replace("{gz}", ".gz")).exists()
    if is_zipped:
        return file_path.with_name(file_path.name.replace("{gz}", ".gz"))
    else:
        return file_path.with_name(file_path.name.replace("{gz}", ""))


def empty_file(file_path):
    """Check if a file is empty

    TODO: Remove when files.path() is removed

    Args:
        file_path (Path):  Path to a file.

    Returns:
        Bool:  Whether path is empty or not.
    """
    if not path_exists(file_path):
        log.error(f"File '{file_path}' does not exist.")

    return False if file_path.stat().st_size > 0 else True


def path_exists(file_path):
    """Check if a path exists

    Unfortunately, Windows throws an error when doing file_path.exists() if the file path contains wildcard characters
    like *. Thus, we wrap this in a check on whether the file path syntax is correct before calling file_path.exists.
    If the file path contains non-path characters, the file path can not exist.

    TODO: Remove when files.path() is removed

    Args:
        file_path (Path):  Path to a file.

    Returns:
        Bool:  Whether path exists or not.
    """
    try:
        return file_path.exists()
    except OSError:
        return False


def is_path_zipped(file_path):
    """Indicate whether a path is to a gzipped file or not

    For now, this simply checks whether the path ends in .gz or not.

    TODO: Remove when files.path() is removed

    Args:
        file_path (Path):  Path to a file.

    Returns:
        Boolean:   True if path is to a gzipped file, False otherwise.
    """
    try:
        file_name = file_path.name  # Assume file_path is Path-object
    except AttributeError:
        file_name = file_path  # Fall back to file_path being string
    return file_name.endswith(".gz")


def print_file(file_path):
    """Print the contents of a file to the console

    TODO: Only used as fallback if editor package is not installed. Can probably be removed?

    Args:
        file_path (Path):  Path to a file.
    """
    log.info(f"Printing contents of {file_path}")
    with open_path(file_path, mode="r") as fid:
        for line in fid:
            print(line.rstrip())


def delete_file(file_path):
    """Deletes a file

    Does not delete directories

    TODO: Only used by _kalman. Can maybe use pathlib.Path.unlink directly?

    Args:
        file_path(Path):     Path to a file
    """
    try:
        file_path.unlink()
        log.debug(f"Deleted {file_path}")
    except OSError:
        log.warn(f"Unable to delete {file_path}")


def publish_files(publish=None):
    """Publish files to specified directories

    TODO: Not moved! Seems Where-specific?

    The publish string should list file_keys specified in files.conf. Each file_key needs to have a field named publish
    specifying a directory the file should be copied to.

    Args:
        publish (String):   List of file_keys that will be published.
    """
    if not config.where.files.publish.bool:
        return

    publish_list = config.tech.get("files_to_publish", value=publish).list
    for file_key in publish_list:
        try:
            source = path(file_key)
        except KeyError:
            log.error(f"File key '{file_key}' in publish configuration is unknown. Ignored")
            continue
        if not source.exists():
            try:
                log.error(f"File '{source}' (file key='{file_key}') does not exist, and can not be published")
            except KeyError:
                log.error(f"File key='{file_key}' has incomplete filename information and can not be published")
            continue

        try:
            destinations = config.files[file_key].publish.replaced.as_list(convert=pathlib.Path)
        except AttributeError:
            log.error(f"File key '{file_key}' does not specify 'publish' directory in file configuration. Ignored")
            continue

        # Copy file to destinations
        for destination in destinations:
            log.info(f"Publishing {file_key}-file {source} to {destination}")
            destination.mkdir(parents=True, exist_ok=True)
            shutil.copy(source, destination)


def download_file(file_key, file_vars=None, create_dirs=True):
    """Download a file from the web and save it to disk

    TODO: Remove when files.path() is removed

    Use pycurl (libcurl) to do the actual downloading. Request might be nicer for this, but turned out to be much
    slower (and in practice unusable for bigger files) and also not really supporting ftp-downloads.

    Args:
        file_key (String):   File key that should be downloaded.
        file_vars (Dict):    File variables used to find path from file_key.
        create_dirs (Bool):  Create directories as necessary before downloading file.
    """
    if (
        not config.where.files.download_missing.bool
        or "url" not in config.files[file_key]
        or not config.files[file_key].url.str
    ):
        return None

    file_path = path(file_key, file_vars=file_vars, download_missing=False)
    if file_path.exists():
        return None
    if create_dirs:
        file_path.parent.mkdir(parents=True, exist_ok=True)

    file_url = url(file_key, file_vars=file_vars)
    file_path = file_path.with_name(file_url.name)
    log.info(f"Download {file_key} from '{file_url}' to '{file_path}'")
    with builtins.open(file_path, mode="wb") as fid:
        c = pycurl.Curl()
        c.setopt(c.URL, file_url)
        c.setopt(c.WRITEDATA, fid)
        try:
            c.perform()
            if not (200 <= c.getinfo(c.HTTP_CODE) <= 299):
                raise pycurl.error()
        except pycurl.error:
            log.error(f"Problem downloading file: {c.getinfo(c.EFFECTIVE_URL)} ({c.getinfo(c.HTTP_CODE)})")
            if file_path.exists():  # Print first 10 lines to console
                head_of_file = f"Contents of '{file_path}':\n" + "\n".join(file_path.read_text().split("\n")[:10])
                print(console.indent(head_of_file, num_spaces=8))
                file_path.unlink()
            log.warn(f"Try to download '{file_url}' manually and save it at '{file_path}'")
        else:
            log.info(f"Done downloading {file_key}")
        finally:
            c.close()
    return file_path


def glob_paths(file_key, file_vars=None, is_zipped=None):
    """Find all filepaths matching a filename pattern

    Using pathlib.Path.glob() here is not trivial because we need to split into a base directory to start searching
    from and a pattern which may include directories. With glob.glob() this is trivial. The downside is that it only
    returns strings and not pathlib.Paths.
    """
    import sys

    caller = sys._getframe(1)
    func_name = caller.f_code.co_name
    file_name = caller.f_code.co_filename
    line_num = caller.f_lineno
    log.dev(
        f"{file_name} ({line_num}) {func_name}: 'lib.files.glob_paths()' is deprecated. Use 'lib.config.files.glob_paths()' instead"
    )

    path_string = str(path(file_key, file_vars, default="*", is_zipped=is_zipped))
    glob_path = pathlib.Path(re.sub(r"\*+", "*", path_string))
    idx = min((i for i, p in enumerate(glob_path.parts) if "*" in p), default=len(glob_path.parts) - 1)
    glob_base = pathlib.Path(*glob_path.parts[:idx])
    glob_pattern = str(pathlib.Path(*glob_path.parts[idx:]))
    return list(glob_base.glob(glob_pattern))


def glob_variable(file_key, variable, pattern, file_vars=None):
    """Find all possible values of variable
    """
    import sys

    caller = sys._getframe(1)
    func_name = caller.f_code.co_name
    file_name = caller.f_code.co_filename
    line_num = caller.f_lineno
    log.dev(
        f"{file_name} ({line_num}) {func_name}: 'lib.files.glob_variable()' is deprecated. Use 'lib.config.files.glob_variable()' instead"
    )

    # Find available paths
    file_vars = dict() if file_vars is None else dict(file_vars)
    file_vars[variable] = "*"
    search_paths = glob_paths(file_key, file_vars)

    # Set up the regular expression
    re_vars = {**file_vars, variable: f"(?P<{variable}>__pattern__)"}
    path_pattern = str(path(file_key, file_vars=re_vars, default=".*")).replace("\\", "\\\\")
    for i in itertools.count():
        # Give unique names to each occurance of variable
        path_pattern = path_pattern.replace(f"<{variable}>", f"<{variable}__{i}>", 1)
        if f"<{variable}>" not in path_pattern:
            break
    re_pattern = re.compile(path_pattern.replace("__pattern__", pattern))

    # Find each match
    values = set()
    for search_path in search_paths:
        match = re_pattern.search(str(search_path))
        if match:
            matches = set(match.groupdict().values())
            if len(matches) > 1:
                log.warn(f"Found multiple values for {variable!r} in {search_path}: {', '.join(matches)}")
            values |= matches
    return values


def _log_file_open(file_path, description="", mode="r"):
    """Write a message to the log about a file being opened

    Args:
        file_path (Path/String):  The path to file being opened.
        description (String):     Description used for logging.
        mode (String):            Same as for the built-in open, usually 'r' or 'w'.
    """
    if description:
        description += " "

    mode_text = "Read {}from {}"
    if "w" in mode:
        mode_text = "Write {}to {}"
        if file_path.is_file():
            mode_text = "Overwrite {}on {}"
    if "a" in mode:
        mode_text = "Append {}to {}"
    log.debug(mode_text.format(description, file_path))
