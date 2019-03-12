"""Where library module for timing running time of functions and code blocks

Description:
------------

The `lib.timer` can be used to log running time of functions and general code blocks. Typically, you will import the
`timer`-class from within the module::

    from where.lib.timer import timer

The timer can then be used in three different ways:

1. As a decorator to time one function::

    @timer('The time to execute some_function was', unit='minutes')
    def some_function(some_argument, some_other_argument=some_value):
        pass

2. As a context manager together with `with` to time a code block::

    with timer('Finish doing stuff in', logger=log.debug) as t:
        do_something()
        do_something_else()

3. With explicit `start`- and `end`-statements::

    t = timer()
    t.start()
    do_something()
    do_something_else()
    t.end()

As can be seen in the examples above, `timer()` may be called with several optional parameters, including the text to
report when the timer ends, which units to show the time in and which logger is used to report the timing. See
`timer.__init__` for more details.

"""

# Standard library imports
from contextlib import ContextDecorator
import time
import sys

# External library imports

# Where imports
from where.lib import log
from where.lib.unit import Unit


class timer(ContextDecorator):
    """Class for timing running time of functions and code blocks.
    """

    def __init__(self, text="Elapsed time:", unit=None, logger=log.time):
        """Set up a new timer

        The text to be shown when logging the timer can be customized. Typically, the value of the timer will be added
        at the end of the string (e.g. 'Elapsed time: 0.1234 seconds'). However, this can be customized by adding a
        '{}' to the text. For example `text='Used {} to run the code'` will produce something like 'Used 0.1234 seconds
        to run the code'.

        Args:
            text (String):      Text used when logging the timer (see above).
            unit (String):      Unit used for logging the timer (Default is seconds).
            logger (Function):  Function used to do the logging.
        """
        super().__init__()
        self._start = None
        self._end = None
        self.text = text if "{}" in text else (text + " {}").strip()
        self.unit_name = "seconds" if unit is None else unit
        self.unit_factor = 1 if unit is None else Unit("seconds", unit)
        self.logger = logger

        # Use midgard instead
        caller = sys._getframe(1)
        func_name = caller.f_code.co_name
        file_name = caller.f_code.co_filename
        line_num = caller.f_lineno
        log.dev(f"{file_name} ({line_num}) {func_name}: where.lib.timer is deprecated, use midgard.dev.timer instead")

    @staticmethod
    def timer():
        """Get current value of timer

        Using the built-in `time.perf_counter` to do the timing.

        Returns:
            Float:  Current value of timer.
        """
        return time.perf_counter()

    def start(self):
        """Start the timer
        """
        self._start = self.timer()
        self._end = None

    def end(self):
        """End the timer and log the time elapsed

        Returns:
            Float:  The time elapsed in terms of `self.unit_name` (default seconds).
        """
        self._end = self.timer()
        time_elapsed = self.elapsed()
        self._start = None

        return time_elapsed

    def elapsed(self):
        """Log the time elapsed

        Can be used explicitly to log the time since a timer started without ending the timer.

        Returns:
            Float:  The time elapsed in terms of `self.unit_name` (default seconds).
        """
        if self._start is None:
            log.warn(f"The timer is not running. See `help({self.__module__})` for information on how to start one.")
            return

        timer_end = self.timer() if self._end is None else self._end
        time_elapsed = self.unit_factor * (timer_end - self._start)
        self._log(time_elapsed)

        return time_elapsed

    def _log(self, time_elapsed):
        """Do the actual logging of elapsed time

        Args:
            time_elapsed (Float):  The time elapsed in terms of `self.unit_name` (default seconds).
        """
        time_text = f"{time_elapsed:.4f} {self.unit_name}"
        if self.logger:
            self.logger(self.text.format(time_text))

    def __enter__(self):
        """Start the timer as a context manager
        """
        self.start()

    def __exit__(self, *exc):
        """End the timer and log the time elapsed as a context manager
        """
        self.end()
