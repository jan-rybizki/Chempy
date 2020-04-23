"""
Provide logged base Class.
"""

import io
import re
import copy
import os
import time
import logging
import contextlib
import inspect

from fortranfile import FortranReader
import utils

from datetime import datetime, timedelta
from byte2human import byte2human
from time2human import time2human
from version2human import version2human

class Timer(object):
    """
    A basic timer class.
    """
    def __init__(self, start = True):
        self.reset(start)
    def _get_time(self):
        return datetime.now()
    def start(self):
        """
        start timer
        """
        assert not self._running
        assert self._time_offset == timedelta()
        self._running = True
        self._time = self._get_time()
    def stop(self):
        """
        stop (pause) timer
        """
        assert self._running
        self._time_offset = self()
        self._running = False
    def reset(self, start = True):
        """
        reset timer
        """
        self._time_offset = timedelta()
        self._running = False
        if start:
            self.start()
    def resume(self):
        """
        resume timer
        """
        assert not self._running
        self._running = True
        self._time = self._get_time()
    def __call__(self):
        """
        get current time value
        """
        if self._running:
            return self._time_offset + self._get_time() - self._time
        else:
            return self._time_offset

class Timed(object):
    """
    Object to provide timers.

    Default timer 'name' is None.
    """
    def __init__(self,
                silent = True):
        """
        Constructor to set up timing.

        Parameters:
        silent = True - no output.

        Better to use setup_timed() instead.
        """
        self.setup_timer(silent)

    def setup_timer(self,
                    timer = None,
                    silent = True):
        """
        Set up one default timer 'None'
        """
        self._timers = dict()
        self._timers[None] = Timer()
        if timer is not None:
            self._timers[timer] = Timer()
        self._current_timer = timer

    def has_timer(self, name = None):
        """
        Return True if timer exists.
        """
        return name in self._timers

    def new_timer(self, name = None, start = True):
        """
        Add a new timer.  Raise if it timer exists.
        """
        assert not self._timers.haskey(name), 'timer {} already exisits'.format(name)
        self._timers[name] = Timer(start)

    def add_timer(self, name = None, start = True):
        """
        Add a new timer.  Reset if it timer exists.
        """
        self._timers[name] = Timer(start)

    def start_timer(self, name = None):
        """
        Start timer.  Reset if it timer exists.
        """
        if name in self._timers:
            self._timers[name].start
        else:
            self._timers[name] = Timer(True)
    def get_timer(self, name = None):
        """
        Retrive current time of timer
        """
        return self._timers[name]()
    def get_timer_human(self, name = None):
        """
        Retrive current time of timer in human-readable format (str).
        """
        return time2human(self._timers[name]())
    def finish_timer(self, name = None):
        """
        End timer and retrive current time of timer
        """
        return self._timers.pop(name)()
    def timer_is_running(self, name = None):
        """
        Retrive current time status
        """
        return self._timers[name]._running
    def stop_timer(self, name = None):
        """
        Stop given timer
        """
        self._timers[name].stop()
    def reset_timer(self, name = None):
        """
        Reset named timer
        """
        self._timers[name].reset()
    def resume_timer(self, name = None):
        """
        Start or resume given timer
        """
        self._timers[name].resume()
    def restart_timer(self, name = None):
        """
        restart/reinitialize named timer
        """
        self._timers[name].restart()
    def switch_timer(self, name = None):
        """
        Switch to new current timer
        """
        assert self._current_timer is not None
        self._current_timer.stop()
        self._current_timer = self._timers.setdefault(name, Timer())
        self._current_timer.resume()
    def set_timer(self, name = None):
        """
        Set a current timer.
        """
        self._current_timer = self._timers.setdefault(name, Timer())
    def replace_timer(self, name = None, timer = None):
        """
        Replace a time in data base.
        """
        assert timer is not None
        self._timers[name] = timer


class Logged(Timed):
    """
    Provide logging setup.
    """
    def __init__(self, **kwargs):
        """
        Constructor to set up logging.

        Parameters:
        silent = True - no output.

        Better to use setup_logger() instead.
        """
        kwargs.setdefault('silent', True)
        self.setup_logger(**kwargs)

    @contextlib.contextmanager
    def timeenv(self, **kwargs):
        """
        Context manager for logging with timing default.

        Specifically, you may specify a named timer that can be used
        otherwise as well and is returned in the with statement (timer
        name, not timer object).  This way you can have nested timers,
        should need be.
        """
        kw = dict(kwargs)
        message = kw.pop('message', None)
        timing  = kw.pop('timing' , True)
        timer   = kw.setdefault('timer' , None)
        self.setup_logger(**kw)
        yield timer
        self.close_logger(
            message = message,
            timing  = timing)

    @contextlib.contextmanager
    def logenv(self, **kwargs):
        """
        Context manager for logging.

        Set timing to False to disable printout.
        """
        kw = dict(kwargs)
        message = kw.pop('message', None)
        timing  = kw.pop('timing' , False)
        self.setup_logger(**kw)
        yield
        self.close_logger(
            message = message,
            timing  = timing,
            )

    # Since we can't be sure the __init__ method was called, this was
    # implemented as a property
    @property
    def logger_count(self):
        """
        Logger count level.
        """
        try:
            count = self._logger_count
        except AttributeError:
            count = 0
            self._logger_count = count
        return count
    @logger_count.setter
    def logger_count(self, count):
        self._logger_count = count
    @logger_count.deleter
    def logger_count(self):
        del self._logger_count

    def _info(self, message):
        """
        One-time logging of messages.
        """
        if self.logger_count <= 0:
            with self.logenv(silent = False,
                             timing = False):
                self.logger.info(message)
        else:
            self.logger.info(message)

    def setup_logger(self,
                     silent = None,
                     logfile = None,
                     level = None,
                     format = None,
                     timer = None,
                     filename = False,
                     pathname = False,
                     funcname = False,
                     linenumber = False,
                     process = False,
                     name = None,
                     classname = None,
                     ):
        """
        Set up logger for output.

        Parameters:
        silent = True - no output.

        Use:
        self.logger.info('My string')
        to log output data

        TODO:
        add setups for file handler and level
        """
        self.logger_count += 1
        self.logger_logfile = None
        if self.logger_count == 1:
            if not '_timers' in self.__dict__:
                self.setup_timer()
            self._log_timers = []
            if not silent in (None, True, False):
                level = silent
                silent = False
            if silent is None:
                if level is None:
                    try:
                        silent = self.logger_silent
                    except AttributeError:
                        silent = True
                else:
                    silent = False
            self.logger_silent = silent
            if level is None:
                try:
                    level = self.logger_level
                except AttributeError:
                    level = logging.INFO

            # TODO if current silence level - or file - is different
            # from root level: install new logger  - maybe we should
            # never recycle root logger ...

            self.logger_level = level
            self.logger = logging.getLogger(self.__class__.__name__)
            root_logger = logging.getLogger()
            if len(root_logger.handlers) == 0 and len(self.logger.handlers) == 0:
                if format is None:
                    info = ''
                    if name is not None:
                        info = name
                        if classname is True:
                            info = '%(name)s-' + info
                    elif classname is not False:
                        info = '%(name)s'
                    if funcname:
                        info += '.%(funcName)s'
                    if filename:
                        info = '%(filename)s.' + info
                    elif pathname:
                        info = '%(pathname)s.' + info
                    if linenumber:
                        info += ':%(lineno)d'
                    if process:
                        info = '%(process)d:' + info
                    formatter = logging.Formatter(' [{}] %(message)s'.format(info))
                elif format is 'UTC':
                    formatter = utils.UTCFormatter('%(asctime)s%(msecs)03d %(nameb)12-s %(levelname)s: %(message)s',
                                                   datefmt = '%Y%m%d%H%M%S')
                else:
                    raise Exception('Logger Format Not Recognized')
                if silent is True:
                    self.logger_handler = logging.NullHandler()
                elif logfile is not None:
                    self.logger_logfile = logfile
                    self.logger_handler = logging.FileHandler(logfile, 'w')
                else:
                    self.logger_handler = logging.StreamHandler()
                self.logger_handler.setFormatter(formatter)
                self.logger.addHandler(self.logger_handler)
                self.logger.setLevel(level)
            else:
                self.logger.setLevel(level)
                self.logger_handler = None
        else:
            self._log_timers += [(self._log_timer, self._timers[self._log_timer])]
        self.add_timer(timer)
        self._log_timer = timer

    @staticmethod
    def _format_timing(timing, time):
        """
        Format timing string.

        If it contains a {.*} expression, assume this is for use for
        formatting.
        """
        if timing is False:
            return None
        elif timing in (True, None):
            timing = 'Runtime:'
        stime = time2human(time)
        if len(re.findall(r'({.*})', timing)) == 1:
            s = timing.format(stime)
        else:
            s = '{:s} {:s}.'.format(timing, stime)
        return s

    def logger_timing(self,
                      timing = None,
                      timer = None,
                      finish = False):
        """
        Print logger timing info.
        """
        if timer is not None and timing is None:
            timing = timer + ':'
        timing = self._format_timing(timing, self.get_timer(timer))
        if timing is not None:
            self.logger.info(timing)
        if finish and timer is not None:
            self.finish_timer(timer)

    def close_logger(self,
                     timing = None,
                     message = None):
        """
        Reset logger if it was changed.
        """
        if timing is not None:
            self.logger_timing(timing = timing)
        if message is not None:
            self.logger.info(message)
        self.logger_count -= 1
        assert self.logger_count >=0, 'logger count is {}'.format(self.logger_count)
        if self.logger_count == 0:
            if self.logger_handler is not None:
                self.logger.removeHandler(self.logger_handler)
            if self.logger_logfile is not None:
                logging.shutdown()
        else:
            self.finish_timer(self._log_timer)
            self._log_timer, timer = self._log_timers.pop()
            self.replace_timer(name = self._log_timer, timer = timer)


    def close_timer(self,
                    timer = None,
                    timing = None):
        """
        End/close timer, return time, optionally log message (parameter 'timing').

        This should not be used to close logger timer(s).
        """
        try:
            log_timers = self._log_timers
        except:
            log_timers = []
        if self._timers[timer] in log_timers:
            self.logger.error('Cannot close logger timer {}'.format(timer))
            time = self.get_timer(timer)
            return time
        try:
            time = self.finish_timer(timer)
        except KeyError:
            time = 0.
        if timing is not None:
            self.logger.info(self._format_timing(timing, time))
        return time


# probably the following should go into their own class
    def logger_file_info(self, f):
        """
        Log file information.
        """
        if isinstance(f, FortranReader):
            filename = f.filename
            filesize = f.filesize
            self.logger.info('Loading {:s} ({:s})'.format(
                filename,
                byte2human(filesize)))
            if f.compressed:
                filesize = f.stat.st_size
                self.logger.info('Compressed file size: {:s} (modulo {:s} for gz)'.format(
                byte2human(filesize),
                byte2human(2**32)))
        elif isinstance(f, io.IOBase):
            stat = os.fstat(f.fileno())
            filename = f.name
            filesize = stat.st_size
            self.logger.info('Loading {:s} ({:s})'.format(
                filename,
                byte2human(filesize)))


    def logger_load_info(self,
                         nvers,
                         ncyc0,
                         ncyc1,
                         nmodels,
                         time = None):
        """
        Log information on load process and data.
        """
        if time is None:
           time = self.get_timer()
        self.logger.info('         version {:>9s}'.format(version2human(nvers)))
        self.logger.info('first model read {:>9d}'.format(int(ncyc0)))
        self.logger.info(' last model read {:>9d}'.format(int(ncyc1)))
        self.logger.info(' num models read {:>9d}'.format(int(nmodels)))
        self.logger.info('  data loaded in {:>9s}'.format(time2human(time)))

class logged(Logged, contextlib.ContextDecorator):
    def __init__(self, **kwargs):
        self.kw = dict(kwargs)
        self.message = self.kw.pop('message', None)
        self.timing  = self.kw.pop('timing' , False)

    def __enter__(self):
        self.setup_logger(**self.kw)

    def __exit__(self):
        self.close_logger(
            message = self.message,
            timing  = self.timing,
            )

class timed(Logged, contextlib.ContextDecorator):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.timer   = self.kw.setdefault('timer' , None)

def static_logging(func):
    func = staticmethod(func)
    class logged(Logged):
        """
        class to provied logging to sttaic functions
        """
        def __call__(self, *args, **kwargs):
            return func(*args, **kwargs)
    logged.__dict__.update(func.__dict__)
    logged.__dict__['method'] = func.__name__
    if func.__doc__ is not None:
        logged.__doc__ = func.__doc__ + '\n' + logged.__doc__
    logged.__name__ = logged.__name__
    logged.__module__ = getattr(func, '__module__')
    return logged
