#logConfig.py
#Description:   log file manager for HEAT
#Engineer:      T Looby
#Date:          20240321
import os
import logging
from logging.handlers import RotatingFileHandler


def setup_logging(logfile_path=None, level=logging.INFO, format='%(message)s'):
    """
    sets up logger.  
    enables future calls of this function to change the log file path.

    user can change the log file location like this:
        from logConfig import setup_logging
        setup_logging(logfile_path=<newLogFileName>)

    """
    # Disable werkzeug logging if using Flask
    logFlask = logging.getLogger('werkzeug')
    logFlask.disabled = True

    if logfile_path is None:
        logfile_path = os.getenv("logFile", "HEATlog.txt")

    if os.path.isfile(logfile_path):
        print("Deleting old logfile...")
        os.remove(logfile_path)

    # Clear existing handlers
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        if isinstance(handler, logging.FileHandler):
            handler.close()  # Close the handler to flush and release the file
            root_logger.removeHandler(handler)

    # Set up new handler
    file_handler = RotatingFileHandler(logfile_path, maxBytes=1024*1024*5, backupCount=5)
    file_handler.setFormatter(logging.Formatter(format))

    root_logger.setLevel(level)
    root_logger.addHandler(file_handler)
    return