#logConfig.py
#Description:   log file manager for HEAT
#Engineer:      T Looby
#Date:          20240321
import os
import logging
import logging.handlers

# Global variable to hold the reference to the memory handler
BUFFER_HANDLER = None

def setup_logging(logfile_path=None, level=logging.INFO, format='%(message)s'):
    """
    Sets up logger.
    
    Mode 1 (Startup/Buffering): 
        If logfile_path is None, logs are held in RAM (MemoryHandler).
        This prevents race conditions during parallel startup.
        
    Mode 2 (File Logging):
        If logfile_path is provided, logs are written to disk.
        If a buffer exists from Mode 1, it is flushed to this new file.

    Usage:
        # 1. Startup (Buffers logs)
        setup_logging(logfile_path=None)
        
        # 2. Later, once unique path is known (Flushes buffer -> File)
        setup_logging(logfile_path='/path/to/unique/job/heat.log')
    """
    global BUFFER_HANDLER

    # Disable werkzeug logging if using Flask
    logFlask = logging.getLogger('werkzeug')
    logFlask.disabled = True
    
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    formatter = logging.Formatter(format)

    # --- SCENARIO 1: BUFFERING MODE (No path provided) ---
    if logfile_path is None:
        # Only set up buffer if we don't already have one and aren't already logging to file
        if BUFFER_HANDLER is None and not any(isinstance(h, logging.FileHandler) for h in root_logger.handlers):
            # Capacity: 10,000 lines. target=None means "hold in RAM"
            BUFFER_HANDLER = logging.handlers.MemoryHandler(capacity=10000, target=None)
            BUFFER_HANDLER.setFormatter(formatter)
            root_logger.addHandler(BUFFER_HANDLER)
            print("Logger initialized in BUFFER mode (holding logs in RAM).")
        return

    # --- SCENARIO 2: FILE LOGGING MODE (Path provided) ---
    # If we are here, we are switching to a real file.
    # Clean up old FileHandlers
    if os.path.isfile(logfile_path):
        try:
            print("Deleting old logfile...")
            os.remove(logfile_path)
        except OSError:
            pass # Handle case where file might be locked or race condition persists slightly

    for handler in root_logger.handlers[:]:
        if isinstance(handler, logging.FileHandler):
            handler.close()
            root_logger.removeHandler(handler)

    # Create the new File Handler
    # Try RotatingFileHandler (POSIX) -> FileHandler (S3) -> stderr (last resort)
    file_handler = None
    try:
        file_handler = logging.handlers.RotatingFileHandler(
            logfile_path,
            maxBytes=1024*1024*5,
            backupCount=5
        )
    except OSError:
        try:
            file_handler = logging.FileHandler(logfile_path, mode='w')
        except OSError:
            print(f"WARNING: Cannot write log file to {logfile_path}. Logging to stderr only.")

    if file_handler is not None:
        file_handler.setFormatter(formatter)

    # FLUSH BUFFER
    # If we were buffering, dump everything into the file handler (or discard if no file)
    if BUFFER_HANDLER is not None:
        if file_handler is not None:
            BUFFER_HANDLER.setTarget(file_handler)
            BUFFER_HANDLER.flush()
            print(f"Logger buffer flushed to {logfile_path}")
        root_logger.removeHandler(BUFFER_HANDLER)
        BUFFER_HANDLER.close()
        BUFFER_HANDLER = None

    # Attach file handler if we got one, otherwise add stderr so logs aren't lost
    if file_handler is not None:
        root_logger.addHandler(file_handler)
    elif not any(isinstance(h, logging.StreamHandler) for h in root_logger.handlers):
        stderr_handler = logging.StreamHandler()
        stderr_handler.setFormatter(formatter)
        root_logger.addHandler(stderr_handler)
    return

def get_current_log_path():
    """
    Returns the absolute path of the current log file.
    Returns None if we are currently buffering or only logging to console.
    """
    root_logger = logging.getLogger()
    
    for handler in root_logger.handlers:
        # Check for FileHandler (RotatingFileHandler inherits from this too)
        if isinstance(handler, logging.FileHandler):
            return handler.baseFilename
            
    return None