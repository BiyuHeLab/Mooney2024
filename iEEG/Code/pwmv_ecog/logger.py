from datetime import datetime
import logging
from pathlib import Path

import mne

logger = logging.getLogger(__name__)


def set_top_FileHandler(fHdlr: logging.FileHandler, setMNE=True):
    """
    Add fHdlr to the library top logger. Optionally direct MNE to log to the same file.

    Args:
        fHdlr (logging.FileHandler): FileHandler to use.
        setMNE (bool, optional): Whether to set mne.set_log_file to fHdlr.basename. Defaults to
            True.
    """
    logging.getLogger("pwmv_ecog").addHandler(fHdlr)
    if setMNE:
        set_MNE_logfile(fHdlr.baseFilename)


def set_MNE_logfile(logFName: str | Path):
    """Set MNE to log to a file.

    Args:
        logFName (str | Path): File path at which to create the log file.
    """
    logger.debug(f"Setting MNE to log to {logFName}.")
    mne.set_log_file(
        logFName,
        output_format="%(asctime)s - %(levelname)s@%(name)s: %(message)s",
        overwrite=False,
    )


def get_dated_FileHandler(logBasename: str | Path, **kwargs) -> logging.FileHandler:
    """Set up logging to a log file on disk. Convenience wrapper for get_FileHandler to datestamp
    the log filename.

    Keyword arguments passed through to get_FileHandler.

    Args:
        log_fname (str | Path): Full path including base log filename at which to create the log
            file. A timestamp suffix will be added before any file extension.

    Returns:
        logging.FileHandler: The set up handler.
    """
    logBasename = Path(logBasename)
    logFname = logBasename.with_stem(logBasename.stem + f"-{datetime.utcnow()}".replace(" ", "_"))
    return get_FileHandler(logFname, **kwargs)


def get_FileHandler(log_fname, debug_filelogging=True, overwrite=False) -> logging.FileHandler:
    """Set up logging to a log file on disk.

    Args:
        log_fname (str): File path at which to create the log file.
        debug_filelogging (bool, optional): If True(default), emit maximal messages to
            log. When False, emits one step down, omitting many large log records mostly
            used for debugging.
        overwrite (bool, optional): Determines whether to overwrite any preexisting log
            file. Default False.

    Returns:
        logging.FileHandler: The set up handler.

    """
    Path(Path(log_fname).parent).mkdir(parents=True, exist_ok=True)
    logfile_handler = logging.FileHandler(log_fname, mode="w" if overwrite else "a")
    logfile_handler.setLevel(logging.DEBUG if debug_filelogging else logging.INFO)
    logfile_format = logging.Formatter(
        "{asctime} - {levelname:8}@{name}: {message}", style="{", datefmt="%Y%m%d %H%M%S"
    )
    logfile_handler.setFormatter(logfile_format)
    return logfile_handler
