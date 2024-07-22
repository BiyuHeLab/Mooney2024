import logging
import sys

logger = logging.getLogger("pwmv_ecog")
streamhdlr = logging.StreamHandler(stream=sys.stdout)
streamhdlr.setLevel(logging.INFO)
logger.addHandler(streamhdlr)
logger.setLevel(logging.DEBUG)

MAXWORKERS = 22
