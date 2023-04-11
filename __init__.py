import os
import time

__version__ = '1.0.0'
PKG_DIR = os.path.abspath(os.path.dirname(__file__))

# this is intended to provide a per-script-execution identifier to help
# troubleshoot serverside issues.
TIME_ID = int(time.time())
