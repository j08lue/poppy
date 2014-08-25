import subprocess
import traceback

def get_ulimitn(default=1e3):
    """Try to get the maximum number of open files on the system. Works with Unix."""
    try:
        maxn = int(subprocess.check_output(['ulimit -n'],shell=True))
    except:
        maxn = default
        traceback.print_exc()
    return maxn
