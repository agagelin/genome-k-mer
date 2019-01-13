"""
ASCII codes to color outputs.
"""
import re

def custom_len(values):
    """
    Function to use in ProgressBar if using ascii codes.
    """
    pattern = r"\033\[\d*m"
    ascii_char = re.compile(pattern)
    values = ascii_char.sub("", values)
    return len(values)

# -- FOREGROUND -- #

# Dark foregrounds
FG_BLACK = "\033[30m"
FG_RED = "\033[31m"
FG_GREEN = "\033[32m"
FG_YELLOW = "\033[33m"
FG_BLUE = "\033[34m"
FG_MAGENTA = "\033[35m"
FG_CYAN = "\033[36m"
FG_WHITE = "\033[37m"

# Light foregrounds
FG_GRAY = "\033[90m"
LFG_RED = "\033[91m"
LFG_GREEN = "\033[92m"
LFG_YELLOW = "\033[93m"
LFG_BLUE = "\033[94m"
LFG_MAGENTA = "\033[95m"
LFG_CYAN = "\033[96m"
LFG_WHITE = "\033[97m"

# -- BACKGROUND -- #

# Dark backgrounds
BG_BLACK = "\033[40m"
BG_RED = "\033[41m"
BG_GREEN = "\033[42m"
BG_YELLOW = "\033[43m"
BG_BLUE = "\033[44m"
BG_MAGENTA = "\033[45m"
BG_CYAN = "\033[46m"
BG_WHITE = "\033[47m"

# Light backgrounds
BG_GRAY = "\033[100m"
LBG_RED = "\033[101m"
LBG_GREEN = "\033[102m"
LBG_YELLOW = "\033[103m"
LBG_BLUE = "\033[104m"
LBG_MAGENTA = "\033[105m"
LBG_CYAN = "\033[106m"
LBG_WHITE = "\033[107m"

# -- STOP -- #
END_COLOR = "\033[0m"


