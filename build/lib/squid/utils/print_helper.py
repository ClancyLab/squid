# -*- coding: utf-8 -*-

# System imports
import re
import sys
from subprocess import Popen, PIPE
# Squid imports
from squid import constants
from squid.utils.cast import is_numeric


def color_set(s, c):
    '''
    Colourize a string for linux terminal output.

    **Parameters**

        s: *str*
            String to be formatted.
        c: *str*
            Colour or format for the string, found in constants.COLOUR.

    **Returns**

        s: *str*
            Coloured or formatted string.
    '''
    return constants.COLOR[c] + str(s) + constants.COLOR['ENDC']


colour_set = color_set


def strip_color(s):
    '''
    Remove colour and/or string formatting due to linux escape sequences.

    **Parameters**

        s: *str*
            String to strip formatting from.

    **Returns**

        s: *str*
            Unformatted string.
    '''
    for c in constants.COLOUR:
        col = constants.COLOUR[c]
        while col in s:
            s = s.replace(col, '')
    return s


strip_colour = strip_color


def spaced_print(sOut, delim=['\t', ' '], buf=4):
    '''
    Given a list of strings, or a string with new lines, this will
    reformat the string with spaces to split columns.  Note, this
    only works if there are no headers to the input string/list of strings.

    **Parameters**

        sOut: *str* or *list, str*
            String/list of strings to be formatted.
        delim: *list, str*
            List of delimiters in the input strings.
        buf: *int*
            The number of spaces to have between columns.

    **Returns**

        spaced_s: *str*
            Appropriately spaced output string.
    '''
    s_len = []
    if type(sOut) == str:
        sOut = sOut.split('\n')
    if type(delim) == list:
        delim = ''.join([d + '|' for d in delim])[:-1]
    # Get the longest length in the column
    for i, s in enumerate(sOut):
        s = re.split(delim, s)
        for j, sss in enumerate(s):
            ss = strip_colour(sss)
            try:
                # This makes the part of the list for each column the
                # longest length
                s_len[j] = len(ss) if len(ss) > s_len[j] else s_len[j]
            except:
                # If we are creating a new column this happens
                s_len.append(len(ss))
    # Now we add a buffer to each column
    for i in range(len(s_len)):
        s_len[i] += buf

    # Compile string output
    for i, s in enumerate(sOut):
        s = re.split(delim, s)
        for j, ss in enumerate(s):
            s[j] = ss + ''.join([' '] * (s_len[j] - len(strip_colour(ss))))
        sOut[i] = ''.join(s)

    return '\n'.join(sOut)


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1,
                     length=20, fill='+', buf=None, pad=False):
    '''
    NOTE! THIS IS COPIED FROM STACK OVERFLOW (with minor changes),
    USER Greenstick
    Link: https://stackoverflow.com/a/34325723

    Call in a loop to create terminal progress bar.

    **Parameters**

        iteration: *int*
            Current iteration.
        total: *int*
            Total number of iterations.
        prefix: *str, optional*
            Prefix for the loading bar.
        suffix: *str, optional*
            Suffix for the loading bar.
        decimals: *int, optional*
            Positive number of decimals in percent complete
        length: *int, optional*
            Character length of the loading bar.
        fill: *str, optional*
            Bar fill character.
        pad: *bool, optional*
            Whether to pad the right side with spaces until terminal width.
    '''
    if buf is not None:
        if not hasattr(printProgressBar, "buf"):
            setattr(printProgressBar, "buf", buf)
        elif printProgressBar.buf < 0:
            printProgressBar.buf = buf
        else:
            printProgressBar.buf -= 1
            return

    assert "stty" not in prefix and "stty" not in suffix,\
        "Don't have 'stty' in prefix or suffix."

    assert is_numeric(total),\
        "Error - total is not a numerical value."
    assert float(total) > 0,\
        "Error - total is not greater than 0."
    assert is_numeric(iteration),\
        "Error - iteration is not a numerical value."
    assert float(iteration) >= 0,\
        "Error - iteration is less than 0."
    assert float(iteration) <= total,\
        "Error - iteration is larger than total."

    percent = ("{0:." + str(decimals) + "f}").format(
        100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    bar = '\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix)
    if pad:
        try:
            p = Popen("stty size".split(), stdout=PIPE, stderr=PIPE)
            cols, stderr = p.communicate()
            if stderr.strip() is not "":
                cols = 300
            else:
                cols = int(cols.strip().split()[-1])
        except IndexError:
            # This is a backup in the case that this code is
            # run on the queue or something.
            cols = 300
        cols = max([cols, len(bar) + 1])  # Ensure we are always long enough
        bar = bar + "".join([" " for i in range(cols - len(bar))])
    sys.stdout.write(bar)
    # Print New Line on Complete
    if iteration == total:
        sys.stdout.write('\n')

    sys.stdout.flush()


def bytes2human(n):
    '''
    Convert n bytes (as integer) to a human readable string.  Code was found
    online at activestate (see references).

    **Parameters**

        n: *int*
            The number of bytes.

    **Returns**

        n_in_str: *str*
            The bytes in string format.

    **References**

        - http://code.activestate.com/recipes/578019
    '''
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n
