# -*- coding: utf-8 -*-
"""
The Linux Helper module contains functionality to aid linux users to automate some tasks.

- :func:`color_set`
- :func:`colour_set`
- :func:`printProgressBar`
- :func:`strip_color`
- :func:`strip_colour`
- :func:`spaced_print`

------------

"""

# System imports
import os
import re
import sys
# Squid imports
import constants


def color_set(s, c):
    """
    Colourize a string for linux terminal output.

    **Parameters**

        s: *str*
            String to be formatted.
        c: *str*
            Colour or format for the string, found in constants.COLOUR.

    **Returns**

        s: *str*
            Coloured or formatted string.
    """
    return constants.COLOR[c] + str(s) + constants.COLOR['ENDC']


colour_set = color_set


def strip_color(s):
    """
    Remove colour and/or string formatting due to linux escape sequences.

    **Parameters**

        s: *str*
            String to strip formatting from.

    **Returns**

        s: *str*
            Unformatted string.
    """
    for c in constants.COLOUR:
        col = constants.COLOUR[c]
        while col in s:
            s = s.replace(col, '')
    return s


strip_colour = strip_color


def spaced_print(sOut, delim=['\t', ' '], buf=4):
    """
    Given a list of strings, or a string with new lines, this will reformat the string with spaces
    to split columns.  Note, this only works if there are no headers to the input string/list of strings.

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
    """
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
                # This makes the part of the list for each column the longest length
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


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=20, fill='+', buf=None, pad=False):
    """
    NOTE! THIS IS COPIED FROM STACK OVERFLOW (with minor changes), USER Greenstick
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
    """
    if buf is not None:
        if not hasattr(printProgressBar, "buf"):
            setattr(printProgressBar, "buf", buf)
        elif printProgressBar.buf < 0:
            printProgressBar.buf = buf
        else:
            printProgressBar.buf -= 1
            return

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    bar = '\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix)
    if pad:
        try:
            cols = int(os.popen('stty size', 'r').read().split()[-1])
        except IndexError:
            cols = 300  # This is a backup in the case that this code is run on the queue or something.
        cols = max([cols, len(bar) + 1])  # Ensure we are always long enough
        bar = bar + "".join([" " for i in range(cols - len(bar))])
    sys.stdout.write(bar)
    # Print New Line on Complete
    if iteration == total:
        sys.stdout.write('\n')

    sys.stdout.flush()
