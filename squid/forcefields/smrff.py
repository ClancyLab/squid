def _remove_comments(s):
    '''
    A simple function to remove comments in a string.  Note, we don't care
    about whitespace, so it also removes that.

    **Parameters**

        s: *str*
            A string to remove comments from.

    **Returns**

        s_no_comments: *str*
            A string without comments.
    '''
    s = s.strip()
    if s == '':
        return ''
    elif '#' in s and s[0] != '#':
        return s[:s.index('#')]
    elif s[0] == '#':
        return ''
    return s


def parse_pfile(fname):
    '''
    This function will, given a smrff parameter file, will parse it by
    removing comments, trailing whitespaces, and empty lines.

    **Parameters**

        fname: *str*
            The name of the parameter file to be parsed.

    **Returns**

        parsed: *str*
            A parsed string of said parameter file.
    '''
    # Parse the input file to clean out comments, empty lines, and
    # trailing whitespace
    parsed = map(
        lambda x: x.strip(),
        open(fname, 'r').read().strip().split("\n"))
    parsed = map(_remove_comments, parsed)
    return '\n'.join([r for r in parsed if r != ''])
