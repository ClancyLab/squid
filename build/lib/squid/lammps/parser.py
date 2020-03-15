'''
The LMP_Parser code is code from the pizza.py toolkit
(www.cs.sandia.gov/~sjplimp/pizza.html) developed by Steve Plimpton
(sjplimp@sandia.gov) with some additional warning interspersed through
the thermo output.

- :class:`LMP_Parser`

------------

'''
# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.
# Imports and external programs
import re
import sys
import glob
from os import popen
# log tool

oneline = "Read LAMMPS log files and extract thermodynamic data"

docstr = '''
l = log("file1")                     read in one or more log files
l = log("log1 log2.gz")              can be gzipped
l = log("file*")                     wildcard expands to multiple files
l = log("log.lammps",0)              two args = store filename, but don't read

    incomplete and duplicate thermo entries are deleted

time = l.next()                      read new thermo info from file

    used with 2-argument constructor to allow reading thermo incrementally
    return time stamp of last thermo read
    return -1 if no new thermo since last read

nvec = l.nvec                        # of vectors of thermo info
nlen = l.nlen                        length of each vectors
names = l.names                      list of vector names
t,pe,... = l.get("Time","KE",...)    return one or more vectors of values
l.write("file.txt")          write all vectors to a file
l.write("file.txt","Time","PE",...)  write listed vectors to a file

    get and write allow abbreviated (uniquely) vector names
'''

# History
#   8/05, Steve Plimpton (SNL): original version
#   12/15, Yaset Acevedo: Allows for warnings interspersed through
#                         the thermo output

# ToDo list

# Variables
#   nvec = # of vectors
#   nlen = length of each vector
#   names = list of vector names
#   ptr = dictionary, key = name, value = index into data for which column
#   data[i][j] = 2d array of floats, i = 0 to # of entries, j = 0 to nvecs-1
#   style = style of LAMMPS log file, 1 = multi, 2 = one, 3 = gran
#   firststr = string that begins a thermo section in log file
#   increment = 1 if log file being read incrementally
#   eof = ptr into incremental file for where to start next read

try:
    tmp = PIZZA_GUNZIP
except NameError:
    PIZZA_GUNZIP = "gunzip"


# Class definition
class LMP_Parser(object):
    '''
    Class object to assist in parsing lammps outputs.

    **Parameters**

        list: *str*
            Path to the lammps log file that is to be parsed. Note,
            several files can be included in this string as long as
            they are separated by spaces.
        read_all: *int, optional*
            If this is set to 0, don't read in the whole file upon
            initialization. This lets you use the :func:`next`
            functionality.

    **Contains**

        nvec: *int*
            Number of vectors.
        nlen: *int*
            Length of each vector.
        names: *list, str*
            List of vector names.
        ptr: *dict*
            Dictionary corresponding the thermo keys to which column
            of the output they reside in. ptr[thermo_key] = which
            column this data is in
        data: *list, list, float*
            Raw data from file, organized into 2d array.
        style: *int*
            What style the LAMMPs log file is in. 1 = multi, 2 = one, 3 = gran
        firststr: *str*
            String that begins a thermo section in log file.
        increment: *int*
            1 if log file being read incrementally
        eof: *int*
            ptr into incremental file for where to start next read

    **Returns**

        This :class:`LMP_Parser` object.
    '''
    # --------------------------------------------------------------------

    def __init__(self, *list):
        self.nvec = 0
        self.names = []
        self.ptr = {}
        self.data = []
        # Stores when log file was last modified in seconds
        self.last_modified = 'Null'

        # flist = list of all log file names

        words = list[0].split()
        self.flist = []
        for word in words:
            self.flist += glob.glob(word)
        if len(self.flist) == 0 and len(list) == 1:
            raise Exception("No log file specified")

        if len(list) == 1:
            self.increment = 0
            self.read_all()
        else:
            if len(self.flist) > 1:
                raise Exception("Can only incrementally read one log file")
            self.increment = 1
            self.eof = 0

    # --------------------------------------------------------------------
    # read all thermo from all files
    '''
    Read all the thermo data from all the files.

    **Returns**

        None.
    '''

    def read_all(self):
        self.read_header(self.flist[0])
        if self.nvec == 0:
            raise Exception("Log file has no values")

        # read all files
        for file in self.flist:
            self.read_one(file)

        # sort entries by timestep, cull duplicates
        self.data.sort(self.compare)
        self.cull()
        self.nlen = len(self.data)

    # --------------------------------------------------------------------
    def next(self):
        '''
        Read the next line of thermo information from the file.
        Note, this is used when two arguments are passed during initialization.

        **Returns**

            timestep: *int*
                The timestep of the parsed thermo output.
        '''
        if not self.increment:
            raise Exception("cannot read incrementally")

        if self.nvec == 0:
            try:
                open(self.flist[0], 'r')
            except IOError:
                return -1
            except NameError:
                return -1
            self.read_header(self.flist[0])
            if self.nvec == 0:
                return -1

        self.eof = self.read_one(self.flist[0], self.eof)
        return int(self.data[-1][0])

    # --------------------------------------------------------------------
    def get(self, *keys):
        '''
        Read specific values from thermo output.

        **Parameters**

            keys: *str*
                Which thermo outputs you want by ID.  Not, this is as many
                requests as you want.  ex. l.get("Time", "KE", ...)

        **Returns**

            vecs: *list*
                Desired outputs.
        '''
        if len(keys) == 0:
            raise Exception("no log vectors specified")

        map = []
        for key in keys:
            if key in self.ptr:
                map.append(self.ptr[key])
            else:
                count = 0
                for i in range(self.nvec):
                    if self.names[i].find(key) == 0:
                        count += 1
                        index = i
                if count == 1:
                    map.append(index)
                else:
                    raise Exception("unique log vector %s not found" % key)

        vecs = []
        for i in range(len(keys)):
            vecs.append(self.nlen * [0])
            for j in range(self.nlen):
                vecs[i][j] = self.data[j][map[i]]

        if len(keys) == 1:
            return vecs[0]
        else:
            return vecs

    # --------------------------------------------------------------------
    def write(self, filename, *keys):
        '''
        Write parsed vectors to a file.

        **Parameters**

            filename: *str*
                The name of the file you want to dump all your outputs to.

            keys: *str, optional*
                Which specific vectors you want output to the file.
                ex. >>> l.write("file.txt","Time", "KE", ...)

        **Returns**

            None
        '''
        if len(keys):
            map = []
            for key in keys:
                if key in self.ptr:
                    map.append(self.ptr[key])
                else:
                    count = 0
                    for i in range(self.nvec):
                        if self.names[i].find(key) == 0:
                            count += 1
                            index = i
                    if count == 1:
                        map.append(index)
                    else:
                        raise Exception(
                            "unique log vector %s not found" % key)
        else:
            map = range(self.nvec)

        f = open(filename, "w")
        for i in range(self.nlen):
            for j in range(len(map)):
                print >>f, self.data[i][map[j]],
            print >>f
        f.close()

    # --------------------------------------------------------------------
    def compare(self, a, b):
        if a[0] < b[0]:
            return -1
        elif a[0] > b[0]:
            return 1
        else:
            return 0

    # --------------------------------------------------------------------
    def cull(self):
        i = 1
        while i < len(self.data):
            if self.data[i][0] == self.data[i - 1][0]:
                del self.data[i]
            else:
                i += 1

    # --------------------------------------------------------------------
    def read_header(self, file):
        str_multi = "----- Step"
        str_one = "Step "

        if file[-3:] == ".gz":
            txt = popen("%s -c %s" % (PIZZA_GUNZIP, file), 'r').read()
        else:
            txt = open(file).read()

        if txt.find(str_multi) >= 0:
            self.firststr = str_multi
            self.style = 1
        elif txt.find(str_one) >= 0:
            self.firststr = str_one
            self.style = 2
        else:
            return

        if self.style == 1:
            s1 = txt.find(self.firststr)
            s2 = txt.find("\n--", s1)
            pattern = "\s(\S*)\s*="
            keywords = re.findall(pattern, txt[s1:s2])
            keywords.insert(0, "Step")
            i = 0
            for keyword in keywords:
                self.names.append(keyword)
                self.ptr[keyword] = i
                i += 1

        else:
            s1 = txt.find(self.firststr)
            s2 = txt.find("\n", s1)
            line = txt[s1:s2]
            words = line.split()
            for i in range(len(words)):
                self.names.append(words[i])
                self.ptr[words[i]] = i

        self.nvec = len(self.names)

    # --------------------------------------------------------------------
    def read_one(self, *list):

        # if 2nd arg exists set file ptr to that value
        # read entire (rest of) file into txt

        file = list[0]
        if file[-3:] == ".gz":
            f = popen("%s -c %s" % (PIZZA_GUNZIP, file), 'rb')
        else:
            f = open(file, 'rb')

        if len(list) == 2:
            f.seek(list[1])
        txt = f.read()
        if file[-3:] == ".gz":
            eof = 0
        else:
            eof = f.tell()
        f.close()

        start = last = 0
        while not last:

            # chunk = contiguous set of thermo entries (line or multi-line)
            # s1 = 1st char on 1st line of chunk
            # s2 = 1st char on line after chunk
            # set last = 1 if this is last chunk in file, leave 0 otherwise
            # set start = position in file to start looking for next chunk
            # rewind eof if final entry is incomplete

            s1 = txt.find(self.firststr, start)
            s2 = txt.find("Loop time of", start + 1)

            # ACE EDIT***************************
            s3 = -1
            # Check if there is an warning message on the line prior
            # to "Loop time of"
            s3 = txt.find("WARNING", s1, s2)
            # Change s2 to reflect new position if appropriate
            if s3 < s2 and s3 > s1:
                s2 = s3
            # ***********************************

            if s2 == -1:
                s2 = txt.find("ERROR", start + 1)
            if s1 >= 0 and s2 >= 0 and s1 < s2:
                if self.style == 2:
                    s1 = txt.find("\n", s1) + 1
            elif s1 >= 0 and s2 >= 0 and s2 < s1:
                s1 = 0
            elif s1 == -1 and s2 >= 0:
                last = 1
                s1 = 0
            elif s1 >= 0 and s2 == -1:
                last = 1
                if self.style == 1:
                    s2 = txt.rfind("\n--", s1) + 1
                else:
                    s1 = txt.find("\n", s1) + 1
                    s2 = txt.rfind("\n", s1) + 1
                    eof -= len(txt) - s2
            elif s1 == -1 and s2 == -1:     # found neither
                                            # could be end-of-file section
                                            # or entire read was one chunk

                if txt.find("Loop time of", start) == start:
                    eof -= len(txt) - start
                    break
                if txt.find("ERROR", start) == start:
                    eof -= len(txt) - start
                    break

                # entire read is a chunk
                last = 1
                s1 = 0
                if self.style == 1:
                    s2 = txt.rfind("\n--", s1) + 1
                else:
                    s2 = txt.rfind("\n", s1) + 1
                    eof -= len(txt) - s2
                    if s1 == s2:
                        break

            chunk = txt[s1:s2 - 1]
            start = s2

            # ACE EDIT***************************
            # Correct start position if s3 was utilized
            if s3 != -1:
                s3 = txt.find("Loop time of", start + 1)
                start = s3
            # ***********************************

            # split chunk into entries
            # parse each entry for numeric fields, append to data
            if self.style == 1:
                sections = chunk.split("\n--")
                pat1 = re.compile("Step\s*(\S*)\s")
                pat2 = re.compile("=\s*(\S*)")
                for section in sections:
                    word1 = [re.search(pat1, section).group(1)]
                    word2 = re.findall(pat2, section)
                    words = word1 + word2
                    self.data.append(map(float, words))
            else:
                lines = chunk.split("\n")
                for line in lines:
                    words = line.split()
                    self.data.append(map(float, words))

            # print last timestep of chunk
            sys.stdout.flush()

        return eof
