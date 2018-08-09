# System imports
import os
import time
import getpass
from subprocess import Popen, PIPE

# Squid imports
from squid import print_helper

KNOWN_USERS = {"ns728": "Nikita",
               "hch54": "Henry",
               "bas348": "Blaire",
               "sb2326": "Vineeth",
               "rfh66": "Ryan",
               "yma3": "Ace",
               "te79": "Taha",
               "afh72": "Angela",
               "aec253": "Aron",
               "jds429": "Jonathon",
               "mr937": "Mardochee",
               "jy634": "Jee Won",
               "awr66": "Andrew",
               "sh864": "Spencer",
               "pbd44": "Philippe",
               "jw598": "Jingyang",
               "hj382": "Haili",
               "gmc225": "Greg",
               "ovr4": "Seun",
               "msm329": "Mia",
               "wg257": "Leo"}
FLIPPED_USERS = {v: k for k, v in KNOWN_USERS.items()}


class bc:
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'


def get_colour(usr, lvl, USER):
    if a[0] == USER:
        return bc.BLUE
    if a[0].strip() in FLIPPED_USERS and FLIPPED_USERS[a[0].strip()] == USER:
        return bc.BLUE
    if lvl < 20:
        return bc.GREEN
    if lvl < 70:
        return bc.YELLOW
    return bc.RED


def get_colour_2(lvl):
    if lvl < 200:
        return bc.GREEN
    if lvl < 700:
        return bc.YELLOW
    return bc.RED


USER = getpass.getuser()
while 1:
    # Get input from jlist as a string
    p = Popen(['jlist', '-all'], stdout=PIPE)
    output = p.stdout.read().split('\n')[4:-2]

    # Get user data
    usr_cnt = {}
    for d in output:
        usr = d.split()[2]
        if usr in usr_cnt:
            usr_cnt[usr][0] += 1
        else:
            usr_cnt[usr] = [1, 0, 0]
        if "RUNNING" in d:
            usr_cnt[usr][1] += 1
        else:
            usr_cnt[usr][2] += 1

    # Map to list for sort
    out_dat = []
    for a in usr_cnt:
        out_dat.append((a, usr_cnt[a][0], usr_cnt[a][1], usr_cnt[a][2]))
    out_dat = sorted(out_dat, key=lambda a: a[1])[::-1]

    # Convert known usernames
    for i, a in enumerate(out_dat):
        if a[0].strip() in KNOWN_USERS:
            out_dat[i] = (KNOWN_USERS[a[0].strip()], a[1], a[2], a[3])

    sOut = 'User\tTotal\tRunning\tPending \n'
    # Old output for dictionary
    for a in out_dat:
        sOut += (a[0] +
                 '\t' +
                 get_colour(a[0], a[1], USER) +
                 str(a[1]) +
                 bc.ENDC +
                 '\t' +
                 get_colour(a[0], a[2], USER) +
                 str(a[2]) +
                 bc.ENDC +
                 '\t' +
                 get_colour(a[0], a[3], USER) +
                 str(a[3]) +
                 bc.ENDC +
                 '\n')

    a, b, c = sum([x[1] for x in out_dat]), sum([x[2] for x in out_dat]), sum([x[3] for x in out_dat])

    sOut += 'Total:\t%s%d%s\t%s%d%s\t%s%d%s \n' % (get_colour_2(a), a, bc.ENDC, get_colour_2(b), b, bc.ENDC, get_colour_2(c), c, bc.ENDC)

    sOut = print_helper.spaced_print(sOut, ['\t'], buf=4).split('\n')

    while len(sOut) > 0 and sOut[-1].strip() == '':
        sOut = sOut[:-1]

    column_compare = sOut[1].split()
    sOut = '\n'.join([sOut[0][:-5]] +
                     [''.join(['-'] * (len(sOut[0])))] +
                     sOut[1:-1] +
                     [''.join(['-'] * (len(sOut[0])))] +
                     [sOut[-1]])

    os.system('echo "' + sOut + '"')

    time.sleep(20)
    os.system('clear')
