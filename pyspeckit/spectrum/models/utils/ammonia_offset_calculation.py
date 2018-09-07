from __future__ import print_function
import re
import numpy as np
from astropy import units as u
from astropy import constants
from astroquery.splatalogue import Splatalogue, utils

# Query splatalogue, keeping all of the line strength columns
# Both Lovas and CDMS/JPL can be used
nh3 = Splatalogue.query_lines(20*u.GHz, 40*u.GHz, chemical_name=' NH3 ',
                              show_upper_degeneracy=True,
                              line_strengths=['ls1','ls2','ls3','ls4'])

numbers = {1:'one',
           2:'two',
           3:'three',
           4:'four',
           5:'five',
           6:'six',
           7:'seven',
           8:'eight',
           9:'nine',}

tbls = {}


for line in (1,2,3,4,5,6,7,8,9):

    reline = re.compile('^{n}\({n}\).*-{n}'.format(n=line))

    tbl = utils.minimize_table(nh3[np.array([bool(reline.search(x))
                                                  if bool(x) else False
                                                  for x in
                                                  nh3['Resolved QNs']],
                                                 dtype='bool')],
                                    columns=['Species', 'Chemical Name', 'Resolved QNs',
                                             'Freq-GHz', 'Meas Freq-GHz',
                                             'Log<sub>10</sub> (A<sub>ij</sub>)',
                                             'CDMS/JPL Intensity',
                                             'Linelist',
                                             'E_U (K)', 'Upper State Degeneracy'])

    if len(tbl) == 0:
        pass

    # Select only TopModel lines from CDMS/JPL
    tbls[line] = tbl[tbl['Linelist'] == 'TopModel']

for par in ('tau_wts','voff_lines','aval','freq'):
    print(par)
    for line in (1,2,3,4,5,6,7,8): # 9 not available

        tbl = tbls[line]

        degeneracyline = tbl['Upper State Degeneracy']
        intensityline = 10**tbl['CDMSJPL_Intensity']

        main = np.argmax(intensityline)

        centerline = tbl['Freq'][main]
        voff_linesline = np.array((centerline-tbl['Freq'])/centerline) * constants.c

        aval = (10**tbl['log10_Aij']).sum()
        weightline = intensityline/intensityline.sum()

        if par == 'freq':
            print("'{n}{n}': {f},".format(n=numbers[line], f=centerline))
        elif par == 'voff_lines':
            print("'{n}{n}': [{v}],".format(n=numbers[line],
                                            v=", ".join(str(x)
                                                        for x in voff_linesline.to(u.km/u.s).value)))
        elif par == 'tau_wts':
            #print "'{n}{n}': {d},".format(n=numbers[line], d=np.array(degeneracyline))
            print("'{n}{n}': [{d}],".format(n=numbers[line],
                                            d=", ".join(str(x) for x in weightline)))
        elif par == 'aval':
            print("'{n}{n}': {d:e},".format(n=numbers[line], d=aval))
