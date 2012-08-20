interactive=False
import subprocess
versnum = subprocess.Popen(["hg","id","--num"],stdout=subprocess.PIPE).communicate()[0].strip().strip("+")
savedir = "tests_%s/" % versnum
import os
if not os.path.exists(savedir):
    os.mkdir(savedir)

# This step became necessary on August 17th, 2012.  It was never necessary before then.
# It doesn't make any sense.  I didn't make ANY changes!  WHY did the behavior change?!
import numpy as np
np.seterr(all='ignore')

from matplotlib.pyplot import ion,ioff
if interactive:
    ion()
else:
    ioff()

print "*****test_fits.py*****"
execfile('test_fits.py',{'interactive':interactive,'savedir':savedir})
print "*****test_hr2421.py*****"
execfile('test_hr2421.py',{'interactive':interactive,'savedir':savedir})
print "*****test_nh3.py*****"
execfile('test_nh3.py',{'interactive':interactive,'savedir':savedir})
print "*****test_sdss.py*****"
execfile('test_sdss.py',{'interactive':interactive,'savedir':savedir})
print "*****test_txt.py*****"
execfile('test_txt.py',{'interactive':interactive,'savedir':savedir})
print "*****test_measurements.py*****"
execfile('test_measurements.py',{'interactive':interactive,'savedir':savedir})
print "*****simple_fit_example.py*****"
execfile('simple_fit_example.py',{'interactive':interactive,'savedir':savedir})
print "*****simple_fit_interactive.py*****"
execfile('simple_fit_interactive.py',{'interactive':interactive,'savedir':savedir})
print "*****alberto_example.py*****"
execfile('alberto_example.py',{'interactive':interactive,'savedir':savedir})

print "*****test_formaldehyde_radex.py*****"
execfile('test_formaldehyde_radex.py',{'interactive':interactive,'savedir':savedir})
print "*****test_formaldehyde.py*****"
execfile('test_formaldehyde.py',{'interactive':interactive,'savedir':savedir})

print "*****vega_echelle.py*****"
execfile('vega_echelle_example.py',{'interactive':interactive,'savedir':savedir})

print "Success!  Or at least, no exceptions..."

print "Running comparison"
execfile('compare_images.py')

