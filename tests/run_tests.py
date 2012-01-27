interactive=False
import subprocess
versnum = subprocess.Popen(["hg","id","--num"],stdout=subprocess.PIPE).communicate()[0].strip().strip("+")
savedir = "tests_%s/" % versnum
import os
if not os.path.exists(savedir):
    os.mkdir(savedir)

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

print "Success!  Or at least, no exceptions..."

print "Running comparison"
execfile('compare_images.py')
