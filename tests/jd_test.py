# grab the GBT reader from pyspeckit
from pyspeckit.spectrum.readers import gbt

# declare our session name and read the file
A029 = gbt.GBTSession('AGBT11B_029_01.raw.acs.fits')

# let's see what's in it
print A029

# reduce one of the targets
ob1 = A029.reduce_target('003918.9+402158.4')
# you can also access ob1 as A029.target['003918.9+402158.4'] or simply A029['003918.9+402158.4']

# now see what it contains
print ob1.spectra.keys()

# Average polarizations & feeds for each IF
# (this may be poorly named)
ob1.average_IFs()
# again, check on the contents
print ob1.spectra.keys()

# now plot IF0
ob1['if0'].plotter()

# convert to velocity units
ob1['if0'].xarr.convert_to_unit('km/s')
# and replot
ob1['if0'].plotter()

# crop the spectrum and replot
ob1['if0'].crop(-600,-525,units='km/s')
ob1['if0'].plotter(reset_xlimits=True)
ob1['if0'].specfit()

# Take advantage of python!
import pylab
pylab.figure() # create a new figure
pylab.imshow(ob1.blocks['A9ON1'].data)

# Start getting a bit tricky...
# Reduce ALL spectra (this can take time)
A029.reduce_all()
print A029

# plot all the reduced IF0's in different windows
for target_name in A029.targets:
    A029['target_name']['if0'].plotter()

