# grab the GBT reader from pyspeckit
from pyspeckit.spectrum.readers import gbt

# declare our session name and read the file
A029 = gbt.GBTSession('AGBT11B_029_01.raw.acs.fits')

# let's see what's in it
print A029

# reduce one of the targets
ob1 = A029.reduce_target('003918.9+402158.4')

# now see what it contains
print ob1.spectra

# Average polarizations & feeds for each IF
# (this may be poorly named)
ob1.average_IFs()
# again, check on the contents
print ob1.spectra

# now plot IF0
ob1['if0'].plotter()

# convert to velocity units
ob1['if0'].xarr.convert_to_unit('km/s')
# and replot
ob1['if0'].plotter()

# crop the spectrum and replot
ob1['if0'].crop(-600,-525,units='km/s')
ob1['if0'].plotter(figure=figure(),reset_xlimits=True)

