
# Import the ALFALFA reader (deals with the specific structure type assumed for
# an ALFALFA src file)
from pyspeckit.spectrum.readers import alfalfa

# Read the full contents of an ALFALFA .src file
HI1459 = alfalfa.read_alfalfa_file('HI145947.9+152515_1500+15d.src')
# print out the contents of the dictionary
print HI1459

# make a shortcut for the source
src = HI1459['145947.9+152515']

# we can plot the whole object, but it has 6 spectra, so we're going to get 6 overlapping spectra
src.plotter()

# it may be more useful to plot each item in the "block"
src.ploteach()

# or you can just plot an individual one (this will not open a new window; it will overlay on a previous plot)
src[0].plotter(color='r',clear=False)



