
<html>
    <head>
        <link rel="stylesheet" type="text/css" href="style.css">
        <link href="prettify.css" type="text/css" rel="stylesheet">
        <link rel="icon" type="image/gif" href="images/logo.ico">
        <script type="text/javascript" src="prettify.js"></script>
        <title> PySpecKit - Examples </title>
        <h2> Examples </h2>
    </head>
    <br>
    <br>
<body onload="prettyPrint()"> 
    <center>
    <table width=800>
        <tr><td><br> <b> <a href="sdss_example.html"> Optical Fitting - AGN from SDSS </a></b></td></tr>
        <tr><td><br> <b> <a href="nh3_example.html"> Radio Fitting - NH<sub>3</sub> example </a></b></td></tr>
        <tr><td><br> <b> <a href="h2co_example.html"> Radio Fitting - H<sub>2</sub>CO RADEX example </a></b></td></tr>
        <tr><td><br> <b> <a href="doublet_example.html"> Doublet Fitting - SII </a></b></td></tr>
        <tr><td>
                <br><br>
                <h3>Simple Radio Fitting: HCO+ example</h3>
            <code class="prettyprint">
              <!--<? readfile("../tests/simple_fit_example.py"); ?>-->
import pyspeckit

# load a FITS-compliant spectrum
spec = pyspeckit.Spectrum('10074-190_HCOp.fits')
# The units are originally frequency (check this by printing spec.xarr.units).
# I want to know the velocity.  Convert!
# Note that this only works because the reference frequency is set in the header
spec.xarr.frequency_to_velocity()
# Default conversion is to m/s, but we traditionally work in km/s
spec.xarr.convert_to_unit('km/s')
# plot it up!
spec.plotter()
# Subtract a baseline (the data is only 'mostly' reduced)
spec.baseline()
# Fit a gaussian.  We know it will be an emission line, so we force a positive guess
spec.specfit(negamp=False)
# Note that the errors on the fits are larger than the fitted parameters.
# That's because this spectrum did not have an error assigned to it.  
# Let's use the residuals:
spec.specfit.plotresiduals()
# Now, refit with error determined from the residuals:
# (we pass in guesses to save time / make sure nothing changes)
spec.specfit(guesses=spec.specfit.modelpars)

# Save the figures to put on the web....
spec.plotter.figure.savefig("simple_fit_example_HCOp.png")
spec.specfit.residualaxis.figure.savefig("simple_fit_example_HCOp_residuals.png")

# Also, let's crop out stuff we don't want...
spec.crop(-100,100)
# replot after cropping (crop doesn't auto-refresh)
spec.plotter()
# replot the fit without re-fitting
spec.specfit.plot_fit()
# show the annotations again
spec.specfit.annotate()
spec.plotter.figure.savefig("simple_fit_example_HCOp_cropped.png")
            </code>
      <tr><td><center><img align=center width=800 src="images/simple_fit_example_HCOp.png" title="Sample HCO+ spectrum fitted with a gaussian"></center></tr></td>
      <tr><td><center><img align=center width=800 src="images/simple_fit_example_HCOp_residuals.png" title="Residuals of the gaussian fit from the previous figure"></center></tr></td>
      <tr><td><center><img align=center width=800 src="images/simple_fit_example_HCOp_cropped.png" title="A zoomed-in, cropped version of the spectrum.  With the 'crop' command, the excess data is discarded."></center></tr></td>
    
    </table>
    <?php include 'navbar.php';?>
    </center>
    
  </body>
<script type="text/javascript">

  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-6253200-7']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();

</script>
  
</html>
