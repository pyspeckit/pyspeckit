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
        <tr><td>
                <br><br>
                <h3>Radio Fitting: N<sub>2</sub>H+ example</h3>
            <code class="prettyprint">
import pyspeckit

# The ammonia fitting wrapper requires a dictionary specifying the transition name
# (one of the four specified below) and the filename.  Alternately, you can have the
# dictionary values be pre-loaded Spectrum instances
filenames = {'oneone':'G032.751-00.071_nh3_11_Tastar.fits',
    'twotwo':'G032.751-00.071_nh3_22_Tastar.fits',
    'threethree':'G032.751-00.071_nh3_33_Tastar.fits',
    'fourfour':'G032.751-00.071_nh3_44_Tastar.fits'}

# Fit the ammonia spectrum with some reasonable initial guesses.  It is
# important to crop out extraneous junk and to smooth the data to make the
# fit proceed at a reasonable pace.  
spdict1,spectra1 = pyspeckit.wrappers.fitnh3.fitnh3tkin(filenames,crop=[0,80],tkin=18.65,tex=4.49,column=15.5,fortho=0.9,verbose=False,smooth=6)
        </code>
    </td></tr>
    <tr><td>
    <img align=center width=800 src=images/n2hp_ophA_fit.png title="Fit to the 15 hyperfine components of N2H+ simultaneously.  The (0)'s indicate that this is the 0'th velocity component being fit">
    </td></tr>
    </table>
    <?php include 'navbar.php';?>
    </center>

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
    
    
  </body>
</html>
