import ImageChops
import numpy
import pylab

def equal(im1, im2):
    return ImageChops.difference(im1, im2).getbbox() is None

def rmsdiff(im1, im2):
    "Calculate the root-mean-square difference between two images"

    diffim = ImageChops.difference(im1, im2)

    return (numpy.array(diffim)**2).sum()**0.5

def plot_difference_image(im1, im2, savename=None):
    """
    In 3 subplots, show the input, output, and difference images

    For histogram normalization reasons, the values 0 and 255 are excluded.
    The 'alpha' map is also excluded.  Ideally, we'd use log=True in these, but
    there's a bug in matplotlib
    (https://github.com/matplotlib/matplotlib/issues/196) that makes this too
    ugly to be useful
    """
    im1arr = numpy.array(im1) 
    im2arr = numpy.array(im2) 
    diffim = ImageChops.difference(im1,im2)
    diffarr = numpy.array(diffim)
    xmin,ymin,xmax,ymax = diffim.getbbox()
    ignore_diff = ( (diffarr == 0).sum(axis=2) == 4 )[ymin:ymax,xmin:xmax]

    pylab.figure(1)
    pylab.clf()
    pylab.subplot(2,3,1)
    pylab.imshow(im1arr[ymin:ymax,xmin:xmax,:])
    pylab.gca().invert_yaxis()
    pylab.subplot(2,3,2)
    pylab.imshow(im2arr[ymin:ymax,xmin:xmax,:])
    pylab.gca().invert_yaxis()
    pylab.subplot(2,3,3)
    diffarr[ymin:ymax,xmin:xmax,3] = 255*(diffarr[ymin:ymax,xmin:xmax,:3].sum(axis=2) > 0)
    pylab.imshow(diffarr[ymin:ymax,xmin:xmax,:])
    pylab.gca().invert_yaxis()
    
    pylab.subplot(2,3,4)
    ignore_im1 = ( (im1arr == 255).sum(axis=2) == 4 )[ymin:ymax,xmin:xmax]
    #pylab.hist(im1arr[ymin:ymax,xmin:xmax,3][True-ignore_im1],edgecolor='black',color='black',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist(im1arr[ymin:ymax,xmin:xmax,0][True-ignore_im1],edgecolor='red',color='red',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist(im1arr[ymin:ymax,xmin:xmax,1][True-ignore_im1],edgecolor='green',color='green',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist(im1arr[ymin:ymax,xmin:xmax,2][True-ignore_im1],edgecolor='blue',color='blue',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')

    pylab.subplot(2,3,5)
    ignore_im2 = ( (im2arr == 255).sum(axis=2) == 4 )[ymin:ymax,xmin:xmax]
    #pylab.hist(im2arr[ymin:ymax,xmin:xmax,3][True-ignore_im2],edgecolor='black',color='black',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist(im2arr[ymin:ymax,xmin:xmax,0][True-ignore_im2],edgecolor='red',color='red',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist(im2arr[ymin:ymax,xmin:xmax,1][True-ignore_im2],edgecolor='green',color='green',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist(im2arr[ymin:ymax,xmin:xmax,2][True-ignore_im2],edgecolor='blue',color='blue',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')

    pylab.subplot(2,3,6)
    #pylab.hist((diffarr)[ymin:ymax,xmin:xmax,3][True-ignore_diff],edgecolor='black',color='black',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist((diffarr)[ymin:ymax,xmin:xmax,0][True-ignore_diff],edgecolor='red',color='red',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist((diffarr)[ymin:ymax,xmin:xmax,1][True-ignore_diff],edgecolor='green',color='green',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')
    pylab.hist((diffarr)[ymin:ymax,xmin:xmax,2][True-ignore_diff],edgecolor='blue',color='blue',alpha=0.25,bins=pylab.linspace(1,254,33),histtype='stepfilled')

    if savename is not None:
        pylab.savefig(savename)


if __name__ == "__main__":
    import optparse
    import glob
    from PIL import Image

    parser=optparse.OptionParser()
    parser.add_option('--version',help='Version number to compare to',default='521')
    parser.add_option('--currentversion',help='Current version',default='tip')
    parser.set_usage("%prog [options]")
    parser.set_description(
    """
    Compare the images!
    """)

    options,args = parser.parse_args()

    if options.currentversion == 'tip':
        import subprocess
        currentversion = subprocess.Popen(["hg","id","--num"],stdout=subprocess.PIPE).communicate()[0].strip().strip("+")
    
    dir1 = "tests_%s/" % (options.version)
    dir2 = "tests_%s/" % (currentversion)
    filelist1 = [s.split('/')[1] for s in glob.glob("%s*png" % dir1) if "compare.png" not in s]
    filelist2 = [s.split('/')[1] for s in glob.glob("%s*png" % dir2) if "compare.png" not in s]

    print "Comparing %i files in %s to %i files in %s" % (len(filelist1), options.version, len(filelist2), currentversion)
    if int(options.version) < 521:
        print "WARNING: Before version 521, the figures may have been displayed\
        in an interactive (ion()) window, which made them smaller, resulting in\
        big differences."
    for fn in filelist1:
        if fn in filelist2:
            im1 = Image.open(dir1+fn)
            im2 = Image.open(dir2+fn)
            if not equal(im1,im2):
                RMS = rmsdiff(im1,im2)
                print "%s differs from version %s to %s.  RMS: %f" % (fn, options.version, currentversion, RMS)
                if RMS > 1000:
                    plot_difference_image(im1,im2,dir2+fn.replace(".png","_compare.png"))
            else:
                print "%s OK in both %s and %s" % (fn, options.version, currentversion)
        else:
            print "%s exists in %s, but not in %s" % (fn, options.version, currentversion)
