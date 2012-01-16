import ImageChops
import numpy

def equal(im1, im2):
    return ImageChops.difference(im1, im2).getbbox() is None

def rmsdiff(im1, im2):
    "Calculate the root-mean-square difference between two images"

    diffim = ImageChops.difference(im1, im2)

    return (numpy.array(diffim)**2).sum()**0.5

if __name__ == "__main__":
    import optparse
    import glob
    from PIL import Image

    parser=optparse.OptionParser()
    parser.add_option('--version',help='Version number to compare to',default='449')
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
    filelist1 = [s.split('/')[1] for s in glob.glob("%s*png" % dir1)]
    filelist2 = [s.split('/')[1] for s in glob.glob("%s*png" % dir2)]

    print "Comparing %i files in %s to %i files in %s" % (len(filelist1), options.version, len(filelist2), currentversion)
    for fn in filelist1:
        if fn in filelist2:
            im1 = Image.open(dir1+fn)
            im2 = Image.open(dir2+fn)
            if not equal(im1,im2):
                print "%s differs from version %s to %s.  RMS: %f" % (fn, options.version, currentversion, rmsdiff(im1,im2))
            else:
                print "%s OK in both %s and %s" % (fn, options.version, currentversion)
        else:
            print "%s exists in %s, but not in %s" % (fn, options.version, currentversion)
