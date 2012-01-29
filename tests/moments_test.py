from pyspeckit.spectrum.moments import moments#,moments2,moments3
from pyspeckit.spectrum.models.inherited_gaussfitter import gaussian
import numpy as np
import pylab

def plot_and_test(args, xarr=np.linspace(-25,25,1000), moments=moments, nsigcut=None):
    arrargs = np.array([0] + list(args))

    gg = gaussian(xarr, *args)

    noise = np.random.randn(1000)

    moments_nonoise = np.array(moments(xarr,gg,veryverbose=True))

    fmtstring = "h: %6.2f a: %6.2f x: %6.2f w: %6.2f "
    print "No noise"
    print "Moments:    "+fmtstring % tuple(moments_nonoise) 
    print "Difference: "+fmtstring % tuple(moments_nonoise-arrargs)
    pylab.subplot(2,2,1)
    pylab.plot(xarr,gg,'k')
    pylab.plot(xarr,gaussian(xarr,*moments_nonoise[1:]),'r')
    print "diff^2: %0.2f" % (((gg-gaussian(xarr,*moments_nonoise[1:]))**2).sum())

    def plottest_withnoise(noiseamp, spnum, xarr=xarr, gg=gg, noise=noise):
        fitspec = gg+noiseamp*noise
        moments_noise1 = np.array(moments(xarr,fitspec,veryverbose=True,nsigcut=nsigcut))
        print "Noise = %i%% of peak" % (round(noiseamp/args[0]*100))
        print "Moments:    "+fmtstring % tuple(moments_noise1) 
        print "Difference: "+fmtstring % tuple(moments_noise1-arrargs)
        pylab.subplot(2,2,spnum)
        pylab.plot(xarr,fitspec,'k')
        pylab.plot(xarr,gaussian(xarr,*moments_noise1[1:]),'r')
        print "chi^2: %0.2f" % ( (((fitspec-gaussian(xarr,*moments_noise1[1:]))**2) / ((noise*noiseamp)**2)).sum() / len(noise) ),
        print "correct chi^2: %0.2f" % ( (((fitspec-gg)**2) / ((noise*noiseamp)**2)).sum() / len(noise) )

    for noiseamp,spnum in zip((0.5,2,5),(2,3,4)):
        plottest_withnoise(noiseamp, spnum)



if __name__=="__main__":

    pylab.figure(1); pylab.clf()
    args = (5.1, 8.1, 0.2)
    plot_and_test(args)
    pylab.figure(2); pylab.clf()
    plot_and_test(args,nsigcut=1)
    pylab.figure(3); pylab.clf()
    plot_and_test(args,nsigcut=2)

    """
    print
    pylab.figure(2)
    pylab.clf()
    args = (11.5, 3.1, 0.7)
    plot_and_test(args)

    print
    pylab.figure(3)
    pylab.clf()
    args = (5.1, 8.1, 0.2)
    plot_and_test(args, moments=moments2)

    print
    pylab.figure(4)
    pylab.clf()
    args = (11.5, 3.1, 0.7)
    plot_and_test(args, moments=moments2)

    print
    pylab.figure(5)
    pylab.clf()
    args = (5.1, 8.1, 5.2)
    plot_and_test(args)

    print
    pylab.figure(6)
    pylab.clf()
    args = (5.1, 8.1, 0.2)
    plot_and_test(args, moments=moments3)

    print
    pylab.figure(7)
    pylab.clf()
    args = (11.5, 3.1, 0.7)
    plot_and_test(args, moments=moments3)

    print
    pylab.figure(8)
    pylab.clf()
    args = (5.1, 8.1, 5.2)
    plot_and_test(args, moments=moments3)
    """


    pylab.draw()
    pylab.ion()
    pylab.show()
