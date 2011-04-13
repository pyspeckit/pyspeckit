import spectrum
from spectrum.readers.read_class import class_to_obsblocks
from numpy import *
import glob
import os
try:
    from agpy import print_timing
except:
    print_timing = lambda(x): x

@print_timing
def fit_source(sp,debug=False,autorefresh=False,refit=False):

    sp.plotter.autorefresh=autorefresh
    sp.plotter(figure=1,title=sp.header.get('OBJECT'))
    sp.baseline(order=2)
    if debug: import pdb; pdb.set_trace()
    sp.xarr.frequency_to_velocity(center_frequency=sp.header.get('RESTFREQ'),velocity_units='km/s')
    sp.plotter(xmin=-100,xmax=150,figure=1)
    sp.specfit(negamp=False,limitedmin=[True,True,False,True])
    sp.specfit.seterrspec(usestd=True)
    print sp.header.get('OBJECT'),sp.specname
    print sp.specfit.guesses,sp.specfit.modelpars,sp.specfit.modelerrs,sp.specfit.errspec.mean()
    sp.specfit(negamp=False,limitedmin=[True,True,False,True])
    print sp.specfit.guesses,sp.specfit.modelpars,sp.specfit.modelerrs,sp.specfit.errspec.mean()
    if debug: pdb.set_trace()
    #raw_input('Wait')
    if (sp.specfit.modelpars[0] <= 2*sp.specfit.modelerrs[0] 
        or sp.specfit.modelpars[2] > 10
        or (sp.specfit.modelpars[2] > 5 and 
            (sp.specfit.modelpars[1] < -90 or sp.specfit.modelpars[1] > 140) ) ):
        print "Fit was declared bad!"
        sp.specfit.clear()
        sp.plotter(reset_ylimits=True)
    else:
        sp.baseline(excludefit=True,order=2)
        sp.specfit(negamp=False,limitedmin=[True,True,False,True])
        sp.plotter(reset_ylimits=True,clear=False)
        if (sp.specfit.modelpars[0] == 0 or sp.specfit.modelpars[2] == 0
            or sp.specfit.modelpars[2] > 10
            or (sp.specfit.modelpars[2] > 5 and 
                (sp.specfit.modelpars[1] < -90 or sp.specfit.modelpars[1] > 140) ) ):
            print "Fit was declared bad after a refit!"
            sp.specfit.clear()
    if sp.specfit.modelpars[1] > 90:
        sp.plotter(ymax=sp.plotter.ymax*1.3,clear=False)
    print sp.specfit.guesses,sp.specfit.modelpars,sp.specfit.modelerrs,sp.specfit.errspec.mean()
    if debug: pdb.set_trace()
    #raw_input("Hold up")
    sp.plotter.refresh()
    savename = "%s_%s_gaussfit.png" % (sp.specname,sp.header['LINE'].strip().replace("+","p"))
    print "Saving fit in %s" % savename
    sp.plotter.savefig(savename)
    return sp

if __name__ == "__main__":
    refit = True
    if refit is False:
        fout_hcop = open('HHT_2011_hcop_bestfits.txt','a')
        fout_n2hp = open('HHT_2011_n2hp_bestfits.txt','a')
        for filename in glob.glob("class*smt"):
            n2hp = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HU','SMT-F1M-VU'],line=['N2HP(3-2)','N2H+(3-2)'])
            hcop = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HL','SMT-F1M-VL'],line=['HCOP(3-2)','HCO+(3-2)'])
            for sp in hcop+n2hp:
                sp = fit_source(sp.average(),refit=False)
                if sp is not None:
                    if 'N2HP' in sp.header.get('LINE'):
                        print >>fout_n2hp,"".join(["%20s" % s 
                            for s in [sp.header.get('OBJECT')]+sp.specfit.modelpars+sp.specfit.modelerrs])
                    elif 'HCO' in sp.header.get('LINE'):
                        print >>fout_hcop,"".join(["%20s" % s 
                            for s in [sp.header.get('OBJECT')]+sp.specfit.modelpars+sp.specfit.modelerrs])
                    else:
                        print "ERROR: line is ",sp.header.get('LINE')

    else:
        fout_hcop = open('HHT_2011_hcop_bestfits.txt','w')
        fout_n2hp = open('HHT_2011_n2hp_bestfits.txt','w')
        print >>fout_hcop,"".join(["%20s" % s for s in ("Source_Name","amplitude","center","width","amp_err","cen_err","wid_err")])
        print >>fout_n2hp,"".join(["%20s" % s for s in ("Source_Name","amplitude","center","width","amp_err","cen_err","wid_err")])
        for filename in glob.glob("class*smt"):
            #n2hp = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HU','SMT-F1M-VU'],line=['N2HP(3-2)','N2H+(3-2)'])
            hcop = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HL','SMT-F1M-VL'],line=['HCOP(3-2)','HCO+(3-2)'])
            print "Found %i spectra in hcop" % (len(hcop))
            for sp in hcop: #+n2hp:
                spn = fit_source(sp.average(),refit=True)
                if 'N2HP' in spn.header.get('LINE'):
                    print >>fout_n2hp,"".join(["%20s" % s 
                        for s in [spn.header.get('OBJECT')]+spn.specfit.modelpars+spn.specfit.modelerrs])
                elif 'HCO' in spn.header.get('LINE'):
                    print >>fout_hcop,"".join(["%20s" % s 
                        for s in [spn.header.get('OBJECT')]+spn.specfit.modelpars+spn.specfit.modelerrs])
                else:
                    print "ERROR: line is ",spn.header.get('LINE')

    fout_hcop.close()
    fout_n2hp.close()
