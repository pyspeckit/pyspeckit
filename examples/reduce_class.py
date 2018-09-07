from __future__ import print_function
from spectrum.readers.read_class import class_to_obsblocks
from numpy import *
import glob
import os
try:
    from agpy import print_timing
except:
    print_timing = lambda x: x

@print_timing
def fit_source(sp,debug=False,autorefresh=False,refit=False):

    sp.plotter.autorefresh=autorefresh
    sp.plotter(figure=1,title=sp.header.get('OBJECT'))
    sp.baseline(order=2)
    if debug: import pdb; pdb.set_trace()
    sp.xarr.frequency_to_velocity(center_frequency=sp.header.get('RESTFREQ'),velocity_unit='km/s')
    sp.plotter(xmin=-100,xmax=150,figure=1)
    sp.specfit(negamp=False,limitedmin=[True,True,False,True])
    sp.specfit.seterrspec(usestd=True)
    print(sp.header.get('OBJECT'),sp.specname)
    print(sp.specfit.guesses,sp.specfit.modelpars,sp.specfit.modelerrs,sp.specfit.errspec.mean())
    sp.specfit(negamp=False,limitedmin=[True,True,False,True])
    print(sp.specfit.guesses,sp.specfit.modelpars,sp.specfit.modelerrs,sp.specfit.errspec.mean())
    if debug: pdb.set_trace()

    if (sp.specfit.modelpars[0] <= 2*sp.specfit.modelerrs[0] 
        or sp.specfit.modelpars[2] > 10
        or (sp.specfit.modelpars[2] > 5 and 
            (sp.specfit.modelpars[1] < -90 or sp.specfit.modelpars[1] > 140) ) ):
        print("Fit was declared bad! Model pars:",sp.specfit.modelpars)
        sp.specfit.modelpars = [0 for P in sp.specfit.modelpars]
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
                print("Fit was declared bad after a refit!  Model Pars:",sp.specfit.modelpars)
                sp.specfit.modelpars = [0 for P in sp.specfit.modelpars]
                sp.specfit.clear()
    if sp.specfit.modelpars[1] > 90:
        sp.plotter(ymax=sp.plotter.ymax*1.3,clear=False)
    print(sp.specfit.guesses,sp.specfit.modelpars,sp.specfit.modelerrs,sp.specfit.errspec.mean())
    if debug: pdb.set_trace()

    sp.plotter.refresh()
    savename = "%s_%s_gaussfit.png" % (sp.specname,sp.header['LINE'].strip().replace("+","p"))
    print("Saving fit in %s" % savename)
    sp.plotter.savefig(savename)
    return sp

if __name__ == "__main__":
    import re
    import sys
    import optparse 

    parser=optparse.OptionParser()
    parser.add_option("--hcop","--hco","-H",help="Do HCOP",default=None,action="store_true")
    parser.add_option("--n2hp","-N",help="Do N2HP",default=None,action="store_true")
    parser.add_option("--both","-b",help="Do both HCOP and N2HP",default=None,action="store_true")
    parser.add_option("--refit","-r",help="Refit?",default=True)
    parser.add_option("--interactive","-i",help="Interactive?",default=False,action="store_true")
    options,args = parser.parse_args()
    if len(args) == 1:
        GLOBSTRING = args[0]
        if "class" not in GLOBSTRING:
            import pdb; pdb.set_trace()
        filelist = glob.glob(GLOBSTRING)
    else:
        filelist = args

    if options.both:
        options.n2hp = True
        options.hcop = True

    R = re.compile("class_([0-9]*).smt")
    refit = options.refit
    if refit is False:
        print("Refit disabled.")
        #for filename in filelist:
        #    sessionnum = R.search(filename).groups()[0]
        #    fout_hcop = open('HHT_2011_hcop_bestfits_%s.txt' % sessionnum,'a')
        #    fout_n2hp = open('HHT_2011_n2hp_bestfits_%s.txt' % sessionnum,'a')
        #    n2hp = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HU','SMT-F1M-VU'],line=['N2HP(3-2)','N2H+(3-2)','N2HP','N2HO','N2DP','N2HO(3-2)','N2DP(3-2)'])
        #    hcop = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HL','SMT-F1M-VL'],line=['HCOP(3-2)','HCO+(3-2)'])
        #    for sp in hcop+n2hp:
        #        sp = fit_source(sp.average(),refit=False)
        #        if sp is not None:
        #            if any(X in sp.header.get('LINE') for X in ('N2HP','N2HO','N2DP')):
        #                print >>fout_n2hp,"".join(["%20s" % s 
        #                    for s in [spn.header.get('OBJECT'),spn.header.get('OBSNUM')]+spn.specfit.modelpars+spn.specfit.modelerrs+[spn.specfit.residuals.std()]])
        #            elif 'HCO' in sp.header.get('LINE'):
        #                print >>fout_hcop,"".join(["%20s" % s 
        #                    for s in [spn.header.get('OBJECT'),spn.header.get('OBSNUM')]+spn.specfit.modelpars+spn.specfit.modelerrs+[spn.specfit.residuals.std()]])
        #            else:
        #                print "ERROR: line is ",sp.header.get('LINE')
        #                import pdb; pdb.set_trace()
        #    fout_hcop.close()
        #    fout_n2hp.close()
        #    del hcop 
        #    del n2hp

    else:
        for filename in filelist:
            sessionnum = R.search(filename).groups()[0]
            if options.interactive:
                hcop = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HL','SMT-F1M-VL'],line=['HCOP(3-2)','HCO+(3-2)','HCOP','HCO+'])
                n2hp = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HU','SMT-F1M-VU'],line=['N2HP(3-2)','N2H+(3-2)','N2HP','N2HO','N2DP','N2HO(3-2)','N2DP(3-2)'])
            else:
                if options.hcop:
                    fout_hcop = open('HHT_2011_hcop_bestfits_%s.txt' % sessionnum,'w')
                    print("".join(["%20s" % s for s in ("Source_Name","scannum","amplitude","center","width","amp_err","cen_err","wid_err","RMS")]), file=fout_hcop)
                    hcop = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HL','SMT-F1M-VL'],line=['HCOP(3-2)','HCO+(3-2)','HCOP','HCO+'])
                    for sp in hcop:
                        spn = fit_source(sp.average(),refit=True)
                        if 'HCO' in spn.header.get('LINE'):
                            print("".join(["%20s" % s 
                                for s in [spn.header.get('OBJECT'),spn.header.get('OBSNUM')]+spn.specfit.modelpars+spn.specfit.modelerrs+[spn.specfit.residuals.std()]]), file=fout_hcop)
                        else:
                            print("ERROR: line is ",spn.header.get('LINE'))
                            import pdb; pdb.set_trace()
                    fout_hcop.close()
                    del hcop 

                if options.n2hp:
                    fout_n2hp = open('HHT_2011_n2hp_bestfits_%s.txt' % sessionnum,'w')
                    print("".join(["%20s" % s for s in ("Source_Name","scannum","amplitude","center","width","amp_err","cen_err","wid_err","RMS")]), file=fout_n2hp)
                    n2hp = print_timing(class_to_obsblocks)(filename,telescope=['SMT-F1M-HU','SMT-F1M-VU'],line=['N2HP(3-2)','N2H+(3-2)','N2HP','N2HO','N2DP','N2HO(3-2)','N2DP(3-2)'])
                    #print "Found %i spectra in hcop" % (len(hcop))
                    for sp in n2hp:
                        spn = fit_source(sp.average(),refit=True)
                        if any(X in sp.header.get('LINE') for X in ('N2HP','N2HO','N2DP')):
                            print("".join(["%20s" % s 
                                for s in [spn.header.get('OBJECT'),spn.header.get('OBSNUM')]+spn.specfit.modelpars+spn.specfit.modelerrs+[spn.specfit.residuals.std()]]), file=fout_n2hp)
                        else:
                            print("ERROR: line is ",spn.header.get('LINE'))
                            import pdb; pdb.set_trace()
                    fout_n2hp.close()
                    del n2hp
