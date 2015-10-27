import pyspeckit
# neet matplotlib so we can make mouse-click events from the script
import matplotlib 
import os
# list of annotations so we can clear them
annotations = []
excesslines = []

# get the data from http://cdsarc.u-strasbg.fr/ftp/cats/II/179/sp/hr2421.fit
import urllib2
url = urllib2.urlopen('http://cdsarc.u-strasbg.fr/ftp/cats/II/179/sp/hr2421.fit')
with open('hr2421.fit','wb') as outfile:
    outfile.write(url.read())

# Load the spectrum
sp = pyspeckit.Spectrum('hr2421.fit')

# Plot a particular spectral line
sp.plotter(xmin=4700,xmax=5000)

# Need to fit the continuum first
sp.baseline(interactive=True, subtract=False)

# Left-click to select the fitting region
event1 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,157,316,button=1)
annotations.append( sp.plotter.axis.annotate("First click\n(button 1)", xy=(event1.xdata,event1.ydata),  xytext=(event1.xdata+20,event1.ydata+5e-11),
        textcoords='data', xycoords='data', ha='center',
        va='bottom', arrowprops=dict(arrowstyle="->",
            connectionstyle='arc,rad=0.5', color='green')) )
annotations.append (sp.plotter.axis.annotate("", xy=(4714,9.54e-10),  xytext=(event1.xdata+20,event1.ydata+8e-11),
        textcoords='data', xycoords='data', ha='center', va='bottom',
        arrowprops=dict(arrowstyle="->", connectionstyle='arc,rad=0.5',
            color='green')) )
annotations.append( sp.plotter.axis.annotate("(only the horizontal\n position matters)",xy=(event1.xdata+10,event1.ydata),xycoords='data') )
excesslines.append( sp.plotter.axis.vlines(4714,sp.plotter.ymin.value,sp.plotter.ymax.value,linestyle='--',color='green') )
sp.plotter.refresh()
sp.plotter.savefig('interactive_example_hr2421_baseline_firstclick.png', bbox_inches=None)
event2 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,850,280,button=1)
# event2.xdata and .ydata is None here
event2.xdata = 850
event2.ydata = 850
annotations.append( sp.plotter.axis.annotate("Second click\n(button 1)", xy=(event2.xdata,event2.ydata),  xytext=(event2.xdata-20,event2.ydata+2.5e-11),
        textcoords='data', xycoords='data', ha='center',
        va='bottom', arrowprops=dict(arrowstyle="->",
            connectionstyle='arc,rad=0.5', color='green')) )
annotations.append( sp.plotter.axis.annotate("", xy=(event2.xdata,8.0e-10),  xytext=(event2.xdata-20,event2.ydata+6e-11),
        textcoords='data', xycoords='data', ha='center', va='bottom',
        arrowprops=dict(arrowstyle="->", connectionstyle='arc,rad=0.5',
            color='green')) )
excesslines.append( sp.plotter.axis.vlines(event2.xdata,sp.plotter.ymin.value,sp.plotter.ymax.value,linestyle='--',color='green'))
sp.plotter.refresh()
sp.plotter.savefig('interactive_example_hr2421_baseline_secondclick.png', bbox_inches=None)

sp.baseline.event_manager(event1)
sp.baseline.event_manager(event2)
sp.plotter.savefig('interactive_example_hr2421_baseline_secondclick_highlight.png', bbox_inches=None)

event3 = matplotlib.backend_bases.KeyEvent('button_press_event', sp.plotter.axis.figure.canvas,x=425,y=316,key='x')
annotations.append( sp.plotter.axis.annotate("Third click\n(key 'x')", xy=(event3.xdata,event3.ydata),  xytext=(event3.xdata-20,event3.ydata+5e-11),
        textcoords='data', xycoords='data', ha='center',
        va='bottom', arrowprops=dict(arrowstyle="->",
            connectionstyle='arc,rad=0.5', color='red')) )
annotations.append (sp.plotter.axis.annotate("", xy=(4825,8.67e-10),  xytext=(event3.xdata-20,event3.ydata+8e-11),
        textcoords='data', xycoords='data', ha='center', va='bottom',
        arrowprops=dict(arrowstyle="->", connectionstyle='arc,rad=0.5',
            color='red')) )
annotations.append( sp.plotter.axis.annotate("(to exclude the line\nfrom the baseline fit)",xy=(event3.xdata-55,event3.ydata-5e-11),xycoords='data',color='red') )
excesslines.append( sp.plotter.axis.vlines(4825,sp.plotter.ymin.value,sp.plotter.ymax.value,linestyle='--',color='red') )
sp.plotter.refresh()
sp.plotter.savefig('interactive_example_hr2421_baseline_thirdclick.png', bbox_inches=None)

event4 = matplotlib.backend_bases.KeyEvent('button_press_event', sp.plotter.axis.figure.canvas,x=645,y=280,key='x')
event4.xdata = 900
event4.ydata = 900
annotations.append( sp.plotter.axis.annotate("Fourth click\n(key 'x')", xy=(event4.xdata,event4.ydata),  xytext=(event4.xdata+20,event4.ydata+2.5e-11),
        textcoords='data', xycoords='data', ha='center',
        va='bottom', arrowprops=dict(arrowstyle="->",
            connectionstyle='arc,rad=0.5', color='red')) )
annotations.append( sp.plotter.axis.annotate("", xy=(4905,8.0e-10),  xytext=(event4.xdata+20,event4.ydata+6e-11),
        textcoords='data', xycoords='data', ha='center', va='bottom',
        arrowprops=dict(arrowstyle="->", connectionstyle='arc,rad=0.5',
            color='red')) )
excesslines.append( sp.plotter.axis.vlines(4905,sp.plotter.ymin.value,sp.plotter.ymax.value,linestyle='--',color='red'))
sp.plotter.refresh()
sp.plotter.savefig('interactive_example_hr2421_baseline_fourthclick.png', bbox_inches=None)

sp.baseline.event_manager(event3)
sp.baseline.event_manager(event4)
sp.plotter.savefig('interactive_example_hr2421_baseline_fourthclick_highlight.png', bbox_inches=None)

annotations.append( sp.plotter.axis.annotate("Fifth click - perform the fit\n(button 3;\nbutton 2 will subtract)", 
    xy=(4865, 1e-9),  xytext=(4865, 1e-9), textcoords='data',
    xycoords='data', ha='center', va='bottom') )
event5 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,787,223,button=3)
sp.baseline.event_manager(event5)
sp.plotter.savefig('interactive_example_hr2421_baseline_fifthclick_fit.png', bbox_inches=None)

for LC in excesslines:
    if LC in sp.plotter.axis.collections:
        sp.plotter.axis.collections.remove(LC)
for AN in annotations:
    if AN in sp.plotter.axis.texts:
        sp.plotter.axis.texts.remove(AN)
sp.plotter.refresh()

# ********************************************
# Start up an interactive line-fitting session
# ********************************************
sp.specfit(interactive=True)

# Left-click to select the fitting region
event1 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,257,316,button=1)
annotations.append( sp.plotter.axis.annotate("First click\n(button 1)", xy=(event1.xdata,event1.ydata), xytext=(event1.xdata+20,event1.ydata),
        textcoords='data', xycoords='data', ha='center',
        va='bottom', arrowprops=dict(arrowstyle="->",
            connectionstyle='arc,rad=0.5', color='green')) )
annotations.append (sp.plotter.axis.annotate("", xy=(event1.xdata,9.26e-10),  xytext=(event1.xdata+20,event1.ydata+5e-11),
        textcoords='data', xycoords='data', ha='center', va='bottom',
        arrowprops=dict(arrowstyle="->", connectionstyle='arc,rad=0.5',
            color='green')) )
annotations.append( sp.plotter.axis.annotate("(only the horizontal\n position matters)",xy=(event1.xdata+2,event1.ydata-1e-10),xycoords='data') )
excesslines.append( sp.plotter.axis.vlines(event1.xdata,sp.plotter.ymin,sp.plotter.ymax,linestyle='--',color='green') )
sp.plotter.refresh()
sp.plotter.savefig('figures/interactive_example_hr2421_firstclick.png', bbox_inches=None)

event2 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,732,280,button=1)
annotations.append( sp.plotter.axis.annotate("Second click\n(button 1)", xy=(event2.xdata,event2.ydata),  xytext=(event2.xdata+20,event2.ydata+3e-11),
        textcoords='data', xycoords='data', ha='center',
        va='bottom', arrowprops=dict(arrowstyle="->",
            connectionstyle='arc,rad=0.5', color='green')) )
annotations.append( sp.plotter.axis.annotate("", xy=(4943,8.29e-10),  xytext=(event2.xdata+20,event2.ydata+6e-11),
        textcoords='data', xycoords='data', ha='center', va='bottom',
        arrowprops=dict(arrowstyle="->", connectionstyle='arc,rad=0.5',
            color='green')) )
excesslines.append( sp.plotter.axis.vlines(4943,sp.plotter.ymin,sp.plotter.ymax,linestyle='--',color='green'))
sp.plotter.refresh()
sp.plotter.savefig('figures/interactive_example_hr2421_secondclick.png', bbox_inches=None)

sp.specfit.event_manager(event1)
sp.specfit.event_manager(event2)
sp.plotter.savefig('figures/interactive_example_hr2421_secondclick_highlight.png', bbox_inches=None)

for LC in excesslines:
    if LC in sp.plotter.axis.collections:
        sp.plotter.axis.collections.remove(LC)
for AN in annotations:
    if AN in sp.plotter.axis.texts:
        sp.plotter.axis.texts.remove(AN)
sp.plotter.refresh()

event3 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,523,194,button=2)
annotations.append( sp.plotter.axis.annotate("Third click\n(button 2)", xy=(event3.xdata,event3.ydata),  xytext=(event3.xdata+40,event3.ydata-5e-11),
        textcoords='data', xycoords='data', ha='center',
        va='bottom', arrowprops=dict(arrowstyle="->",
            connectionstyle='arc,rad=0.5', color='orange', shrinkB=5)) )
sp.specfit.event_manager(event3)
sp.plotter.savefig('figures/interactive_example_hr2421_thirdclick.png', bbox_inches=None)

event4 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,485,264,button=2)
annotations.append( sp.plotter.axis.annotate("Fourth click\n(button 2)", xy=(event4.xdata,event4.ydata),  xytext=(event4.xdata-20,event4.ydata-5e-11),
        textcoords='data', xycoords='data', ha='center',
        va='bottom', arrowprops=dict(arrowstyle="->",
            connectionstyle='arc,rad=0.5', color='orange', shrinkB=5)) )
sp.specfit.event_manager(event4)
sp.plotter.savefig('figures/interactive_example_hr2421_fourthclick.png', bbox_inches=None)

model = sp.specfit.Registry.multifitters['gaussian'].n_modelfunc(sp.specfit.guesses)(sp.xarr) + sp.baseline.basespec
sp.plotter.axis.plot(sp.xarr,model,color='b')
annotations.append( sp.plotter.axis.annotate("The guessed model", 
    xy=(4770,8.2e-10),  xytext=(4770,8.2e-10), textcoords='data',
    xycoords='data', ha='center', va='bottom', color='blue') )
sp.plotter.refresh()
sp.plotter.savefig('figures/interactive_example_hr2421_gaussmodelguess.png', bbox_inches=None)

for LC in excesslines:
    if LC in sp.plotter.axis.collections:
        sp.plotter.axis.collections.remove(LC)
for AN in annotations:
    if AN in sp.plotter.axis.texts:
        sp.plotter.axis.texts.remove(AN)
annotations.append( sp.plotter.axis.annotate("Fifth click - perform the fit\n(button 3)", 
    xy=(4855, 9.1e-10),  xytext=(4855,9.1e-10), textcoords='data',
    xycoords='data', ha='center', va='bottom') )
event5 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,787,223,button=3)
sp.specfit.event_manager(event5)
sp.plotter.savefig('figures/interactive_example_hr2421_fifthclick_fit.png', bbox_inches=None)

#sp.plotter.figure.savefig(savedir+'hr2421_interactive_selectregion.png', bbox_inches=None)
#sp.plotter.figure.savefig(savedir+'hr2421_interactive_guesses.png', bbox_inches=None)
#event5 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,611,247,button=3)
#sp.specfit.event_manager(event5)
