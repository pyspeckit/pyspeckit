import os
import shutil
import urllib
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

# pyfits depends on numpy.  pywcs hard-depends on numpy (won't install without it)
install_order = ('python','numpy','pyfits','pywcs','matplotlib','boto','Imaging','ordereddict')#,'lmfit')
print_order = ('python','pywcs','numpy','matplotlib','pyfits')

NUMPY_VERSIONS = ['1.4.1', '1.5.1', '1.6.1']
PYTHON_VERSIONS = ['2.6', '2.7']
PYFITS_VERSIONS = ['2.4.0', '3.0.1', '3.0.2', '3.0.3', '3.0.4', '3.0.5']
PYWCS_VERSIONS = ['1.8.1-4.4.4', '1.9-4.4.4', '1.10-4.7', '1.11-4.8.2']
MATPLOTLIB_VERSIONS = ['1.1.0'] # matplotlib 1.0.0 doesn't exist @ pypi, follows nonstandard path @ sourceforge
                                # matplotlib 1.0.1 simply wouldn't install; it died on libpng if I used the 1.1.0 make.osx or died on a bad tarball if I used the given version
BOTO_VERSIONS = ['2.2.2']
PIL_VERSIONS = ['1.1.7']
# http://effbot.org/downloads/Imaging-1.1.7.tar.gz
#LMFIT_VERSIONS = ['0.5']
# https://github.com/newville/lmfit-py/tarball/0.5
ORDEREDDICT_VERSIONS = ['1.1']

pwd = os.getcwd()
os.putenv('PATH','/usr/bin:'+os.getenv('PATH'))

package_dict = OrderedDict({
        'python':PYTHON_VERSIONS,
        'pywcs':PYWCS_VERSIONS,
        'numpy':NUMPY_VERSIONS,
        'matplotlib':MATPLOTLIB_VERSIONS,
        'pyfits':PYFITS_VERSIONS,
        'boto':BOTO_VERSIONS,
        'Imaging':PIL_VERSIONS,
        #'lmfit':LMFIT_VERSIONS,
        'ordereddict':ORDEREDDICT_VERSIONS,
        })

mpl_deps_path = '/Users/adam/repos/jenkins_tests/mpl_dependencies/'
osx_ver = '10.6'
matplotlib_flags = (
    'PKG_CONFIG_PATH={mpl}/lib/pkgconfig'.format(mpl=mpl_deps_path) + 
    'MACOSX_DEPLOYMENT_TARGET=10.6' + 
    'CFLAGS="-arch i386 -arch x86_64 -I${mpl}/include -I${mpl}/include/freetype2 -isysroot /Developer/SDKs/MacOSX{osx}.sdk"'.format(mpl=mpl_deps_path,osx=osx_ver)+
    'LDFLAGS="-arch i386 -arch x86_64 -L${mpl}/lib -syslibroot,/Developer/SDKs/MacOSX${osx}.sdk"'.format(mpl=mpl_deps_path,osx=osx_ver) )

matplotlib_command = "PATH={pybinpath}:$PATH make -j 6 -f make.osx PYVERSION={pyver} PREFIX={pyrootpath} fetch deps mpl_install >> logfile"

package_install_prefixes = {'numpy':'',
        'pyfits':'',
        'pywcs':'',
        'boto':'',
        'matplotlib':'',
        'Imaging':'',
        #'lmfit':'',
        'ordereddict':'',
        }

package_install_flags = {'numpy':' --fcompiler=g95',
        'pyfits':'',
        'pywcs':'',
        'boto':'',
        'matplotlib':'',
        'Imaging':'',
        #'lmfit':'',
        'ordereddict':'',
        }

package_import_names = {'numpy':'numpy',
        'pyfits':'pyfits',
        'pywcs':'pywcs',
        'boto':'boto',
        'matplotlib':'matplotlib',
        'Imaging':'PIL',
        #'lmfit':'lmfit',
        'ordereddict':'ordereddict',
        }

package_url = {'numpy':'http://sourceforge.net/projects/numpy/files/NumPy/{0}/numpy-{0}.tar.gz',
        'pyfits':'http://pypi.python.org/packages/source/p/pyfits/pyfits-{0}.tar.gz',
        'pywcs':'http://stsdas.stsci.edu/astrolib/pywcs-{0}.tar.gz',
        'matplotlib':'http://pypi.python.org/packages/source/m/matplotlib/matplotlib-{0}.tar.gz',
        'boto':'http://boto.googlecode.com/files/boto-{0}.tar.gz',
        'Imaging':'http://effbot.org/downloads/Imaging-{0}.tar.gz',
        #'lmfit':'https://github.com/newville/lmfit-py/tarball/{0}',
        'ordereddict':'http://pypi.python.org/packages/source/o/ordereddict/ordereddict-{0}.tar.gz'
        # 'matplotlib':'http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.0/matplotlib-{0}.tar.gz',
        }

# download and untar the packages
for package_name,versions in package_dict.items():
    for vers in versions:
        package_path = '{0}-{1}'.format(package_name,vers)
        #if os.path.exists(package_path):
        #    shutil.rmtree(package_path)

        tarfile = package_path+'.tar.gz'
        if package_name in package_url:
            tar_url = package_url[package_name].format(vers)
            if not os.path.exists(tarfile):
                urllib.urlretrieve(tar_url,tarfile)
            if not os.path.exists(package_path):
                os.system('tar xzvf %s' % tarfile)
                os.system('chmod -R +w ' + package_path)


# MAJOR HACK - matplotlib 1.0.1 make.osx simply does not work!  So, get rid of the pos...
shutil.copy('matplotlib-1.1.0/make.osx','matplotlib-1.0.1/make.osx')

def cross_product(ss,row=[],level=0):
    """ from http://code.activestate.com/recipes/159975-one-liner-to-generate-the-cross-product-of-an-arbi/"""
    return (len(ss)>1 
            and reduce(lambda x,y:x+y,[cross_product(ss[1:],row+[i],level+1) for i in ss[0]]) 
            or [row+[i] for i in ss[0]])

# make a set of lists of install versions, always sorted by install order
crossprod_versions = cross_product([[(name,value)  for value in package_dict[name]] for name in install_order])

# add Python versions (don't build them)
#package_dict = OrderedDict({'python':PYTHON_VERSIONS}.items()+package_dict.items()) # DO include python vers. in path (but don't try to unpack it)
packages_versions = ([["{0}{1}".format(name,value)  for value in package_dict[name]] for name in print_order if name is not 'boto'])
paths = ["-".join(names) for names in cross_product(packages_versions)]

for versions,targetpath in zip(crossprod_versions,paths):

    for target,vers in versions:
        print target,vers
        if target == 'python':
            versions.remove((target,vers))
            python_version = vers

    if python_version in ['3.1', '3.2']:
        if numpy_version == '1.4.1':
            continue
        if versions == ('numpy','1.4.1'):
            continue
    # this should have been eliminated at the root level...
    if not "python{0}".format(python_version) in targetpath:
        continue

    if 'boto' in targetpath: # this is purely a hack / sanity check: I don't want to accidentally create a bunch of undesired directories
        raise
    if targetpath.count('python') > 1:
        import pdb; pdb.set_trace()

    fullpath = "{0}/{1}/".format(pwd, targetpath)
    print fullpath

    # Skip if already exists
    #if os.path.exists(targetpath):
    #    continue

    # Generate virtual environment
    if not os.path.exists(targetpath):
        os.system('virtualenv-{pv} --no-site-packages --distribute {targetpath}'.format(targetpath=targetpath,pv=python_version))

    python = '{cwd}/{targetpath}/bin/python{pv}'.format(targetpath=targetpath,pv=python_version,cwd=pwd)
    interpreter = 'CC=gcc '+python
    #import pdb; pdb.set_trace()


    for target,vers in versions:
        #if target =='pywcs' and vers == '1.11-4.8.2':
        #    interpreter = interpreter.replace('CC=gcc ', '')
        importname = package_import_names[target]
        tryvalue = os.system(python + " -c 'import {0}'".format(importname))
        if tryvalue != 0:
            print "{0}{2} was not installed (import returned {1}).  Installing now.".format(importname,tryvalue,vers)
            os.chdir('{0}-{1}'.format(target,vers))
            if target == 'matplotlib':
                os.system(matplotlib_command.format(pybinpath=fullpath+'/bin/', pyrootpath=fullpath, pyver=python_version ))
            else:
                os.system(package_install_prefixes[target] +" "+ interpreter + ' setup.py build'+package_install_flags[target]+' >> logfile')
                os.system(package_install_prefixes[target] +" "+ interpreter + ' setup.py install >> logfile')
                os.system('chmod -R +w *')
            os.chdir('../')
        tryvalue2 = os.system(python + " -c 'import {0}'".format(importname))
        if tryvalue2 != 0:
            print "{0}{2} failed its install (import returned {1}).  Re-installing now, after cleaning.".format(importname,tryvalue,vers)
            os.chdir('{0}-{1}'.format(target,vers))
            os.system(package_install_prefixes[target] +" "+ interpreter + ' setup.py clean >> logfile')
            os.system('rm -r build/')
            if target == 'matplotlib':
                os.system(matplotlib_command.format(pybinpath=fullpath+'/bin/', pyrootpath=fullpath, pyver=python_version ))
            else:
                os.system(package_install_prefixes[target] +" "+ interpreter + ' setup.py build'+package_install_flags[target]+' >> logfile')
                os.system(package_install_prefixes[target] +" "+ interpreter + ' setup.py install >> logfile')
            os.system('chmod -R +w *')
            os.chdir('../')
        tryvalue3 = os.system(python + " -c 'import {0}'".format(importname))
        if tryvalue3 != 0:
            print "ERROR: {0} failed to install after a clean!".format(importname)
            import pdb; pdb.set_trace()
    
    # independent matplotlib check
    trypylab = os.system(python + " -c 'import pylab'")
    if trypylab != 0:
        # matplotlib failure
        import pdb; pdb.set_trace()


    print (interpreter+" -c 'import {0}'").format(','.join(package_import_names.values()))
    if (os.system((interpreter+" -c 'import {0}'").format(','.join(package_import_names.values())))) != 0:
        import pdb; pdb.set_trace()
        #raise Exception("Install failed!")

#os.system('rm -r pywcs-{wv}'.format(**versions))
#os.system('rm -r pyfits-{fv}'.format(**versions))
#os.system('rm -r matplotlib-{mv}'.format(**versions))
#os.system('rm -r numpy-{nv}'.format(**versions))



