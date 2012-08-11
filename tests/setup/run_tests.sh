pyexe=/Users/adam/repos/jenkins_tests/python$python-pywcs$pywcs-numpy$numpy-matplotlib$matplotlib-pyfits$pyfits/bin/python
if [ -d private ]; then
    cd private
    hg up
    cd ..
else
    hg clone ssh://hg@bitbucket.org/pyspeckit/pyspeckit-private private
fi
$pyexe private/boto_download.py --includestr fits
$pyexe setup.py install
cd tests
$pyexe run_tests.py
#$pyexe private/boto_upload.py --includesubstrdir tests_

