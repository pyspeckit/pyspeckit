import os

files = [
    'doublet_example.html.php',
    'examples.html.php',
    'h2co_example.html.php',
    'n2hp_example.html.php',
    'nh3_example.html.php',
    'sdss_example.html.php',
    ]

for fn in files:
    os.system('php %s > %s' % (fn,fn.replace(".php","")))
