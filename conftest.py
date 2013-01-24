collect_ignore = ["setup.py", "debug_test.py",
        "build",
        "*ipython*",
        "doc/sphinxext",
        "doc/_build",
        "doc/_static",
        "doc/_templates",
        "doc/conf.py",
        "private",
        "tests/setup/run_tests.sh",
        "tests/setup/setup_jenkins_virtual.py", # downloads - don't do it.
        "tests/setup/test_imports.py",
        "tests/*ipython*",
        ".hg",
        ".git"
        ]

import warnings
warnings.filterwarnings("ignore",
    category=RuntimeWarning)

#    WARNING: RuntimeWarning: invalid value encountered in divide [image_registration.fft_tools.convolve_nd]
