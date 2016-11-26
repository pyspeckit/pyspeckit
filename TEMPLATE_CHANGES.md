This file keeps track of changes between tagged versions of the Astropy
package template, for the benefit of affiliated package maintainers. It can
be removed in affiliated packages.

The changes below indicate what file the change was made in so that these can
be copied over manually if desired.

1.1.3 (unreleased)
------------------

- Updated the setup.cfg to include the version number. [#129]

- Removed Python 2.6 tests from travis.yml file as astropy 1.2 no longer supports Python 2.6 [#183]

- Updated ``setup.cfg`` and ``setup.py`` so that the install requirements
  are defined in ``setup.cfg``. [#208]

1.1.2 (2016-07-02)
------------------

- Updated .travis.yml to show usage of SETUP_XVFB ci-helpers option. [#177]

- Updated astropy-helpers to v1.2. [#180]

- Fixed import of ``configparser`` in ``docs/conf.py``. [#180]

1.1.1 (2016-06-03)
------------------

- Fixed the import of configparser on Python 3.5. [#172]

- Updated ``ez_setup.py`` to the latest version. [#174]

1.1 (2016-06-01)
----------------

- Fixed the import of example_mod.py in __init__.py to work properly on
  Python 3. [#167]

- Fixed the version import in conftest.py to work properly if version.py
  hasn't been generated yet. [#167]

- Switch to using container-based builds on Travis. [#133]

- Entry points are now defined in the ``setup.cfg`` file. [#130]

- Expanded the default .gitignore file. [#146]

- Switch to using ci-helpers on Travis, and add example AppVeyor config file.
  [#140]

- Updated astropy-helpers to v1.1.2. [#147]

- Remove ``adjust_compiler`` from ``setup.py`` (this is now dealt with in
  astropy-helpers). [#154]

- Update ``MANIFEST.in`` to avoid bundling temporary files in astropy-helpers
  when releasing packages. [#154]

- Updated ``ez_setup.py`` to the latest version. [#135]

- Catch ``KeyError`` when setting up custom test headers in ``conftest.py``.
  [#143]

- Update ``.travis.yml`` to include testing with the LTS release of Astropy,
  remove testing against Python 2.6, and update Numpy versions [#148, #168]

- Add an example of how to include data in ``example_subpkg/setup_package.py``.
  [#158]

- Update ReadTheDocs domains. [#166]

- Updated Sphinx Makefile. [#171]

- Updated AppVeyor config to run on Python 2.7 and 3.5. [#170]

In summary, the following files should be updated by affiliated packages:

- ``.gitignore``
- ``MANIFEST.in``
- ``docs/Makefile``
- ``ez_setup.py``
- ``packagename/conftest.py``
- ``setup.py``
- ``setup.cfg`` (if using entry points)

In addition, astropy-helpers should be updated to v1.1.2 and the
``ah_bootstrap.py`` file should be updated to match the version in
``astropy-helpers``.

Finally, the following files have been updated, but don't necessarily need to
be updated since they will likely be heavily customized in affiliated packages:

- ``.travis.yml``
- ``appveyor.yml``

1.0 (2015-05-31)
----------------

- The instructions for the documentation have now been clarified to indicate
  that packages do not *have* to include the documentation in a sub-folder of
  ``docs`` (see updated note in ``docs/index.rst``). [#123]

- Updated ``setup.cfg`` to enable ``doctest_plus`` by default.

- Updated ``.travis.yml`` to:

  - Update apt-get package list

  - Add ``jinja2`` as a dependency to be installed with conda [#114]

  - Drop Python 3.2 testing [#114]

  - Drop Numpy 1.5 testing, and use Numpy 1.9 as a baseline [#114]

- Updated ``MANIFEST.in`` to:

  - Recursively include *.pyx, *.c, and *.pxd files

  - Globally exclude *.pyc and *.o files

  - Include ``CHANGES.rst``

- Update ``docs/conf.py`` to import Sphinx extensions from
  ``astropy_helpers`` instead of ``astropy``. [#119]

- Added 'Powered by Astropy badge' to ``README.rst``. [#112]

- Show how to add and remove packages from pytest header in
  ``packagename/conftest.py``, and show how to show the package version
  instead of the astropy version in the top line.

- Minor documentation change in ``packagename/_astropy_init.py``. [#110]

- Use setuptools entry_points for command line scripts (change in
  ``setup.py``). [#116]

- Updated ``astropy-helpers`` and ``ah_bootstrap.py`` to v1.0.2.

- Remove requires and provides from setup.py. [#128]

0.4.1 (2014-10-22)
------------------

- Changed order of exclusion in MANIFEST.in, excluding *.pyc *after* including
  astropy-template

- Updated astropy-helpers to v0.4.3

0.4 (2014-08-14)
----------------

- Initial tagged version, contains astropy-helpers v0.4.1
