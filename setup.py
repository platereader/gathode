"""
GATHODE: Growth Analysis Tool for High-throughput Optical Density Experiments

GATHODE is a software package for analysing time series of optical
density measurements that were recorded with the help of a plate
reader. It allows to extract growth parameters such as maximal growth
rate, lag-time and growth yield.
"""

import sys
from setuptools import setup, find_packages

packagedir = 'platereader'
version_py = packagedir+'/_version.py'

def getVersion():
    """ get the version from the file _version.py """
    try:
        fh=open(version_py, 'r')
        version=fh.read().strip().split('=')[-1].replace("'",'').lstrip()
        fh.close()
    except:
        return None

    return version

def genericSetupDict(app='both'):
    genericSetupOpts=dict(
        name = app,
        version = getVersion(),
        packages = [packagedir, packagedir+'/parser'],
        # NOTE PyQt4 is not installable via pip, so this dependency is not listed here
        install_requires = ["numpy","scipy","matplotlib"],

        author = "Nils Christian",
        author_email = "nils.christian@uni.lu",
        url = "https://platereader.github.io/",
        license = "AGPL",
        classifiers=[
            'Intended Audience :: Science/Research',
            'Environment :: Console',
            'Environment :: MacOS X',
            'Environment :: Win32 (MS Windows)',
            'Environment :: X11 Applications :: Qt',
            'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            ],
        keywords = "plate reader, optical density, growth curve",
        )
    if app == 'both':
        genericSetupOpts['name'] = 'GATHODE'
        genericSetupOpts['description'] = """Growth and Chronological Life Span Analysis Tools for High-throughput
Optical Density Experiments (GATHODE/CATHODE)"""
        genericSetupOpts['long_description'] = """
The Growth Analysis Tool for High-throughput Optical Density
Experiments (GATHODE) is a software package for analysing time series
of optical density measurements that were recorded with the help of a
plate reader. It allows to extract growth parameters such as maximal
growth rate, lag-time and growth yield.
The Chronological Life Span is defined as the time cells can survive
in a non-dividing state. The Chronological life span Analysis Tool for
High-throughput Optical Density Experiments (CATHODE) uses multiple
output files of GATHODE to analyse this survival."""

    elif app == 'GATHODE':
        genericSetupOpts['description'] = "Growth Analysis Tool for High-throughput Optical Density Experiments (GATHODE)"
        genericSetupOpts['long_description'] = """
The Growth Analysis Tool for High-throughput Optical Density
Experiments (GATHODE) is a software package for analysing time series
of optical density measurements that were recorded with the help of a
plate reader. It allows to extract growth parameters such as maximal
growth rate, lag-time and growth yield."""

    elif app == 'CATHODE':
        genericSetupOpts['description'] = "Chronological life span Analysis Tool for High-throughput Optical Density Experiments (CATHODE)"
        genericSetupOpts['long_description'] = """
The Chronological Life Span is defined as the time cells can survive
in a non-dividing state. The Chronological life span Analysis Tool for
High-throughput Optical Density Experiments (CATHODE) uses multiple
output files of the Growth Analysis Tool for High-throughput Optical
Density Experiments (GATHODE) to analyse this survival."""

    else:
        raise RuntimeError('no such app "'+app+'"')

    return genericSetupOpts

def setupPackage():
    genericSetupOpts=genericSetupDict()
    setup(
        entry_points = {
            'console_scripts': [
                'gathodecli = platereader.odcli:odCommandlineInterface',
                'cathodecli = platereader.clscli:clsCommandlineInterface',
                ],
            'gui_scripts': [
                'gathode = platereader.odgui:gui_main',
                'cathode = platereader.clsgui:clsgui_main',
                ]
            },
        **genericSetupOpts
        )

def setupOsXapp(app):
    setupOpts=genericSetupDict(app)
    setupOpts['options']={
        'py2app': { 'argv_emulation': True },
        'plist': {},
        }
    if app == 'GATHODE':
        setupOpts['app']=['platereader/odgui.py']
    elif app == 'CATHODE':
        setupOpts['app']=['platereader/clsgui.py']
    setup(
        **setupOpts
        )

if __name__ == "__main__":
    if len(sys.argv)<=1 or sys.argv[1] == 'version':
        print(getVersion())
    elif len(sys.argv)>1 and sys.argv[1] == 'py2app':
        app='GATHODE'
        if len(sys.argv) > 2:
            app=sys.argv.pop()
        setupOsXapp(app)
    else:
        setupPackage()
