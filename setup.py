from distutils.core import setup

setup(
    name='altimpy',
    version='0.1.0',
    author='Fernando Paolo',
    author_email='fspaolo@gmail.com',
    packages=['altimpy', 'altimpy.tests'],
    #scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    #url='https://pypi.python.org/pypi/altimpy',
    url='https://github.com/fspaolo/altimpy',
    download_url='https://github.com/fspaolo/altimpy/tarball/master',
    license='LICENSE.txt',
    description='Set of tools for processing satellite altimetry data',
    long_description=open('README.rst').read(),
    install_requires=[
        "PyTables >= 2.3.1",
        "basemap >= 1.0.7",
        "netcdf4 >= 1.0.8",
    ],
)
