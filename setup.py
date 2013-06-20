from distutils.core import setup

setup(
    name='altimpy',
    version='0.1.0',
    author='Fernando Paolo',
    author_email='fspaolo@gmail.com',
    packages=['altimpy', 'altimpy.test'],
    #scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    url='http://#',
    license='LICENSE.txt',
    description='Set of tools for processing satellite altimetry data',
    long_description=open('README.txt').read(),
    install_requires=[
        "PyTables >= 2.3.1",
    ],
)
