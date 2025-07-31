import os
from setuptools import setup, find_packages

# get long_description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

# list of all scripts to be included with package
scripts = [os.path.join('scripts',f) for f in os.listdir('scripts') if f.endswith('.py') or f.endswith('.sh')]

scripts = [os.path.join('scripts',f) for f in os.listdir('scripts') if not (f[0]=='.' or f[-1]=='~' or os.path.isdir(os.path.join('scripts', f)))]

print('installing scripts:'+str(scripts))
setup(
    name='CS2_proc',
    version='1.0.0.0',
    description='package to estimate elevations based on Cryosat-2 levei-1 data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='',
    author='Ben Smith', 
    author_email='besmith@uw.edu',
    license='MIT',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    keywords='Cryosat-2, altimetry, ice, interferometry',
    packages=find_packages(),
    scripts=scripts,
)
