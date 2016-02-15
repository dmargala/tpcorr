#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup, Command
except ImportError:
    from distutils.core import setup, Command

# Run pre-built py.tests as described at
# https://pytest.org/latest/goodpractises.html#integrating-with-distutils-python-setup-py-test
class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        import sys
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

with open('DESCRIPTION.rst') as f:
    long_description = f.read()

requirements = [
    'numpy>=1.9.2',
    'astropy>=1.0.1',
    'matplotlib>=1.5.1',
    'scipy>=0.15.1',
    # 'bossdata>=0.2.9dev',
    # 'specsim>=0.2.dev283',
    # 'h5py>=2.5.0',
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='tpcorr',
    version='0.1.1',
    description='Throughput correction code for offset fibers in SDSS.',
    long_description=long_description,
    author='tpcorr developers',
    author_email='dmargala@uci.edu',
    url='https://github.com/dmargala/tpcorr',
    packages=[
        'tpcorr',
    ],
    package_dir={'tpcorr':
                 'tpcorr'},
    scripts = [
        'bin/tpcorr',
    ],
    #include_package_data=True,
    #zip_safe=False,
    install_requires=requirements,
    license='MIT',
    keywords='tpcorr',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    cmdclass = {'test': PyTest},
    #test_suite='tests',
    #tests_require=test_requirements
)