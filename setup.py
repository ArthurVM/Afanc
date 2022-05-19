import subprocess
from glob import glob
from os import listdir
from os.path import splitext, basename, dirname, realpath, exists, join
from setuptools import setup, find_packages
from setuptools.command.install import install as DistutilsInstall

class Install(DistutilsInstall):

    def run(self):
        self._install_autodatabase()

    def _install_autodatabase(self):
        """ clone autodatabase repo from github
        """
        currant_dir = dirname(realpath(__file__))

        autodatabase_git_dir = join(currant_dir, "autodatabase")

        if not exists(autodatabase_git_dir):
            subprocess.call("git clone https://github.com/annacprice/autodatabase", cwd=currant_dir, shell=True)

        DistutilsInstall.run(self)

with open("README", 'r') as f:
    long_description = f.read()

setup(
    name='Afanc',
    version='0.1-alpha',
    description='General purpose high-resolution metagenomics disambiguation.',
    license="MIT",
    long_description=long_description,
    author='Arthur V. Morris',
    author_email='morrisa28@cardiff.ac.uk',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
    'Biopython==1.79',
    'numpy==1.20.3',
    'pysam==0.18.0',
    'pandas==1.3.4',
    'scipy==1.7.1'
    ],
    entry_points={
        'console_scripts': [
            'Afanc = Afanc.cli:main',
        ],
    },
    # cmdclass={'install': Install},
)
