from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='MultiMM',  # Package name
    version='1.0.0',  # Version of the software
    description='A tool for chromatin modeling from nucleosomes to chromosomal territories.',
    author='Sebastian Korsak',
    author_email='s.korsak@datascience.edu.pl',
    url='https://github.com/SFGLab/MultiMM',  # GitHub repository URL
    license='GNU General Public License v3.0',
    packages=find_packages(),  # Automatically finds all packages and sub-packages
    install_requires=[  # List your package dependencies here
        'numpy',
        'scipy',
        'pandas',
        'argparse',
        'matplotlib',
        'mdtraj',
        'seaborn',
        'scikit-learn',
        'configparser',
        'typing-extensions',
        'pyBigWig',
        'hilbertcurve',
        'tqdm',
        'pyvista[all]',
        'OpenMM'
    ],
    entry_points={
        'console_scripts': [
            'MultiMM=MultiMM.run:main',  # MultiMM command points to run.py's main function
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Linux Debian, Fedora, Arch, Ubuntu etc.',
    ],
    python_requires='>=3.10',  # Specify Python version compatibility
)
