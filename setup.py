from setuptools import setup, find_packages
from pathlib import Path

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='MultiMM',  # Package name
    version='1.0.15',  # Version of the software
    description='A tool for chromatin modeling from nucleosomes to chromosomal territories.',
    long_description=long_description,
    long_description_content_type='text/markdown',
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
        'OpenMM',
    ],
    entry_points={
        'console_scripts': [
            'MultiMM=simulation.run:main',  # MultiMM command points to run.py's main function
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX :: Linux',  # General OS classifier
    ],
    python_requires='>=3.10',  # Specify Python version compatibility
)
