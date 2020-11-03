"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='metagenlab_libs',

    version='1.0',
    description='Various utilities to generate circos plots, ete3 plots, jvenn plots',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/metagenlab/metagenlab_libs',
    author='MetaGenLab',
    author_email='trestan.pillonel@gmail.com',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # These classifiers are *not* checked by 'pip install'. See instead
        # 'python_requires' below.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    keywords='bacteria genome circos ete3 jvenn',

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    python_requires='>=3.0',

    install_requires=['matplotlib', 
                      'ete3>=3', 
                      'biopython>=1.77', 
                      'argparse']
                      #'PyQt5'],
    
    project_urls={
        'Bug Reports': 'https://github.com/metagenlab/metagenlab_libs/issues',
        'Source': 'https://github.com/metagenlab/metagenlab_libs',
    },
)
