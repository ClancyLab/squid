import setuptools

description = "A set of modules to aid in atomistic and molecular simulations."
long_description = open("README.rst", 'r').read().strip()

setuptools.setup(
    name='clancyLab-squid',
    version='2.0.5',
    author="Clancy Group",
    author_email="ClancyLabJHU@gmail.com",
    description=description,
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/clancylab/squid",
    include_package_data=True,
    install_requires=[
        'pyDOE',
        'numpy',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'chko=squid.console_scripts.chkDFT:chkDFT',
            'chkDFT=squid.console_scripts.chkDFT:chkDFT',
            'scanDFT=squid.console_scripts.scanDFT:scanDFT',
            'pysub=squid.console_scripts.pysub:pysub',
            'procrustes=squid.console_scripts.procrustes:procrustes',
        ]
    },
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
    ],
)
