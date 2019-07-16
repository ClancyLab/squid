import setuptools

description = "A set of modules to aid in atomistic and molecular simulations."
long_description = open("README.rst", 'r').read().strip()

setuptools.setup(
    name='squid',
    version='2.0.0',
    author="Clancy Group",
    author_email="ClancyLabJHU@gmail.com",
    description=description,
    long_description=long_description,
    long_description_content_type="text/rst",
    url="https://github.com/clancylab/squid",
    include_package_data=True,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CC BY-NC 3.0 US",
        "Operating System :: OS Independent",
    ],
)
