import setuptools
from setuptools import setup

setup(
    name='jk.hecke',
    version='0.0.1',
    description='Hecke algebra computations',
    packages=setuptools.find_packages(),
    author='Johan Kahrstrom',
    author_email='johan.kahrstrom@gmail.com',
    install_requires=[
        'conttest>=0.0.8',
        'nose>=1.3.7',
        'numpy>=1.15.2',
        'sortedcontainers>=2.0.5'
    ],
    include_package_data=True
)
