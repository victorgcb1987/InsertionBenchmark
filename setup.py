from setuptools import setup
from os.path import join, dirname

setup(
    name='InsertionBenchmark',
    version='1.0',
    packages=["./"],
    long_description=open(join(dirname(__file__), 'README.md')).read(),
)
