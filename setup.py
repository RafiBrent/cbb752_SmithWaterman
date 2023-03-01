from setuptools import setup, find_packages

setup(
    name='SmithWaterman',
    version='0.1.0',
    packages=find_packages(include=['cbb752_WatermanSmith', 'SmithWaterman.*']),
)
