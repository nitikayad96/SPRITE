from setuptools import setup, find_packages

setup(
    name='sprite',
    version='0.1.0',
    packages=find_packages(),
    scripts=['sprite_utils.py', 'red_pitaya_utils.py']
)