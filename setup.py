from setuptools import setup

setup(
    name='plastic_spectra',
    url='https://github.com/kburmak/plastic_mass_spectra',
    author='Karina Burmak',
    author_email='ksburmak@edu.hse.ru',
    packages=['plastic_spectra'],
    install_requires=['pandas', 'matplotlib', 'scipy'],
    version='0.1',
    license='Karina Burmak',
    description='Python library designed to provide a workflow for plastic mass spectrometry data annotation, including creation of database, fragment ions generation, filling database, plastic annotation, and plotting resulting spectra.',
    long_description=open('README.txt').read(),
)
