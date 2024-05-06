from setuptools import setup

setup(
    name='plastic_spectra',
    url='https://github.com/kburmak/plastic_mass_spectra',
    author='Karina Burmak',
    author_email='ksburmak@edu.hse.ru',
    packages=['plastic_spectra'],
    install_requires=['sqlite3', 'pandas', 'matplotlib', 'scipy'],
    version='0.1',
    license='Karina Burmak',
    description='A mass mass spectrum analysis package, which will help you to find plastics peaks (now PSS and PFAS)',
    long_description=open('README.txt').read(),
)