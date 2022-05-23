import setuptools

readme = open('README.md').read()
requirements = open('requirements.txt').read().splitlines()

setuptools.setup(
    name='darkelf',
    version='0.0.0',
    description='Calculating dark matter scattering and absorption rates with the energy loss functions (ELF)',
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Brian Campbell-Deem, Simon Knapen, Jonathan Kozaczuk, Tongyan Lin and Ethan Villarama',
    url='https://github.com/tongylin/DarkELF',
    packages=setuptools.find_packages(),
    install_requires=requirements,
    package_dir={'darkelf': 'darkelf'},
    package_data={'darkelf': ['data/Al2O3/*', 'data/SiC/*', 'data/C/*', 'data/GaAs/*', 'data/Si/*',
                              'data/SiO2/*', 'data/Al/*', 'data/Ge/*', 'data/ZnS/*', 'data/Xe/*',
                              'data/GaN/*',
                              'examples/data/*']},
    zip_safe=False)
