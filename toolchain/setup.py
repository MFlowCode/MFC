import setuptools

setuptools.setup(
    name='mfc',
    version='0.1',
    license='MIT',
    zip_safe=False,
    packages=setuptools.find_packages(),
    fullname='MFC Toolchain',
    author_email='hberre3@gatech.edu',
    url='https://github.com/MFlowCode/MFC',
    install_requires=[
        'argparse', 'dataclasses', 'typing',
        'pyyaml',   'rich',        'fypp'
    ],
    description='Toolchain for a fully-documented parallel simulation software for multi-component, multi-phase, and bubbly flows.',
)
