from setuptools import setup, find_packages

setup(
    name='my_radio_astronomy_project',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'numpy',
        'matplotlib',
        'astropy'
    ],
    entry_points={
        'console_scripts': [
            'run-app=src.main:main'
        ]
    },
)
