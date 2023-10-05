from setuptools import setup, find_packages

setup(
    name='tiepy',
    version='0.1.0',
    author = "Nikhil Rao",
    author_email = "nikhilr@ucar.edu",
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=[
        'cartopy',
        'matplotlib',
        'numpy',
        'xarray'
    ],
    entry_points={
        'console_scripts': [
            'tiepy=src.tiepy:main',
            'getoptions=src.getoptions:main'
        ]
    }
)
