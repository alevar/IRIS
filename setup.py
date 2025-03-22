from setuptools import setup, find_packages

setup(
    name='iris-av',
    version='1.0.0',
    author='Ales Varabyou',
    author_email='ales.varabyou@jhu.edu',
    url='https://github.com/alevar/IRIS',
    description='IRIS: Detection and Validation Of Chimeric Reads',
    license='GPLv3',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'intervaltree',
        'pyfaidx',
        'setuptools'
    ],
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'iris=iris.core:main',
        ],
    },
)