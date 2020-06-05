import setuptools

setuptools.setup(
    name='hgvsp',
    version='0.4',
    packages=['hgvsp'],
    url='https://github.com/FowlerLab/hgvs-patterns.git',
    license='MIT',
    author='Daniel',
    author_email='esposito.d@wehi.edu.au',
    description=(
        'A python utility containing HGVS Regex patterns to '
        'match a subset of the HGVS standard.'
    ),
    keywords=['hgvs', 'variants', 'regex', 'bioinformatics', 'biology'],
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    test_suite="tests",
)
