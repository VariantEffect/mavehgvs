import setuptools

setuptools.setup(
    name="hgvsp",
    version="0.4",
    packages=["hgvsp"],
    url="https://github.com/VariantEffect/hgvs-patterns.git",
    author="Daniel Esposito and Alan F Rubin",
    author_email="alan.rubin@wehi.edu.au",
    description=(
        "Regular expression-based validation of HGVS variant strings for clinical genetics and genomics applications."
    ),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    test_suite="tests",
)
