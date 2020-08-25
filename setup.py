import setuptools
import sys

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = [
    "fqfa>=1.1.0",
]
# require backported dataclasses in Python 3.6
if sys.version_info.major == 3 and sys.version_info.minor == 6:
    requirements.append("dataclasses")

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
    python_requires=">=3.6",
    install_requires=requirements,
    test_suite="tests",
)
