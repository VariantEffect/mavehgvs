import setuptools
import sys

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = ["fqfa>=1.2.0"]
# fqfa requires backported dataclasses in Python 3.6
if sys.version_info.major == 3 and sys.version_info.minor == 6:
    requirements.append("dataclasses")

setuptools.setup(
    name="mavehgvs",
    version="0.1.0",
    author="Daniel Esposito and Alan F Rubin",
    author_email="alan.rubin@wehi.edu.au",
    description=(
        "Regular expression-based validation of HGVS variant strings for clinical genetics and genomics applications."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/VariantEffect/mavehgvs",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
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
