import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = ["fqfa>=1.2.1"]

setuptools.setup(
    name="mavehgvs",
    version="0.5.0",
    author="MaveDB Developers",
    author_email="alan.rubin@wehi.edu.au",
    description=(
        "Regular expression-based validation of HGVS-style variant strings for Multiplexed Assays of Variant Effect."
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
