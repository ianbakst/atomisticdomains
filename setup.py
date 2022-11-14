import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="atomistic_domains",  # Replace with your own username
    version="0.0.1",
    author="Ian N. Bakst, Ph.D.",
    author_email="ian.n.bakst@gmail.com",
    description="A package to create and manipulate atomistic domains.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    include_package_data=True,
    data_files=[("", [])]
)
