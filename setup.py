from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name = "Chempy",
    version = 0.1,
    description = "chemical evolution modeling",
    long_description = readme(),
    author = "Jan Rybizki",
    author_email = "",
    url = "https://github.com/jan-rybizki/Chempy",
    packages = find_packages(),
    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'Operating System :: OS Independent',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering'
      ],
    zip_safe=False
)
