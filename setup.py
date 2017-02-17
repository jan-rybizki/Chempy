from setuptools import setup, find_packages
import os

result = [os.path.join(dp, f) for dp, dn, filenames in os.walk('Chempy/input') for f in filenames if (os.path.isfile(os.path.join(dp, f)) & ('~' not in f))]

result1 = [os.path.join(dp, f) for dp, dn, filenames in os.walk('Chempy/mcmc') for f in filenames if (os.path.isfile(os.path.join(dp,f)) & ('~' not in f))]

result += result1

for i,item in enumerate(result):
	result[i] = item[7:]

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
    package_dir = {'Chempy' : 'Chempy'},
    package_data = {'Chempy' : result},
    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'Operating System :: OS Independent',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering'
      ],
    zip_safe=False
)
