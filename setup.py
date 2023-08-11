from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Python Wrapper for yezhengSTAT/ADTnorm'

setup(  name="adtnormpy", 
        version=VERSION,
        author="Daniel Caron",
        author_email="<dpc2136@cumc.columbia.edu>",
        description=DESCRIPTION,
        packages=find_packages(),
        url = 'https://github.com/donnafarberlab/ADTnormPy.git',
        keywords=['ADTnorm'],
	    install_requires=['pandas','numpy','anndata','mudata','rpy2'],
        extras_require= {'pytest':['pytest']},
      classifiers= [])