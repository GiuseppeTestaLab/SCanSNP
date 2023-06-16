from setuptools import setup, find_packages


requirements=['pandas>=1.1.5',
'pysam>=0.16.0.1',
'numpy>=1.19.5',
'scanpy>=1.7.2',
'scikit-learn>=0.23.1',
'scipy>=1.5.4',
'rpy2>=3.4.2',
'anndata>=0.7.6']


setup(
	name='SCanSNP',

	url='https://github.com/GiuseppeTestaLab/SCanSNP',

	author='Davide Castaldi',
	author_email='davide.castaldi@fht.org',

	license='GNU GENERAL PUBLIC LICENSE',
	
	entry_points={
		'console_scripts': [
			  'SCanSNP = SCanSNP.SCanSNP:main',
			  ],
		  }, 


	# You can just specify the packages manually here if your project is
	# simple. Or you can use find_packages().
	packages=find_packages(),
	scripts=['SCanSNP/fitdistrplus/fitdist.R', 'SCanSNP/fitdistrplus/quantile.R', 'SCanSNP/ModelFitting/ModelFitting.R'],
	zip_safe=False,
	include_package_data=True,

	# List run-time dependencies here.  These will be installed by pip when
	# your project is installed. For an analysis of "install_requires" vs pip's
	# requirements files see:
	# https://packaging.python.org/en/latest/requirements.html
	
	install_requires=requirements,

)