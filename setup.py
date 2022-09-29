import setuptools


setuptools.setup(
	name="Liftoff",
	version="1.6.3",
	author="Alaina Shumate",
	author_email="ashumat2@jhmi.edu",
	description="A gene annotation mapping tool",
	url="https://github.com/ashumate/Liftoff",
	install_requires=['numpy==1.22.0', 'biopython==1.76', 'gffutils==0.10.1', 'networkx==2.4', 'pysam==0.19.1','pyfaidx==0.5.8','interlap==0.2.6', "ujson==5.4.0","parasail==1.2.1"],
	python_requires='>=3.6',
	packages=['liftoff'],
	entry_points={'console_scripts': ['liftoff = liftoff.run_liftoff:main'], },
)
