import setuptools


setuptools.setup(
	name="Liftoff",
	version="1.1.1",
	author="Alaina Shumate",
	author_email="ashumat2@jhmi.edu",
	description="A gene annotation mapping tool",
	url="https://github.com/ashumate/Liftoff",
	install_requires=['numpy', 'biopython', 'gffutils', 'networkx', 'pysam'],
	python_requires='>=3.0',
	packages=['liftoff'],
	entry_points={'console_scripts': ['liftoff = liftoff.liftoff:main'], },
)
