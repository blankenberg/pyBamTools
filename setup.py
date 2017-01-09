import sys

if sys.version_info < (2, 5):
    print >> sys.stderr, "ERROR: pyBamTools requires python 2.5 or greater"
    sys.exit()

try:
    import setuptools
except ImportError:
    # Automatically download setuptools if not available
    from distribute_setup import use_setuptools
    use_setuptools()

from setuptools import setup, find_packages
from glob import glob

extra = {}
if sys.version_info >= (3,):
    extra['use_2to3'] = True

       
def main():
    setup(  name = "pyBamTools",
            version = "0.0.1",
            packages = find_packages( 'lib' ),
            package_dir = { '': 'lib' },
            scripts = glob( "scripts/*.py" ),
            setup_requires = ['numpy'],
            author = "Daniel Blankenberg",
            author_email = "dan.blankenberg@gmail.com",
            description = "Tools for parsing BAM data",
            url = "http://add_URL_HERE",
            zip_safe = False,
            dependency_links = [],
            **extra )

if __name__ == "__main__":
    main()
