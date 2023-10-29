from setuptools import setup, find_packages

from os import system, chdir
from pathlib import Path

#check if maturin is installed, pip install it if not
try:
    from maturin import __version__
except ImportError:
    system('pip install maturin[patchelf]')
print('\n\nNow to do the setup\n\n')
setup(
    name='polysolver',
    version='1.0',
    packages=find_packages(),
    url='https://github.com/tedjohnson12/astr731_polytrope',
    author='Ted Johnson',
    author_email='ted.johnson@unlv.edu',
)
print('\n\nNow to do the maturin stuff...\n\n')
original_dir = Path.cwd()
rs_path = Path(__file__).parent

chdir(rs_path)
system('maturin build --release')
chdir(original_dir)
