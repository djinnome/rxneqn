import io
from os.path import dirname, join
from setuptools import setup

def get_version(relpath):
  '''Read version info from a file without importing it'''
  for line in io.open(join(dirname(__file__), relpath), encoding='cp437'):
    if '__version__' in line:
      if '"' in line:
        # __version__ = "0.9"
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]


setup(name='rxneqn',
      version=get_version('rxneqn/__init__.py'),
      description='Balance chemical reactions and perform arithmetic on them',
      url='http://github.com/djinnome/rxn_balancer',
      author='Jeremy Zucker',
      author_email='djinnome@gmail.com',
      install_requires=[
       'pandas',
       'periodic'],
      license='MIT',
      packages=['rxneqn'],
      zip_safe=False)
