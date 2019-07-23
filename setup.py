from setuptools import setup
import versioneer
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

setup(
    name='OTSun',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=['otsun',],
    url='https://github.com/bielcardona/OTSun',
    install_requires=['numpy==1.12', 'enum34' , 'autologging'],
    license='MIT',
    author='Gabriel Cardona, Ramon Pujol',
    author_email='gabriel.cardona@uib.es, ramon.pujol@uib.es',
    description='Analizer of sun collectors',
    long_description=long_description,
    long_description_content_type='text/markdown'
)

