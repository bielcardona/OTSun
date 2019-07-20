from setuptools import setup
import versioneer

setup(
    name='OTSun',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=['otsun',],
    url='',
    install_requires=['numpy==1.12', 'enum34' , 'autologging'],
    license='MIT',
    author='Gabriel Cardona, Ramon Pujol',
    author_email='gabriel.cardona@uib.es, ramon.pujol@uib.es',
    description='',
)

