from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='TOVpy',
    version='1.0.0',
    description='Tolman-Oppenheimer-Volkoff solver in Python',
    license="GNU-GPLv3",
    long_description=long_description,
    author='S.Bernuzzi and others',
    author_email='sebastiano.bernuzzi@uni-jena.de',
    url="",
    packages=['tovpy'],  
    install_requires=['sys', 'os', 'shutil',
                      'numpy', 'scipy',
                      'matplotlib'], 
    scripts=[],
)
