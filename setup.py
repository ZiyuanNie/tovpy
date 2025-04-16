from setuptools import setup, find_packages

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
    url="https://github.com/computationalrelativity/tovpy",
    package_dir={'tovpy': 'src/tovpy'},
    packages=find_packages(where='src'),
    package_data={"tovpy": ["eos/*"]},
    install_requires=['sys', 'os', 'shutil',
                      'numpy', 'scipy',
                      'matplotlib'], 
    scripts=[],
)
