import setuptools

setuptools.setup(
    name='BactInspectorMax',
    version='0.0.8',
    description='Package to investigate mash hits against refseq',
    author='Anthony Underwood',
    author_email='au3@sanger.ac.uk',
    license='MIT',
    packages=setuptools.find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'bactinspector = bactinspector.run_bactinspector:main'
        ]
    },
    install_requires=['pandas', 'numpy'],
    test_suite='nose.collector',
    tests_require=['nose'],
    classifiers=[ 
        'Development Status :: 3 - Alpha', 
        'Intended Audience :: Science/Research', 
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'
    ]
)