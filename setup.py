from setuptools import setup
import os

setup(name = 'logan_blaster',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author = 'Pierre Peterlongo',
    py_modules=['logan_blaster'],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'logan_blaster=logan_blaster:main',
        ],
    },
    author_email = 'pierre.peterlongo@inria.fr',
    maintainer = logan_blaster.__author__,
    maintainer_email = 'pierre.peterlongo@inria.fr',
    download_url = 'https://github.com/pierrepeterlongo/logan_blaster',
    url = 'https://github.com/pierrepeterlongo/logan_blaster',
    description = 'Align genomic sequences with Logan contigs or unitigs.',
    long_description = open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
    license = 'GPL-3.0'
    )

