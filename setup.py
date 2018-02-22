from setuptools import setup

setup(name='PGMF',
      version='1.0.0',
      description='predict gene model frames for multiple tissues, sexes, and species',
      url='http://github.com/haiwangyang/PGMF',
      author='Haiwang Yang',
      author_email='haiwangyang@gmail.com',
      license='MIT',
      packages=['PGMF'],
      install_requires=['urllib.request', 'argparse', 'pathlib', 'pyfaidx', 'Bio.Seq', 'pandas'],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False,
      entry_points={
          'console_scripts':
          [
              'downloadbig = PGMF.script.downloadbig:main',
        ],
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
