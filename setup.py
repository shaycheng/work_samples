from setuptools import setup

with open('README.md') as f:
    long_description = f.read()

with open('requirementsSetup.txt') as f:
    requirements = f.read().splitlines()

setup(name='mydemo',
      version='0.0.1',
      description='Demo project',
      long_description=long_description,
      long_description_content_type='text/markdown',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: MIT License',
          'Operating System :: Unix',
          'Operating System :: MacOS',
          'Programming Language :: Python :: 3'
      ],
      url='https://github.com/shaycheng/work_samples',
      author='ShayCheng',
      author_email='shuitingc@gmail.com',
      keywords='Magnetic Resonance Imaging MRI analysis demo',
      packages=['mydemo'],
      install_requires=[],
      include_package_data=True,
      zip_safe=False)
