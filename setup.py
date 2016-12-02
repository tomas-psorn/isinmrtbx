from setuptools import setup

setup(name='isinmrtbx',
      version='0.6.3',
      description='Tools for routine MRI data manipulation',
      author='Tomas Psorn',
      author_email='tomaspsorn@isibrno.cz',
      url='https://www.isibrno.cz',
      packages=['isinmrtbx.datatypes', 'isinmrtbx.inout' ,'isinmrtbx.tools'],
      install_requires=['numpy', 'matplotlib', 'scipy'],
      data_files=[('examples',['isinmrtbx/examples/data/ute/*'])],
      license='GPL',
      zip_safe=False
     )