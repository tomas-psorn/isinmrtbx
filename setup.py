from distutils.core import setup

setup(name='isinmrtbx',
      version='0.5',
      description='Tools for routine MRI data manipulation',
      author='Tomas Psorn',
      author_email='tomaspsorn@isibrno.cz',
      url='https://www.isibrno.cz',
      # packages=['isinmrtbx', 'isinmrtbx.datatypes','isinmrtbx.inout','isinmrtbx.tools'],
      packages=['isinmrtbx.datatypes', 'isinmrtbx.inout' ,'isinmrtbx.tools'],
     )