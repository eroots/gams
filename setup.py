from setuptools import setup


setup(name='gams',
      version='0.0.84',
      description='GUI tool providing various transformations of magnetic data',
      url='http://github.com/eroots/gams',
      author='Eric Roots',
      author_email='eroots087@gmail.com',
      packages=['gams',
                'gams.GUI',
                'gams.utils',
                'gams.resources'],
      long_description=open('readme.rst').read(),
      scripts=[],
      entry_points={'console_scripts': ['gams = gams.GUI.gams_gui:main']},
      install_requires=['numpy',
                        'scipy',
                        'matplotlib',
                        'geoh5py',
                        'scikit-image',
                        'pyqt5'],
      include_package_data=True)
