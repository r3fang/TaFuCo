from distutils.core import setup, Extension
setup(name="foo", version="1.0", \
     ext_modules=[Extension('foo', ['src/foo.c', 'src/uthash.c'])])