import warnings
from distutils.core import setup, Extension
warnings.filterwarnings("ignore")
setup(name="foo", version="1.0", \
     ext_modules=[Extension('foo', ['src/foo.c', 'src/common.c', 'src/index.c', 'src/predict.c'])])