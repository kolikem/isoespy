from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize([
        "src/isoespy/isoespy_de.pyx",
        "src/isoespy/isoespy_ff.pyx",
        "src/isoespy/isoespy_edger.pyx",
        "src/isoespy/isoespy_makefa.pyx",
        "src/isoespy/isoespy_makegtf.pyx",
        "src/isoespy/isoespy_orfpred.pyx",
    ])
)

