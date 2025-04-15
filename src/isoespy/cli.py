# handler function
import sys
from isoespy import isoespy_de
from isoespy import isoespy_ff
from isoespy import isoespy_edger
from isoespy import isoespy_makefa
from isoespy import isoespy_makegtf
from isoespy import isoespy_orfpred

def isoespy_de_entry():
    isoespy_de.main(sys.argv[1:])

def isoespy_ff_entry():
    isoespy_ff.main(sys.argv[1:])

def isoespy_edger_entry():
    isoespy_edger.main(sys.argv[1:])

def isoespy_makefa_entry():
    isoespy_makefa.main(sys.argv[1:])

def isoespy_makegtf_entry():
    isoespy_makegtf.main(sys.argv[1:])

def isoespy_orfpred_entry():
    isoespy_orfpred.main(sys.argv[1:])
