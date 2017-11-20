__all__ = []
from hetio.preproc.BindingDB.dhimmel import preproc  
# import preproc 
from ipdb import set_trace 


def load():
    return preproc.load() 


def test_load():
    from hetio.preproc.BindingDB import dhimmel 
    dataset = dhimmel.load()

    set_trace()