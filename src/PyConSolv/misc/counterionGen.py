from .parametrizerInterface import Parametrizer


class counterionParametrizer(Parametrizer):
    def __init__(self, structurePath: str):
        stype = 'CTI'
        Parametrizer.__init__(self, structurePath, stype)


