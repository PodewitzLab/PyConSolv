from .parametrizerInterface import Parametrizer


class solventParametrizer(Parametrizer):
    def __init__(self, structurePath: str):
        stype = 'SLV'
        Parametrizer.__init__(self, structurePath, stype)