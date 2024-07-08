from collections import OrderedDict


class LimitedDict(OrderedDict):

    def __init__(self, *args, maxlen=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.maxlen = maxlen
        self._ensure_length()

    def __setitem__(self, *args):
        super().__setitem__(*args)
        self._ensure_length()

    def _ensure_length(self):
        if self.maxlen is None:
            return
        while len(self) > self.maxlen:
            self.popitem(last=False)



