
class BSVar:

    def __init__(self, name, cache):
        self.name = name
        self.cache = cache

    def __repr__(self):
        return f"{self.name} = {self.value}"


    def get(self):
        return self.cache.data.get(self.name)

    value = property(get)


    def __iter__(self):
        return self

    def __next__(self):
        next(self.cache)
        return self.get()



