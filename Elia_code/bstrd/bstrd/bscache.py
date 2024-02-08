from time import sleep
from bsread import source, dispatcher
from .bsvar import BSVar
from .prodthread import ProdThread


FIXED_CHANNELS = {
    "pid",
    "toc"
}


class BSCache:

    def __init__(self, timeout=1000, maxsize=100):
        self.timeout = timeout
        self.channels = {}
        self.data = None
        self.pt = ProdThread(self.run, maxsize=maxsize)

    def __repr__(self):
        return str(self.data)


    def __iter__(self):
        return self

    def __next__(self):
        self.pt.start()
        self.data = data = self.pt.get()
        return data


    def run(self, running):
        timeout_counter = 0
        configs = self.channels.values()
        with source(channels=configs, receive_timeout=self.timeout) as src:
            while running.is_set():
                msg = src.receive()
                if msg is None:
                    timeout_counter += 1
                    print("receive timed out", timeout_counter)
                    continue
                data = repack(msg)
                if data:
                    data["toc"] = timeout_counter
                    yield data


    def start(self):
        self.pt.start()

        while self.data is None:
            print("dropping empty data")
            next(self)

        while not self.channels.keys() <= self.data.keys():
            missing = self.channels.keys() - self.data.keys() - FIXED_CHANNELS
            print("dropping data that is missing channel(s):", sorted(missing))
            next(self)


    def stop(self):
        self.pt.stop()


    def get_vars(self, names):
        if not isinstance(names, dict):
            names = {n: {} for n in names}

        new_chans = {}
        for name, kwargs in names.items():
            if name not in FIXED_CHANNELS and name not in self.channels:
                check_availability(name)
                cfg = make_channel_config(name, *kwargs)
                new_chans[name] = cfg

        if new_chans:
            print("add new channels", sorted(new_chans))
            self.add_vars(new_chans)

        return {n: BSVar(n, self) for n in names}


    def get_var(self, name, modulo=None, offset=None):
        if name not in FIXED_CHANNELS and name not in self.channels:
            check_availability(name)
            cfg = make_channel_config(name, modulo, offset)
            print("add new channel", name)
            self.add_var(name, cfg)
        return BSVar(name, self)


    def add_vars(self, chans):
        self.stop()
        self.channels.update(chans)
        self.start()

    def add_var(self, name, cfg):
        self.stop()
        self.channels[name] = cfg
        self.start()

    def rem_vars(self, names):
        self.stop()
        for n in names:
            del self.channels[n]
        self.start()

    def rem_var(self, name):
        self.stop()
        del self.channels[name]
        self.start()

    def clear_vars(self):
        self.stop()
        self.channels.clear()
        # cannot start without any channel

    def flush(self):
        self.pt.queue.clear()



def check_availability(name):
    if not is_available(name):
        raise ValueError(f"channel {name} is not available")

def is_available(name):
    available = get_available_channels()
    return name in available

def get_available_channels():
    channels = dispatcher.get_current_channels()
    return set(ch["name"] for ch in channels)



def make_channel_config(name, modulo=None, offset=None):
    res = {}
    if modulo is not None:
        res["modulo"] = modulo
    if offset is not None:
        res["offset"] = offset
    if not res:
        return name
    res["name"] = name
    return res



def repack(message):
    data     = message.data.data
    pulse_id = message.data.pulse_id

    res = {n: v.value for n, v in data.items()}
#    res = {k: v for k, v in res.items() if v is not None}

    #TODO: should this be a ValueError?
    if any(v is None for v in res.values()):
        return None

    if not res:
        return None

    res["pid"] = pulse_id
    return res



