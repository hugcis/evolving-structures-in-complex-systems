import collections

Options = collections.namedtuple('Options', 'radius timesteps n_hid')

class DataItem:
    def __init__(self, options: Options, epochs,
                 train=None, test=None, fisher=None,
                 lookup_train=None, lookup_test=None):
        self.n_hid = options.n_hid
        self.epochs = epochs
        self.horizon = options.radius
        self.timesteps = options.timesteps
        if train is not None:
            self.train = train
        if test is not None:
            self.test = test
        if fisher is not None:
            self.fisher = fisher
        if lookup_train is not None:
            self.lookup_train = lookup_train
        if lookup_test is not None:
            self.lookup_test = lookup_test

    def is_metadata(self, options: Options, epochs):
        return (self.n_hid == options.n_hid
                and self.horizon == options.radius
                and self.timesteps == options.timesteps)

    def __str__(self):
        return DataItem.make_repr(
            Options(self.horizon, self.timesteps, self.n_hid), self.epochs)

    def __repr__(self):
        return "Data Item: " + str(self)

    @classmethod
    def make_repr(cls, options: Options, epochs):
        return "h{} r{} t{}".format(options.n_hid, options.horizon,
                                    options.timesteps)
