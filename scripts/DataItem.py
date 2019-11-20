class DataItem:
    def __init__(self, n_hid, epochs, horizon, timesteps,
                 train=None, test=None, fisher=None,
                 lookup_train=None, lookup_test=None):
        self.n_hid = n_hid
        self.epochs = epochs
        self.horizon = horizon
        self.timesteps = timesteps
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

    def is_metadata(self, n_hid, epochs, horizon, timesteps):
        return (self.n_hid == n_hid
                and self.horizon == horizon
                and self.timesteps == timesteps)

    def __str__(self):
        return DataItem.make_repr(self.n_hid, self.epochs,
                                  self.horizon, self.timesteps)

    def __repr__(self):
        return "Data Item: " + str(self)

    @classmethod
    def make_repr(cls, n_hid, epochs, horizon, timesteps):
        return "h{} r{} t{}".format(n_hid, horizon, timesteps)
