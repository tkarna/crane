"""
Simple data structures for holding data identified by list of attributes.

Tuomas Karna 2013-11-05
"""

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------


def uniqueList(seq):
    "Return list with all duplicates removed. Preserves ordering."""
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]


def filterTuples(allTuples, queries, exclude=False):
    """A generic routine for filtering a list of tuples.

    Returns all items in allTuples that match the queries.
    Query is a tuple indicating the desired attributes. For example,
    ('a',1,None) would match all cases where the first item is 'a' and the second
    is 1. None is used as a wildcard, so in this case all values for the last
    item are accepted. Setting exclude=False excludes the matching tuples from
    the list.
    """
    olist = list()
    if not exclude:
        for query in queries:
            for tup in allTuples:
                if tupleMatches(tup, query):
                    olist.append(tup)
        olist = uniqueList(olist)
    else:
        olist = uniqueList(allTuples)
        for query in queries:
            for tup in allTuples:
                if tupleMatches(tup, query):
                    olist.remove(tup)
    return list(olist)


def tupleMatches(tup, query):
    """Compares a tuple to a query tuple. Query tuple has the same size as the
    tuple. None is used as a wildcard; Items with None are not tested and any
    value is accepted.
    """
    # find non-None keys for query
    nonZeroIx = [i for i in range(len(query)) if query[i]]
    # compare tuples, here's the one-liner!
    return all(v[0] == v[1]
               for i, v in enumerate(zip(tup, query)) if i in nonZeroIx)

#------------------------------------------------------------------------------
# Classes
#------------------------------------------------------------------------------


class tupleList(object):
    """Stores a list of tuples in a list. Each element is associated
    with a unique keyword.

    Supports filtering for keywords and basic combinations.

    Examples:
    # constructor
    tl = tupleList(['species','color','number_of_legs']
    # adding
    tl.addSample( species='dog', color='red', number_of_legs=3 )
    tl.addSample( ('dog','black',4) ) # equivalent to above
    # filtering
    matching_keys = tl.getKeys( species='dog' )
    matching_keys = tl.getKeys( color='black' )
    # get subset
    tl2 = tl.getSubset( color='black' )

    """

    def __init__(self, keywords, source=None):
        """Create an empty database with the given keywords.

        Params
        ------
        keywords : list of str
            Unique names to identify each member in the tuple
        source : iterable, optional
            If provided, all tuples in the iterable will be added to this object

        """
        if len(keywords) > len(set(keywords)):
            raise Exception('keywords must be unique')
        self.keywords = keywords
        self.data = []
        if source is not None:
            for i in source:
                self.addSample(tup=i)

    def __iter__(self):
        return self.data.__iter__()

    def __len__(self):
        return self.data.__len__()

    def addSample(self, tup=None, **keywords):
        """Adds a tuple in the list if not already present.

        Sample can be given as a tuple or with the keyword arguments.

        Examples:
        tl = tupleList(['species','color','number_of_legs']
        tl.addSample( ('dog','black',3) )
        tl.addSample( species='cat',color='gray' )
        """
        if tup is None and not keywords:
            raise Exception('tuple or keywords must be given')
        if tup is None:
            tup = tuple([keywords.get(k, None) for k in self.keywords])
        if len(tup) != len(self.keywords):
            msg = 'Given tuple has incorrect length: %d != %d' % (
                len(tup), len(self.keywords))
            raise Exception(msg)
        if self.data.count(tup) == 0:
            self.data.append(tup)

    def getKeys(self):
        """Return all tuples in a list"""
        return self.data

    def genKey(self, **kwargs):
        """Given the kwargs, generates a tuple of keywords.
        Missing keywords are replaced by None."""
        return tuple([kwargs.get(k, None) for k in self.keywords])

    def getTuples(self, query=None, exclude=False, **kwargs):
        """Returns a list of all keys that match the query.
        query is a list of dictionaries, with appropriate keywords.
        Alternatively, keywords can be given as argument to this function.
        In this case one keyword can be a list.
        Setting exclude=True negates the query.
        """
        if query is None:
            # kwargs is a single dict, one value may be list
            nListArgs = sum(isinstance(v, list) for v in kwargs.values())
            if nListArgs > 1:
                raise Exception('Only one of kwargs can be a list')
            if nListArgs:
                # construct tuples
                for k in kwargs.keys():
                    if isinstance(kwargs[k], list):
                        break
                tuples = []
                for i in range(len(kwargs[k])):
                    # construct dicts with one item of the list in k
                    d = dict(kwargs)
                    d[k] = kwargs[k][i]
                    tuples.append(self.genKey(**d))
            else:
                tuples = [self.genKey(**kwargs)]
        else:
            # filter is a list of dict
            # convert to list of tuples
            tuples = [self.genKey(**d) for d in query]
        return filterTuples(self.data, tuples, exclude)

    def getSubset(self, query=None, exclude=False, **kwargs):
        """Returns a copy of this object, where members are chosen by the given
        filter."""
        return tupleList(
            self.keywords, self.getTuples(
                query, exclude, **kwargs))

    def update(self, other):
        """Append tuples from other tupleList to this object."""
        for t in other.getTuples():
            self.addSample(tup=t)
