try:
    import pg
except:
    try:
        import _pg
        pg = _pg
    except:
        # could not import pg, database access will be dysfunctional
        pass
import sys
import re
import popen2
import time

import ConfigParser

#import cmop
# class loggerStub:
#'''A logging abstraction to make it easier to decouple this module from the cmop module.'''
# def debug(self, s, level=10):
#cmop.debug(s, level)

# def info(self, s):
# cmop.info(s)

# def profile(self, s):
# cmop.profile(s)

#logger = loggerStub()

DELIM = ';'


class SQLError(Exception):
    pass


class DuplicateError(SQLError):
    pass


def donothing(data, xxx_todo_changeme):
    (e, v, t) = xxx_todo_changeme
    pass


def commas(xs):
    return ', '.join([str(x) for x in xs])

#defaults = cmop.config


def escape(v):
    return v.replace("'", "''")

expr = re.compile("\s*('.*')|(\(.*\))\s*")


def quote(v):
    # format strings for postgres
    if isinstance(v, str) or isinstance(v, unicode):
        # add single quotes to strings,
        # unless they are already there
        if not expr.match(v):
            return "'%s'" % (escape(v),)
        else:
            return v
    if v is None:
        return 'NULL'
    else:
        return v


def equals(a, v):
    # build an SQL equailty expression
    # for attribute = value
    # floats are not necessarily equal
    if isinstance(v, float):
        return "trunc(%s::numeric,2)::numeric=trunc(%s,2)::numeric" % (a, v)
    else:
        return "%s=%s" % (a, v)

identpatt = re.compile('^[_a-zA-Z][_a-zA-Z0-9]*$')


def validIdentifier(name):
    '''Quotes the name if necessary, returning a valid SQL identifier.
 Behavior changed: There seems to be no harm in always quoting field names.
 '''

    return '"%s"' % (name,)
    # if identpatt.match(name):
    #  return name
    # else:
    #  return '"%s"' % (name,)

tblexpr = re.compile('(?:(.*?)\.)?(?:"(.*)"|(.*))')


def tablepair(table):
    '''return the schema and table name'''
    grps = tblexpr.match(table).groups()
    sch = grps[0]
    tab = grps[1] or grps[2]
    if not sch:
        sch = 'public'
    return sch, tab


def fqTableName(table):
    '''Return a fully qualified table name'''
    sch, tab = tablepair(table)
    return '%s."%s"' % (sch, tab)


class DB:
    '''
    A Database abstraction class
    --works with postgres only, but exists to protect higher-level
      apps from a db software change
    --Provides direct query (with results) and command (no results)
    --Provides abstraction of logging, configuration, and errorhandling for the database

    make sure user has permissions in pg_hba.conf on hostname
    '''

    def __init__(self, dbname=None, host=None,
                 user=None, password=None):
        self.dbname = dbname
        self.user = user
        self.host = host
        self.password = password
        self.dbconn = None
        self.connect()
        self.suppressquery = False
        self.schema = {}

    def execCommand(self, qry, retries=1):
        '''Executes a SQL statement, ignoring the results'''
        self.connect()
        #logger.debug('Executing command (%s): \n%s' % (self.host,qry), 7)

        if self.suppressquery:
            # cmop.info(qry)
            return

        try:
            result = self.dbconn.query(qry)
        except pg.ProgrammingError as e:
            if "duplicate" in str(e):
                raise DuplicateError(
                    '%s on %s: %s \n SQL="%s"' %
                    (self.dbname, self.host, e, qry))
            elif retries == 0:
                raise SQLError(
                    '%s on %s: %s \n SQL=%s' %
                    (self.dbname, self.host, e, qry))
            elif 'terminating' in str(e):
                # connection is bad
                self.dbconn = None
                self.execCommand(qry, retries=retries - 1)
            else:
                raise SQLError(
                    '%s on %s: %s \n SQL=%s' %
                    (self.dbname, self.host, e, qry))

    def close(self):
        if self.dbconn:
            self.dbconn.close()
            self.dbconn = None

    def cursor(self, qry, blocksize=100, cursorname='mycursor'):
        '''A generator for tuples for query qry using a cursor'''
        self.begin()
        sql = '''DECLARE %s CURSOR FOR %s''' % (cursorname, qry)
        self.execCommand(sql)
        fetch = '''FETCH %i FROM %s''' % (blocksize, cursorname)
        rs = self.execQuery(fetch)
        while rs:
            for r in rs:
                yield r
            rs = self.execQuery(fetch)
        self.commit()

    def reconnect(self):
        self.close()
        self.connect()

    def execQuery(self, qry, withfields=False, retries=1):
        '''
  Executes a SQL statement, returning the results
        '''
        self.connect()
        #logger.debug('Executing query (%s): \n%s' %(self.host,qry), 7)
        qry = qry.strip()

        if self.suppressquery:
            # cmop.info(qry)
            if withfields:
                return ([], [])
            else:
                return []

        if not len(qry):
            if withfields:
                return ([], [])
            else:
                return []
        try:
            response = self.dbconn.query(qry)
        except pg.ProgrammingError as e:
            if retries == 0:
                raise SQLError(
                    '%s on %s: %s \n SQL=%s' %
                    (self.dbname, self.host, e, qry))
            elif 'terminating' in str(e):
                # connection is broken, so retry
                self.dbconn = None
                return self.execQuery(qry, retries=retries - 1)
            else:
                raise SQLError(
                    '%s on %s: %s \n SQL=%s' %
                    (self.dbname, self.host, e, qry))

        if response:
            if withfields:
                return (response.listfields(), response.getresult())
            else:
                return response.getresult()

    def connect(self):

        dbname = self.dbname
        hostname = self.host
        user = self.user
        password = self.password

        reset = False

        # check connection validity
        # several ways to do this
        # needs to be cheap -- no roundtrip
        if not self.dbconn:
            reset = True
        elif self.dbconn.fileno() == -1:
            reset = True

        if reset:
            #logger.info("Establishing connection to %s:%s using %s" % (hostname, dbname, user))

            self.dbconn = pg.DB(dbname=dbname,
                                host=hostname,
                                port=-1,
                                user=user,
                                passwd=password)

    def CachedPrimaryKey(self, table):
        if table in self.schema:
            keys, attrs = self.schema[table]
        else:
            keys = self.PrimaryKey(table)
            attrs = self.Attributes(table)
            self.schema[table] = keys, attrs
        return keys

    def PrimaryKey(self, table):
        pair = table.split(".")
        if len(pair) == 1:
            sch, tab = 'public', pair[0]
        else:
            sch, tab = pair

        keysql = '''
select a.attname
  from pg_attribute a, pg_constraint c, pg_class t, pg_namespace s
 where contype = 'p'
   and conrelid = t.oid
   and t.relnamespace = s.oid
   and a.attnum = ANY (c.conkey)
   and a.attrelid = t.oid
   and s.nspname = '%s'
   and t.relname = '%s'
'''
        rs = self.execQuery(keysql % (sch, tab))
        return [a[0] for a in rs]

    def Attributes(self, table, namefilter=None):
        attrsql = '''
select attname
  from pg_attribute a,
       pg_class c,
       pg_namespace s
 where a.attrelid = c.oid
   and c.relnamespace = s.oid
   and s.nspname = '%s'
   and c.relname = '%s'
   and a.attnum > 0
   and a.atttypid != 0
'''
        sch, tab = tablepair(table)

        if namefilter:
            attrsql = (attrsql + '''
    and attname NOT LIKE '%s%%'
 ''') % (sch, tab, namefilter)
        else:
            attrsql = attrsql % (sch, tab)

        rs = self.execQuery(attrsql)
        return [a[0] for a in rs]

    def AllTables(self, sch='%'):
        '''Return all tables in all schemas matching regular expression sch.  Return type is a list of tuples each of the form (schema name, table name, table kind), where table kind is 'r' for relations and 'v' for views.'''
        sql = '''
SELECT nspname, relname, relkind
  FROM pg_class c, pg_namespace s
 WHERE s.oid=c.relnamespace
   AND (relkind = 'r' OR relkind = 'v')
   AND nspname NOT LIKE 'pg_%%'
   AND nspname NOT LIKE 'information_schema'
   AND nspname LIKE '%s'
''' % (sch,)
        return self.execQuery(sql)

    def AllSchemas(self):
        sql = '''
SELECT nspname
  FROM pg_namespace
 WHERE nspname NOT LIKE 'pg_%'
   AND nspname NOT LIKE 'information_schema'
'''
        return self.execQuery(sql)

    def QuoteAsString(self, s):
        return "'%s'" % (s,)

    def TableComment(self, tbl):
        if not self.checkTable(tbl):
            raise ValueError("Table %s not found." % (tbl,))

        # postgres only
        sch, tbl = tablepair(tbl)

        sql = '''SELECT obj_description(
        (SELECT c.oid FROM pg_class c, pg_namespace s
         WHERE c.relname='%s' and s.nspname='%s'
           AND c.relnamespace = s.oid), 'pg_class')
         AS comment;''' % (tbl, sch)
        rs = self.execQuery(sql)
        return rs[0][0]

    def SchemaComment(self, sch):
        '''Returns owner, comment pairs for the given schema'''
        # postgres only
        sql = """
SELECT r.usename, oid
  FROM pg_namespace s, pg_user r
 WHERE s.nspowner = r.usesysid
   AND s.nspname='%s'""" % (sch,)

        rs = self.execQuery(sql)
        if rs:
            sql = '''SELECT '%s', obj_description(%s, 'pg_namespace')''' % rs[
                0]
            rs = self.execQuery(sql)
            return rs[0]
        else:
            raise ValueError("Schema %s not found." % (sch,))

    def AttributeComment(self, tblname, attr):
        sch, tbl = tablepair(tblname)
        # postgres only
        sql = """
SELECT col_description(c.oid, attnum)
  FROM pg_namespace s, pg_class c, pg_attribute a
 WHERE s.nspname = '%s' AND c.relname = '%s' AND a.attname = '%s'
   AND c.relnamespace = s.oid
   AND a.attrelid = c.oid
""" % (sch, tbl, attr)

        rs = self.execQuery(sql)
        if rs:
            return rs[0][0]
        else:
            raise ValueError(
                "Attribute %s not found for table %s" %
                (tblname, attr))

    def appendTo(self, tablename, qry):
        '''
  Append a query's results to a table, creating it if it doesn't exist.
       '''
        if self.checkTable(tablename):
            qry = "INSERT INTO %s (%s)"
        else:
            qry = "CREATE TABLE %s AS (%s);" % (tablename, qry)

        self.execCommand(qry)

    def ValuesClause(self, tuples):
        '''Format a list of tuples as a single values clause
  for a multi-tuple insert statement.'''
        rs = tuples
        if not rs:
            #cmop.info("No tuples to transfer.")
            return "()", 0

        def preprow(r):
            vals = ["%s" % (quote(a),) for a in r]
            return "(%s)" % (", ".join(vals),)

        values = ",\n ".join([preprow(r) for r in rs])
        return values

    def TupleExists(self, tablename, attrs, vals):
        conds = []
        tup = dict(zip(attrs, vals))
        # floating point numbers cannot be safely compared for equality
        keyattrs = self.CachedPrimaryKey(tablename)
        for a in keyattrs:
            conds.append(equals(a, quote(tup[a])))
        # If no primary key, allow any insertion
        if not conds:
            return False
        condition = " AND ".join(conds) or "False"
        sql = '''SELECT exists (SELECT %s FROM %s WHERE %s)''' % (
            "*", tablename, condition)
        s = time.time()
        rs = self.execQuery(sql)
        #logger.debug("Tuple Exists: %s (%s)" % (time.time() - s,rs[0][0]), 6)
        return rs[0][0] == 't'

    def IdempotentInsert(self, tablename, attrs, vals):
        exists = self.TupleExists(tablename, attrs, vals)
        if not exists:
            try:
                self.InsertTuple(tablename, attrs, vals)
            except DuplicateError:
                pass
            return True
        else:
            tup = zip(attrs, vals)
            #logger.debug('IdempotentInsert: Tuple exists. (%s)' % (tup,), 7)
            return False

    def Upsert(self, tablename, keyattrs, keyvals, dataattrs, datavals):
        rs = self.TupleExists(tablename, keyattrs, keyvals)
        if not rs:
            self.InsertTuple(
                tablename,
                keyattrs + dataattrs,
                keyvals + datavals)
        else:
            self.Update(tablename, keyattrs, keyvals, dataattrs, datavals)

    def InsertTuple(self, tablename, attributes, values):
        fqt = fqTableName(tablename)
        ins = '''INSERT INTO %s (%s) VALUES (%s)'''
        attrclause = commas([validIdentifier(a) for a in attributes])
        valclause = commas([quote(x) for x in values])
        sql = ins % (fqt, attrclause, valclause)
        self.execCommand(sql)

    def Update(self, tablename, keyattrs, keyvals, attributes, vals):
        '''Update tuples matching keyattrs=keyvals setting attributes = vals'''

        tbl = fqTableName(tablename)

        setters = zip(attributes, vals)
        setexpr = commas(["%s = %s" %
                          (validIdentifier(k), quote(v)) for k, v in setters])

        where = zip(keyattrs, keyvals)
        condition = " AND ".join(
            [equals(validIdentifier(a), quote(v)) for a, v in where])

        upd = '''UPDATE %s SET %s WHERE %s'''
        self.execCommand(upd % (tbl, setexpr, condition))

    def InsertQry(self, tablename, qry):
        insert = '''INSERT INTO %s (%s)''' % (tablename, qry)
        self.execCommand(insert)

    def createTableAs(self, tablename, qry):
        '''
  Drops and recreates a table based on a query's results.
       '''
        fqt = fqTableName(tablename)
        create = "CREATE TABLE %s AS (%s);" % (fqt, qry)

        if self.checkTable(tablename):
            drop = "DROP TABLE %s; " % fqt
            create = drop + create

        self.execCommand(create)

    def dropTable(self, name):
        '''Drops a table idempotently.'''
        if self.checkTable(name):
            fqt = fqTableName(name)
            drop = 'DROP TABLE %s CASCADE;' % fqt
            self.execCommand(drop)

    def begin(self):
        self.execCommand("BEGIN TRANSACTION")

    def rollback(self):
        self.execCommand("ROLLBACK")

    def commit(self):
        self.execCommand("COMMIT")

    def get_attnames(self, table):
        return self.dbconn.get_attnames(table)

    def checkTable(self, table):
        '''
  returns None if table <name> does not exist
        '''
        sch, name = tablepair(table)
        check = '''
select c.relname
  from pg_class c, pg_namespace n
 where c.relname = '%s'
   and c.relnamespace = n.oid
   and n.nspname = '%s'
''' % (name, sch)
        result = self.execQuery(check)
        return result

    def Union(self, queries):
        return "\n UNION \n".join([q.replace(";", "") for q in queries])

    def createTable(self, name, fields, types, key=None,
                    oids=False, quote=False):
        '''
  idempotent table creator.  Returns False if table already existed.
  fields is a sequence of string column names and
  types is sequence of string type names.

  If quote is True, attribute names are quoted to preserve case sensitivity.
        '''
        exists = self.checkTable(name)
        if quote:
            safe = lambda s: '"%s"' % (s,)
        else:
            safe = validIdentifier
        fields = [safe(f) for f in fields]

        fqt = fqTableName(name)
        if not exists:
            sql = "CREATE TABLE " + fqt + "("
            fieldstypes = ['%s %s' % ft for ft in zip(fields, types)]
            sql = sql + ','.join(fieldstypes)
            if key:
                sql = sql + ', PRIMARY KEY (%s)' % (','.join(key),)
            sql = sql + ')'
            if oids:
                sql = sql + ' WITH OIDS'

            self.execCommand(sql)
            return True
        else:
            return False

    def LoadCSV(
            self,
            fname,
            schema='public',
            overwrite=False,
            quote=False,
            onerror=donothing,
            retry_limit=2):
        '''Load comma-separated values from a file with one line of headers.
  A table with the name of the file will be created in schema <schema>.  All fields will be of type text.i

  If <overwrite> is True, then if the table already exists, it will be deleted.

  If quote is True, column names will be quoted to match the file headers exactly.  Otherwise, illegal identifier names will throw an error.

  onerror is a callback for lines that generate errors.  Signature is handle(data, (e, v, t)),  where data is the line that generated the error, and (e,v,t) is an exception tuple.

  The callback may return a corrected line value that will be retried.

  The number of recursive retries can be controlled with the retry_limit argument.
  '''
        conn = self
        f = file(fname)

        raise Exception('This method requires CMOP python lib')
        #hdrs = cmop.parseCSV(f.readline())

        types = ['text'] * len(hdrs)

        tblname = '%s.%s' % (schema, fname)

        if overwrite:
            conn.dropTable(tblname)
        conn.createTable(tblname, hdrs, types, quote=quote)

        tblname = fqTableName(tblname)
        conn.execCommand("GRANT SELECT ON %s TO public" % (tblname,))

        def process(line, cnt=1):

            #fields = cmop.parseCSV(line)
            try:
                conn.InsertTuple(tblname, hdrs, fields)
            except:
                e, v, t = sys.exc_info()
                newline = onerror(line, (e, v, t))
                if newline and cnt < retry_limit:
                    process(newline, cnt + 1)

        for line in f:
            process(line)

# if __name__ == '__main__':
    # print logger.debug("main")
