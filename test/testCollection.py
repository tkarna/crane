"""
Tests data colletion.

All test cases should be independent, so that they can be executed
without any other component in the sofware bundle. This ensures that
modules remain independent with clean interfaces.

Tuomas Karna 2013-11-05
"""
import unittest
from data.collection import tupleList

class testTupleList(unittest.TestCase) :
  """Tests tupleCollection object"""
  def testCreation(self) :
    tl = tupleList(keywords=['species','color','number_of_legs'])
    self.assertTrue(isinstance(tl,tupleList))
    
  def testAdditionAndGetKeys(self) :
    tl = tupleList(keywords=['species','color','number_of_legs'])
    tl.addSample(('cat','black',4))
    tl.addSample(tup=('cat','red',3))
    tl.addSample(tup=('cat','red',3))
    tl.addSample(species='moose',color='blue',number_of_legs=15)
    tl.addSample(species='booze',color='blue')
    data = tl.getKeys()
    self.assertEqual(data,[('cat','black',4),('cat','red',3),
                           ('moose','blue',15),('booze','blue',None)])
    
  def testAdditionFail(self) :
    tl = tupleList(keywords=['species','color','number_of_legs'])
    with self.assertRaises(Exception) as e :
      tl.addSample(('cat','black'))
    self.assertEqual(e.exception.args[0],'Given tuple has incorrect length: 2 != 3')

  def testFiltering(self):
    tl = tupleList(keywords=['species','color','number_of_legs'])
    tl.addSample(('cat','blue',4))
    tl.addSample(tup=('cat','red',3))
    tl.addSample(tup=('cat','red',3))
    tl.addSample(species='moose',color='blue',number_of_legs=15)
    tl.addSample(species='booze',color='blue')
    data = tl.getTuples(species='cat')
    self.assertEqual(data,[('cat','blue',4),('cat','red',3)])
    data = tl.getTuples(color='blue')
    self.assertEqual(data,[('cat','blue',4),('moose','blue',15),('booze','blue',None)])
    data = tl.getTuples(color='blue',species='moose')
    self.assertEqual(data,[('moose','blue',15)])
    data = tl.getTuples(color='blue',species=['cat','moose'])
    self.assertEqual(data,[('cat','blue',4),('moose','blue',15)])
    data = tl.getTuples(color='red',exclude=True)
    self.assertEqual(data,[('cat','blue',4),('moose','blue',15),('booze','blue',None)])
    data = tl.getTuples(species=['moose','booze'],color='blue',exclude=True)
    self.assertEqual(data,[('cat','blue',4),('cat','red',3)])

if __name__ == '__main__':
  """Run all tests"""
  unittest.main()
