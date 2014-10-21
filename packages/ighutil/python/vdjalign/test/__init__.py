import unittest

modules = [
    'test_ssw',
]

def suite():
    s = unittest.TestSuite()
    for m in modules:
        module = __import__(__package__ + '.' + m, fromlist=m)
        s.addTests(module.suite())

    return s
