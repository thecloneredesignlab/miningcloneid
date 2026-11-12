"""Validate the reporting-only overlay, then enter the verified workflow."""
import sys
import unittest
import test_workflow
import workflow

result=unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromModule(test_workflow))
if not result.wasSuccessful():
    sys.exit(1)
workflow.main()
