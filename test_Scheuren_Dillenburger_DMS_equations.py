import Scheuren_Dillenburger_DMS_equations as dms
import unittest

class DMSRK4TestCase(unittest.TestCase):

    def setUp(self) -> None:
        return super().setUp()

    
    def tearDown(self) -> None:
        return super().tearDown()


    def test_f_example_initial_conditions(self):
        initial_conditions = dms.f_example(0, 2)
        self.assertEqual(initial_conditions, 2)


    def test_f_example_rk4(self):
        pass


if __name__ == '__main__':
    unittest.main()