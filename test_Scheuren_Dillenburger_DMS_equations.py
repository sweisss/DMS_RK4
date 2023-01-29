import Scheuren_Dillenburger_DMS_equations as dms
import unittest


class DMS_RK4TestCase(unittest.TestCase):

    def setUp(self) -> None:
        return super().setUp()

    
    def tearDown(self) -> None:
        return super().tearDown()


    def test_f_example_initial_conditions(self):
        initial_conditions = dms.f_example(0, 2)
        self.assertEqual(initial_conditions, 2)


    def test_f_example_rk4(self):
        t0 = 0
        x0 = 2
        h = 0.1
        steps = 20
        data, headers = dms.rk4(dms.f_example, t0, x0, h, steps)
        expected_data =  {
            0: 2,
            0.1: 2.2051708333333333,
            0.2: 2.4214025708506943,
            0.3: 2.6498584970625374, 
            0.4: 2.8918242400806853, 
            0.5: 3.1487206385968376, 
            0.6: 3.4221179620919324, 
            0.7: 3.713751626596776, 
            0.8: 4.025539563292314, 
            0.9: 4.359601413780069, 
            1.0: 4.718279744135164, 
            1.1: 5.10416349005898, 
            1.2: 5.520113867778058, 
            1.3: 5.969293010013833, 
            1.4: 6.455195613621163, 
            1.5: 6.9816839156353785, 
            1.6: 7.553026347779348, 
            1.7: 8.173940256297257, 
            1.8: 8.84963911466892, 
            1.9: 9.585884701724579, 
            2.0: 10.389044767375538
        }
        expected_headers = ['tn', 'xn']
        self.assertEqual(data, expected_data)
        self.assertEqual(headers, expected_headers)

    
    def test_eqn_1(self):
        c_i_0 = 500    # From Scheuren's Excel 
        k = 0.00130809768004509     # From Scheuren's Excel
        c_i = dms.eqn_1(c_i_0, k, 60, 5)
        expected = [
            500.0,
            337.7070468652124,
            228.09209900484552,
            154.05661833642804,
            104.05201125687249,
            70.27819488388847,
            47.466883306511725,
            32.059801970675444,
            21.65364209320064,
            14.625175050342087,
            9.878049352275617,
            6.671773751091645,
            4.506210021668004
        ]
        self.assertEqual(c_i, expected)


    # def test_eqn_4_exists(self):
    #     if hasattr(dms, 'eqn_4'):
    #         return lambda func: func
    #     return unittest.skip("{!r} doesn't have {!r}".format(dms, 'eqn_4'))


if __name__ == '__main__':
    unittest.main()