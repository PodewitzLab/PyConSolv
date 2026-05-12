import unittest
from unittest import mock

import tests  # noqa: F401
from PyConSolv.interfaces.mdengines import MDEngine


class TestMDEngineDispatch(unittest.TestCase):
    def test_amber_dispatch(self):
        with mock.patch('PyConSolv.interfaces.mdengines.amberInterface') as amber:
            with mock.patch('os.chdir'):
                MDEngine('/tmp', engine='amber')
            amber.assert_called_once_with('/tmp')

    def test_gromacs_dispatch(self):
        with mock.patch('PyConSolv.interfaces.mdengines.gromacsInterface') as grom:
            MDEngine('/tmp', engine='gromacs')
            grom.assert_called_once_with('/tmp')

    def test_run_calls_prepare_and_equilibrate(self):
        with mock.patch('PyConSolv.interfaces.mdengines.amberInterface') as amber:
            instance = amber.return_value
            instance.equilibrate.return_value = 1
            eng = MDEngine('/tmp', engine='amber')
            status = eng.run('/tmp', cpus=4)
            self.assertEqual(status, 1)
            instance.prepare.assert_called_once()
            instance.equilibrate.assert_called_once_with(4)


if __name__ == '__main__':
    unittest.main()
