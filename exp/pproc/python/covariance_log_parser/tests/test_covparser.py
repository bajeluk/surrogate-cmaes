from unittest import TestCase
import os
import tempfile
import json

from covparser.covparser import *


class TestCovarianceFunction(TestCase):
    def setUp(self):
        self.cf = CovarianceFunction()

    def test_to_int(self):
        assert self.cf.to_int("@covRQiso") == 0
        assert self.cf.to_int("{@covSEvlen, {@meanSum, {@meanLinear, @meanConst}}}") == 4

    def test_to_str(self):
        assert self.cf.to_str(0) == "@covRQiso"
        assert self.cf.to_str(4) == "{@covSEvlen, {@meanSum, {@meanLinear, @meanConst}}}"


class TestRun(TestCase):
    def setUp(self):
        self.run = Run()

    def test_empty(self):
        assert self.run.first_covariances == []
        assert self.run.length == 0

    def test_initialize(self):
        self.run.setup(1, 2, 3, 10, 1.5e6, 0.23)
        assert self.run.first_covariances == []
        assert self.run.function == 1
        assert self.run.dimension == 2
        assert self.run.instance == 3
        assert self.run.evaluations_count == 10
        assert self.run.dopt == 1.5e6
        assert self.run.time == 0.23

    def test_insert(self):
        self.run.insert_generation(1, 2, 3)
        self.run.insert_generation(1, 2, 4)
        assert self.run.evaluations[-1] == 4
        assert self.run.length == 2


class TestRunEntryMatcher(TestCase):
    def setUp(self):
        self.cov_line1 = "Covariance function: {@covMaterniso, 5}"
        self.cov_line2 = "Covariance function: @covRQiso"
        self.gen_line = "=[DTS]=   99 / 4957(10:0,0) | 6.6e+01 | 6.9e+00 | 0.02 | 0.09 | 0.38 [ *     ] | .  0/371 | 8.89e+00 | 0.00 | 0.00 | 0.05 | 0.00 |"
        self.final_line = "  f3 in 20-D, instance 90: FEs=5008 with 0 restarts, fbest-ftarget=5.7707e+01, elapsed time [h]: 105.65"
        self.matcher = RunEntryMatcher()
        self.cov = CovarianceFunction()

    def test_match_covariance_line(self):
        self.matcher.match_covariance_line(self.cov_line1)
        self.matcher.match_covariance_line(self.cov_line2)

        assert self.matcher.first_covariance == self.cov.to_int("{@covMaterniso, 5}")
        assert self.matcher.second_covariance == self.cov.to_int("@covRQiso")

    def test_match_generation_line(self):
        self.matcher.match_covariance_line(self.cov_line1)
        self.matcher.match_covariance_line(self.cov_line2)
        self.matcher.match_generation_line(self.gen_line)

        assert self.matcher.run.length == 1
        assert self.matcher.run.get_evaluations() == 4957
        assert self.matcher.run.first_covariances[-1] == self.cov.to_int("{@covMaterniso, 5}")
        assert self.matcher.run.second_covariances[-1] == self.cov.to_int("@covRQiso")

    def test_setup_from_final_line(self):
        self.matcher.setup_from_final_line(self.final_line)

        assert self.matcher.run.function == 3
        assert self.matcher.run.dimension == 20
        assert self.matcher.run.instance == 90
        assert self.matcher.run.evaluations_count == 5008
        assert self.matcher.run.dopt == 5.7707e+01
        assert self.matcher.run.time == 105.65


class TestLogParser(TestCase):
    def setUp(self):
        stdout = \
"""
=[DTS]=   12 /   28( 1:0,0) | 2.4e+01 | 6.4e-07 | 0.00 | 0.00 | 0.00 [ ***** ] | .  0/ 29 | 1.04e+01 | 0.00 | 0.00 | 0.05 | 0.00 |
=[DTS]=   13 /   30( 1:0,0) | 2.4e+01 | 6.4e-07 | 0.00 | 0.00 | 0.00 [ ***** ] | .  0/ 29 | 1.04e+01 | 0.00 | 0.00 | 0.05 | 0.00 |
Covariance function: {@covPoly, 'eye', 1}
Covariance function: @covRQiso
=[DTS]=   14 /   31( 1:0,0) | 1.7e+01 | 4.9e-06 | 0.00 | 0.00 | 0.00 [ ***** ] | .  0/ 30 | 8.16e+00 | 0.00 | 0.00 | 0.05 | 0.00 |
Covariance function: {@covMaterniso, 5}
Covariance function: {@covPoly, 'eye', 1}
=[DTS]=   15 /   32( 1:0,0) | 0.0e+00 | 1.3e-05 | 0.00 | 0.00 | 0.00 [ ***** ] | .  0/ 31 | 9.34e+00 | 0.00 | 0.00 | 0.05 | 0.00 |
  f5 in 10-D, instance 90: FEs=33 with 0 restarts, fbest-ftarget=1.0000e-08, elapsed time [h]: 0.36
"""

        self.tmpfile = os.path.join(tempfile.gettempdir(), "tmp.o123")
        # self.tmpfile = tmpfilepath + "/tmp.o123"
        with open(self.tmpfile, 'w') as f:
            f.write(stdout)

        self.log_parser = LogParser(self.tmpfile)

    def test_process_log(self):
        gson = self.log_parser.process_log()
        struct = json.loads(gson)

        assert struct['function'] == 5
        assert struct['dimension'] == 10
        assert struct['instance'] == 90
        assert struct['evaluationsCount'] == 33
        assert struct['generationsCount'] == 4
        assert struct['distanceOptimum'] == 1e-8
        assert struct['time'] == 0.36
        assert struct['distinctFirstCovariances'] == 2
        assert struct['distinctSecondCovariances'] == 2
        assert struct["firstCovariancesHistogram"] == [0, 0, 1, 1, 0]
        assert struct["secondCovariancesHistogram"] == [1, 0, 0, 1, 0]

        assert struct['firstCovariances'] == [-1, -1, 3, 2]
        assert struct['secondCovariances'] == [-1, -1, 0, 3]
        assert struct['evaluations'] == [28, 30, 31, 32]

    def tearDown(self):
        os.remove(self.tmpfile)
