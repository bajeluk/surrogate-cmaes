# Log processor for DTS-CMS-ES
#
# Author: Lukas Bajer
# Date: 2019-06-25
#
# (c) Lukas Bajer, 2019

import re
import json
import sys


class CovarianceFunction:
    cov_strings_array = [
        "@covRQiso",
        "@covSEiso",
        "{@covMaterniso, 5}",
        "{@covPoly, 'eye', 1}",
        "{@covSEvlen, {@meanSum, {@meanLinear, @meanConst}}}"
    ]

    def __init__(self):
        self.cov_strings_dict = {self.cov_strings_array[i]: i for i in range(len(self.cov_strings_array))}

    def fun_count(self):
        return len(self.cov_strings_array)

    def to_int(self, name):
        return self.cov_strings_dict[name]

    def to_str(self, i):
        return self.cov_strings_array[i]


class Run:
    def __init__(self):
        self.function = -1
        self.instance = -1
        self.dimension = -1
        self.first_covariances = []
        self.second_covariances = []
        self.evaluations = []
        self.length = 0
        self.cov = CovarianceFunction()

    def setup(self, func, dimension, instance, evaluations, dopt, time):
        self.function = func
        self.dimension = dimension
        self.instance = instance
        self.evaluations_count = evaluations
        self.dopt = dopt
        self.time = time

    def insert_generation(self, first_covariance, second_covariance, evaluations):
        self.first_covariances.append(first_covariance)
        self.second_covariances.append(second_covariance)
        if len(self.evaluations) > 0:
            assert self.evaluations[-1] <= evaluations
        self.evaluations.append(evaluations)
        self.length += 1

    def get_evaluations(self):
        if self.length > 0:
            return self.evaluations[-1]
        else:
            return -1

    def histogram(self, cov_array):
        d = {}
        for c in cov_array:
            if c not in d:
                d[c] = 0
            d[c] += 1
        for i in range(self.cov.fun_count()):
            if i not in d:
                d[i] = 0
        return [d[cov] for cov in range(self.cov.fun_count())]

    def to_json(self):
        res = {
            "function": self.function,
            "dimension": self.dimension,
            "instance": self.instance,
            "evaluationsCount": self.evaluations_count,
            "generationsCount": len(self.first_covariances),
            "distanceOptimum": self.dopt,
            "time": self.time,
            "distinctFirstCovariances": len(set(self.first_covariances).difference([-1])),
            "distinctSecondCovariances": len(set(self.second_covariances).difference([-1])),
            "firstCovariancesHistogram": self.histogram(self.first_covariances),
            "secondCovariancesHistogram": self.histogram(self.second_covariances),
            "firstCovariances": self.first_covariances,
            "secondCovariances": self.second_covariances,
            "evaluations": self.evaluations
        }
        return json.dumps(res, indent=2)


class RunEntryMatcher:
    def __init__(self):
        #   f3 in 20-D, instance 90: FEs=5008 with 0 restarts, fbest-ftarget=5.7707e+01, elapsed time [h]: 105.65
        self.final_line_matcher = re.compile(r"^  f([0-9]+) in ([0-9]+)-D, instance ([0-9]+): "
                                             "FEs=([0-9]+) with.*-ftarget=([-+0-9.eE]+), "
                                             "elapsed time \[h\]: ([-+0-9.eE]+)")
        # Covariance function: @covRQiso
        self.covariance_matcher = re.compile(r"^Covariance function: (.*)")
        # =[DTS]=   99 / 4957(10:0,0) | 6.6e+01 | 6.9e+00 | 0.02 | 0.09 | 0.38 [ *     ] | .  0/371 | 8.89e+00 | 0.00 | 0.00 | 0.05 | 0.00 |
        self.generation_line_matcher = re.compile(r"^=\[DTS\]= *([0-9]+) */ *([0-9]+).*")
        self.cov = CovarianceFunction()

        self.run = Run()
        self.first_covariance = -1
        self.second_covariance = -1
        self.init_run()

    def init_run(self):
        self.run = Run()
        self.first_covariance = -1
        self.second_covariance = -1

    def init_generation(self):
        self.first_covariance = -1
        self.second_covariance = -1

    def match_covariance_line(self, line):
        match = self.covariance_matcher.match(line)
        cov_str = match.group(1)
        cov_int = self.cov.to_int(cov_str)

        if self.first_covariance != -1:
            self.second_covariance = cov_int
        else:
            self.first_covariance = cov_int

    def match_generation_line(self, line):
        match = self.generation_line_matcher.match(line)
        # generation = int(match.group(1))
        evaluations = int(match.group(2))
        # assert self.first_covariance != -1
        # assert self.second_covariance != -1
        assert evaluations > 0
        self.run.insert_generation(self.first_covariance, self.second_covariance, evaluations)

    def setup_from_final_line(self, line):
        #   f3 in 20-D, instance 90: FEs=5008 with 0 restarts, fbest-ftarget=5.7707e+01, elapsed time [h]: 105.65
        match = self.final_line_matcher.match(line)
        self.run.setup(int(match.group(1)), int(match.group(2)), int(match.group(3)),
                       int(match.group(4)), float(match.group(5)), float(match.group(6)))


class LogParser:
    def __init__(self, filename):
        self.final_line_matcher = re.compile(r"^  f([0-9]+) ")
        self.covariance_line_matcher = re.compile(r"^Covariance function: ")
        self.generation_line_matcher = re.compile(r"^=\[DTS\]= *([0-9]+) */")

        self.runs = []
        self.filename = filename

    def process_log(self):
        res = ""
        with open(self.filename, "r") as f:
            matcher = RunEntryMatcher()
            for line in f:
                if self.covariance_line_matcher.match(line):
                    matcher.match_covariance_line(line)
                    # print(line, flush=True)

                if self.generation_line_matcher.match(line):
                    matcher.match_generation_line(line)
                    matcher.init_generation()
                    # print(line, flush=True)

                if self.final_line_matcher.match(line):
                    matcher.setup_from_final_line(line)
                    # print(line, flush=True)
                    res = res + matcher.run.to_json() + "\n"
                    self.runs.append(matcher.run)
                    matcher.init_run()
        return res


if __name__ == '__main__':
    parser = LogParser(sys.argv[1])
    print(parser.process_log())
