# Reproduction Process

Observations:
   - To repoduce the issue, I first commented out the call to get() in line 1487 and replaced it with a call to flag instead.
   - The terminal blanks and doesn't output anything when the flag() method is called.
   - I expect the computer to know what to do when it is given a function like flag() and output the correct result.

Environment Setup: I did not face any problems during the environment setup process as I have been using.

Steps to Reproduce:
   1. Open the file case_validator.py using the file path toolchain/mfc/case_validator.py
   2. Go to line 1487 in the check_chemistry() method
   3. Replace the value assigned to the chemistry variable with self.flag("chemistry")
   4. Reproduce the issue on your terminal by running the following: python3 toolchain/mfc/gen_case_constraints_docs.py | grep "ANALYZER"

Branch Link: https://github.com/mansibrahman03/MFlowCode/tree/main

# Solution Approach

Implementation:

1. At the moment, case_validator.py reads boolean case flags with the verbose form self.get("name", "F") == "T" in 134 places. To clean this code up I want to implement a map function. However, I am prevented from doing so at the moment because of the hard coding in ast_analyzer.py in the _build_local_param_map() method in line 232: call.func.attr == "get". Thus, when the analyzer sees the word "flag", it doesn't recognize it and blanks out.

2. My solution is to first make sure the analyzer recognizes flag() by adding a conditional statement in _build_local_param_map() in ast_analyzer.py. I will then implement a flag() method in case_validator.py that maps the verbose form to a simple boolean value.

3. Files I will touch include:
   - case_validator.py
   - ast_analyzer.py
  
4. I will verify my approach works by calling flag() instead of get() in line 1487 of case_validator.py. If the analyzer correctly returns a value instead of blank, it successfully recognized the flag() function. I will also test the flag() method itself to ensure it correctly returns true or false.

# Implementation Notes

Implementation Progress: I fixed _build_local_param_map() in ast_analyzer.py to recognize flag() method and added a flag() helper method in case_validator.py.

Challenges Faced: The main challenge was understanding the codebase.

Testing Strategy: I tested the flag method itself to ensure it correctly return a bool for all possible params (T, F, null). I also test to see that the old broken code fails and the new code works. 

Branch Link: https://github.com/mansibrahman03/MFlowCode/tree/main

# PR

PR Link: [Direct link to your submitted pull request](https://github.com/MFlowCode/MFC/pull/1622)

PR Description: I added a flag() method in toolchain/mfc/case_validator.py that functions to replace the verbose idiom "self.get(x,F)==T" repeated in 134 places across the file. I also updated ast_analyzer.py to recognize calls made to flag().

Maintainer Feedback: I am currently waiting for maintainer to respond.

Status: Awaiting review
