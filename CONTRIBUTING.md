# Reproduction Process

To repoduce the issue, I first commented out the call to get() in line 1487 and replaced it with a call to flag instead like so:

chemistry = self.flag("chemistry")

The terminal blanks out seen in the following output:
<summary><b>chemistry</b> (`chemistry`)</summary>

I expect the computer to know what to do when it is given a new function like flag() and output the correct result.

# Solution Approach

1. At the moment, case_validator.py reads boolean case flags with the verbose form self.get("name", "F") == "T" in 134 places. To clean this code up I want to implement a map function. However, I am prevented from doing so at the moment because of the hard coding in ast_analyzer.py in the _build_local_param_map() method in line 232: call.func.attr == "get". Thus, when the analyzer sees the word "flag", it doesn't recognize it and blanks out.

2. My solution is to first make sure the analyzer recognizes flag() by adding a conditional statement in _build_local_param_map() in ast_analyzer.py. I will then implement a flag() method in case_validator.py that maps the verbose form to a simple boolean value.

3. Files I will touch include:
   - case_validator.py
   - ast_analyzer.py
  
4. I will verify my approach works by calling flag() instead of get() in line 1487 of case_validator.py. If the analyzer correctly returns a value instead of blank, it successfully recognized the flag() function. I will also test the flag() method itself to ensure it correctly returns true or false.
