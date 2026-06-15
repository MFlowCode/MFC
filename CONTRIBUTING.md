Part 1: Reproducing the Issue

To repoduce the issue, I first commented out the call to get() in line 1487 and replaced it with a call to flag instead like so:

# chemistry = self.get("chemistry", "F") == "T"
chemistry = self.flag("chemistry")

The terminal blanks out seen in the following output:
<summary><b>chemistry</b> (`chemistry`)</summary>

I expect the computer to know what to do when it is given a new function like flag() and output the correct result.

Part 2: Solution Plan

1. The root cause is the hard coding in ast_analyzer.py in the _build_local_param_map() function in line 232: call.func.attr == "get". Thus, when the analyzer sees the word "flag", it doesn't recognize it and blanks out.

2. My solution is to implement a flag() method that maps calls the get() method