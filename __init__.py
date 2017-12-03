"""Initialize python repository for pyPR.

Notes
-----
11/17/17,CM, Initial Commit
"""
import sys, os 
cwd = os.getcwd()
sys.path.append(cwd + '/pyPR')


# %load_ext autoreload # activates the load 
# %autoreload 2 # Reload all the modules


print('Welcome to pyPR!')

pathdata = '/Users/chris/Documents/Research/VLA/VLA2017/'
print('Set the path to your data if needed using pathdata!')
print(pathdata)
